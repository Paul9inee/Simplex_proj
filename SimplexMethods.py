import sys

nearzero = sys.float_info.epsilon;

def CopyTableau( tableau ):
    dest = {};
    for ks in tableau.keys():
        dest[ks] = tableau[ks];

    return dest;

def IndepVarsIndex( tableau ):
    tableau['x'] = [];
    for idx in tableau.keys():
        xnum = idx.strip('x');
        if xnum.isdigit():
            j = int( xnum );
            tableau['x'].append(j);

    tableau['x'] = sorted( tableau['x'] );
    # Fix size
    rowmax = len( tableau['t'] );
    colmax = len( tableau['x'] );
    sv = str(rowmax)+'x'+str(colmax);    
    tableau['size'] = sv;

def DepVarsIndex( tableau ):
    tableau['t'] = [];
    for idx in tableau.keys():
        tnum = idx.strip('t');
        if tnum.isdigit():
            i = int( tnum );
            tableau['t'].append(i);

    tableau['t'] = sorted( tableau['t'] );
    # Fix size
    rowmax = len( tableau['t'] );
    colmax = len( tableau['x'] );
    sv = str(rowmax)+'x'+str(colmax);    
    tableau['size'] = sv;

def GetValue( tableau, i, j ):
    sij = str(i)+','+str(j);
    if sij in tableau.keys(): return tableau[sij];

    return 0.0;

def GetValueB( tableau, irow ):
    Bi = 'B'+str(irow);
    if Bi in tableau.keys(): return tableau[Bi];

    return 0.0;

def GetValueC( tableau, jcol ):
    Cj = 'C'+str(jcol);
    if Cj in tableau.keys(): return tableau[Cj];

    return 0.0;

def GetValueD( tableau ):
    if 'd' in tableau.keys(): return tableau['d'];

    return 0.0;

def SetValue(tableau, i, j, value):
    sij = str(i)+','+str(j);
    if abs(value) > 4*nearzero:	
        tableau[sij] = value; 
        return;

    if sij in tableau.keys(): del tableau[sij];

def SetValueB(tableau, irow, value):
    Bi = 'B'+str(irow);
    if abs(value) > 4*nearzero:	
        tableau[Bi] = value; 
        return;

    if Bi in tableau.keys(): del tableau[Bi];

def SetValueC(tableau, jcol, value):
    Cj = 'C'+str(jcol);
    if abs(value) > 4*nearzero:	
        tableau[Cj] = value; 
        return;

    if Cj in tableau.keys(): del tableau[Cj];

def SetValueD(tableau, value):
    if abs(value) > 4*nearzero:	
        tableau['d'] = value; 
        return;

    if 'd' in tableau.keys(): del tableau['d'];

def DelRow( tableau, drow ):
    for j in tableau['x']:
        aij = str(drow)+','+str(j);
        tableau.pop(aij,None); 

    Bi = 'B'+str(drow);
    tableau.pop(Bi, None);

    ti = 't'+str(drow);
    tableau.pop(ti, None);

    # fix row entry numbers 
    DepVarsIndex( tableau );

def DelColumn( tableau, dcol ):
    for i in tableau['t']:
        aij = str(i)+','+str(dcol);
        tableau.pop(aij,None); 

    Cj = 'C'+str(dcol);
    tableau.pop(Cj,None);

    xi = 'x'+str(dcol);
    tableau.pop(xi,None);

    # fix column entry numbers 
    IndepVarsIndex( tableau );

# tableau info 

def GetSize(tableau):
    sz = tableau['size'];
    pair = sz.split('x');
    svec = int(pair[0]), int(pair[1]);

    return svec;

def EqualityIndex( tableau ):
    eqlst = [];
    for tvar in tableau.keys():
        if tvar[0] == 't' and tableau[tvar] == '0':
            tnum = tvar.strip('t');
            eqnum = int( tnum );
            eqlst.append(eqnum);

    return eqlst;

def NonCanonicalIndex(tableau):
    nonlst = [];
    for xvar in tableau.keys():
        if xvar[0] == 'x' and tableau[xvar][0] == '0':
            xnum = xvar.strip('x');
            nonnum = int( xnum );
            nonlst.append(nonnum);

    return nonlst;

def AugEquality(tableau, erow):
    est = 't'+str(erow);
    mxrow = max( tableau['t'] );
    # make one row at the end
    for j in tableau['x']:
        vmrj = -GetValue(tableau, erow, j);
        SetValue(tableau, mxrow+1, j, vmrj);
    vb = -GetValueB(tableau, erow);
    SetValueB(tableau, mxrow+1, vb); 

    #  attach +/- at the equality variables
    tableau[est] = est+'+';
    endrownext = 't'+str(mxrow+1);
    tableau[endrownext] = est+'-';
    # re-calculate the size
    SetSizeOfTableau( tableau );
    print('Fix Equality',tableau);

def AugNonCanonical(tableau, nccol):
    ncvar = 'x'+str(nccol);
    mxcol = max( tableau['x'] );
    # make one column at the end
    for i in tableau['t']:
        vinc = -GetValue(tableau, i, nccol);
        SetValue(tableau, i, mxcol+1, vinc);
    vc = -GetValueC(tableau, nccol);
    SetValueC(tableau, mxcol+1, vc); 
    #  attach +/- at the noncanonical variables
    tableau[ncvar] = ncvar+'+';
    xnc = 'x'+str(mxcol+1);
    tableau[xnc] = ncvar+'-';
    # re-calculate the size
    SetSizeOfTableau( tableau );
    print('Fix Equality',tableau);

def SetSizeOfTableau(tableau):
    # sets we need here
    rowset = set();
    colset = set();

    for wd in tableau.keys():
        pair = wd.split(',');
        if pair[0].isdecimal() and pair[1].isdecimal(): 
            irow = int(pair[0]);
            jcol = int(pair[1]);
            rowset.add(irow);
            colset.add(jcol);

    rowmax = len( rowset );
    colmax = len( colset );

    tableau['t'] = sorted( rowset );
    tableau['x'] = sorted( colset );

    sv = str(rowmax)+'x'+str(colmax);    
    tableau['size'] = sv; # size increases

    return rowmax, colmax;

def Augmenting(tableau):
    # equality contraints
    eqidx = EqualityIndex( tableau );
    print('Equality = ', eqidx);
    for i in eqidx: AugEquality(tableau, i);
        
    # noncanonical variables
    nonidx = NonCanonicalIndex( tableau );
    print('Non Canonical Variables = ', nonidx);
    for j in nonidx: AugNonCanonical(tableau, j);

def SaveNonCanonical( tableau, irow ):
        ti = 't'+str(irow);
        if tableau[ti][0] == '0':
            Oti = '0'+ti;
            tableau[Oti] = {}; 
            for j in tableau['x']:
                var = tableau['x'+str(j)];
                tableau[Oti][var] = GetValue( tableau, irow, j );
            tableau[Oti]['B'+str(irow)] = GetValueB( tableau, irow );
            return True;
        else: return False;

def Reducing(tableau):
    # reducing equality contraints
    eqidx = EqualityIndex( tableau );
    print('Equality = ', eqidx);
    for i in eqidx: 
        jp = PivotOnEquality(tableau, i);
        if jp == -1:
            Bi = 'B'+str(i);
            if abs(tableau[Bi]) < 4*nearzero:
                DelRow( tableau, i );
                continue;
            else: return i, jp, 'infeasible';
        else:
            DelColumn( tableau, jp ); # remove one column
        # save non-canonical variable
        if SaveNonCanonical( tableau, i ): DelRow( tableau, i );

    # reducing noncanonical variables
    nonidx = NonCanonicalIndex( tableau );
    print('Non Canonical Variables = ', nonidx);
    for j in nonidx: 
        PivotOnNonCanonical(tableau, j); #HERE

def VarsOnTableau(tableau):
    # initial slack variables(t-vector)
    skey = tableau.keys();

    for i in tableau['t']:
        st = 't'+str(i);
        if not st in skey: tableau[st] = st;

    # initial primitive variables(x-vector)
    for j in tableau['x']: 
        sx = 'x'+str(j); 
        if not sx in skey: tableau[sx] = sx;

    tableau['f'] = -GetValueD(tableau);
    tableau['-1'] = -1;

# Reducing 

def PivotOnEquality( tableau, eqrow ):
    rowvec = [ abs(GetValue( tableau, eqrow, j )) for j in tableau['x'] ]; 
    jp = rowvec.index( max(rowvec) );

    if rowvec[jp] < 4*nearzero: return -1;

    pv = eqrow, jp, GetValue( tableau, eqrow, jp );
    PivotTransform( tableau, pv );

    return jp;

# utilities

def PrintBasicSolution( tableau ):
    GetBasicSolution( tableau );
    bskey = sorted(tableau['basic'].keys()); 
    s = [ (sk, tableau['basic'][sk]) for sk in bskey ];
    print( 'Basic solution = ',s);
    print( tableau['objective'],' = ', tableau['f'] );

def PrintLine( n1, c1, n2, c2 ):
    print('+', end='');
    for j in range(n1): print(c1, end=''); 
    print('+', end='');

    for j in range(n2): print(c2, end=''); 
    print('+');

def PrintTableau(comment, tableau):
    vs = GetSize(tableau);
    mr = vs[0]; nc = vs[1];

    print(comment);
    PrintLine( nc*12+1, '-', 12, '-' );

    for i in tableau['t']:
        print('|', end='');
        for j in tableau['x']:
            print("%12.4e" %GetValue(tableau,i,j),end='');
        print(" |%12.4e|" %GetValueB(tableau,i))

    PrintLine( nc*12+1, '-', 12, '-' );

    print('|', end='');
    for j in tableau['x']:
        print("%12.4e" %GetValueC(tableau,j),end='');
    print(" |%12.4e|" %GetValueD(tableau));

    PrintLine( nc*12+1, '-', 12, '-' );

# Pivot transformation

def PivotTransform(tableau, pivot):
    if pivot[0] < 0 or pivot[1] < 0:
        print('Solution = ',pivot[2]); return;

    ip = pivot[0];
    jp = pivot[1];
    pval = pivot[2];

    tmptbl = {}; # create temporary tableau
    tmptbl['size'] = tableau['size'];

    # itself 
    SetValue(tmptbl, ip, jp, 1.0/pval);

    # row transform
    for j in tableau['x']:
        if j != jp:
            vj = GetValue(tableau, ip, j);
            SetValue(tmptbl, ip, j, vj/pval);
    vip = GetValueB(tableau, ip);
    SetValueB(tmptbl, ip, vip/pval);

    # column transform
    for i in tableau['t']:
        if i != ip:
            vi = GetValue(tableau, i, jp);
            SetValue(tmptbl, i, jp, -vi/pval);
    vjp = GetValueC(tableau, jp);
    SetValueC(tmptbl, jp, -vjp/pval);

    # other parts transform
    for i in tableau['t']:
        if i != ip:
            bi = GetValueB(tableau, i);
            bip = GetValueB(tableau, ip);
            vijp = GetValue(tableau, i, jp);
            vtrns = bi - bip*vijp/pval;
            SetValueB(tmptbl, i, vtrns);
        for j in tableau['x']:
            if i != ip and j != jp:
                vij = GetValue(tableau, i, j);
                vi = GetValue(tableau, ip, j);
                vj = GetValue(tableau, i, jp);
                vtrns = vij - vi*vj/pval;
                SetValue(tmptbl, i, j, vtrns);

    for j in tableau['x']:
        if j != jp:
            cj = GetValueC(tableau, j);
            cjp = GetValueC(tableau, jp);
            vipj = GetValue(tableau, ip, j);
            vtrns = cj - cjp*vipj/pval;
            SetValueC(tmptbl, j, vtrns);

    dv = GetValueD(tableau);
    bip = GetValueB(tableau, ip);
    cjp = GetValueC(tableau, jp);
    vtrns = dv - cjp*bip/pval;
    SetValueD(tmptbl, vtrns);

    # exchange variables
    xv = 'x'+str(jp);
    tv = 't'+str(ip);
    tableau[xv], tableau[tv] = tableau[tv], tableau[xv];

    # reset tableau 
    for i in tableau['t']:
        Bi = 'B'+str(i);
        if Bi in tableau.keys(): del tableau[Bi];
        for j in tableau['x']:
            sk = str(i)+','+str(j);
            if sk in tableau.keys(): del tableau[sk];

    for j in tableau['x']:
        Cj = 'C'+str(j);
        if Cj in tableau.keys(): del tableau[Cj];

    if 'd' in tableau.keys(): del tableau['d'];

    # copy to the original tableau 
    for Aij in tmptbl.keys(): tableau[Aij] = tmptbl[Aij];

    tableau['f'] = tableau['-1']*GetValueD( tableau );

# make tableau

def AssembleAugLP(tableau, b, c, d):
    lenB = len(b);
    lenC = len(c);
    # attach B at the right
    for i in tableau['t']:
        if i < lenB: SetValueB(tableau, i, b[i]);
    # attach C at the bottom
    for j in tableau['x']:
        if j < lenC: SetValueC(tableau, j, c[j]);
    # attach d at the bottom right corner
    SetValueD(tableau, d);

    # Augmented equalities and non-canonical variables 
    Augmenting( tableau );

    # Set variables on the tableau
    VarsOnTableau( tableau );

    # Set indecies of variables 
    IndepVarsIndex( tableau );
    DepVarsIndex( tableau );

def AssembleRedLP(tableau, b, c, d):
    lenB = len(b);
    lenC = len(c);
    # attach B at the right
    for i in tableau['t']:
        if i < lenB: SetValueB(tableau, i, b[i]);
    # attach C at the bottom
    for j in tableau['x']:
        if j < lenC: SetValueC(tableau, j, c[j]);
    # attach d at the bottom right corner
    SetValueD(tableau, d);

    # Augmented equalities and non-canonical variables 
    Reducing( tableau );

    # Set variables on the tableau
    VarsOnTableau( tableau );

    # Set indecies of variables 
    IndepVarsIndex( tableau );
    DepVarsIndex( tableau );

# pivot selection

def NegativeOnB(tableau):
    rev_t = tableau['t'][::-1]; 
    for i in rev_t:
        if GetValueB(tableau, i) < 0.0: return i;

    return -1;

def PositiveOnC(tableau):
    rev_x = tableau['x'][::-1]; 
    for j in rev_x:
        if GetValueC(tableau, j) > 0.0: return j;

    return -1;

def GetNegativeCol(tableau, nrow):
    for j in tableau['x']:
        prepv = GetValue(tableau, nrow, j);
        if prepv < 0.0: return nrow, j, prepv;

    return nrow, -1,'infeasible';

def GetPositiveRow(tableau, ncol):
    for i in tableau['t']:
        prepv = GetValue(tableau, i, ncol);
        if prepv > 0.0: return i, ncol, prepv;

    return -1, ncol,'unbounded';
	
def PivotOnB(tableau):
    negb = NegativeOnB(tableau);
    if negb == -1: return -1, -1,'feasible';

    prepivot = GetNegativeCol(tableau,negb);
    if prepivot[1] == -1: return prepivot;

    prow = prepivot[0]; pcol = prepivot[1]; prevalue = prepivot[2];
    small = GetValueB(tableau, prow)/prevalue;
    starting = tableau['t'].index(negb);
    ending = len(tableau['t']);
    for k in range(starting+1,ending):
        i = tableau['t'][k];
        pval = GetValue(tableau, i, pcol);
        if pval > 0.0: 
            ratio = GetValueB(tableau, i)/pval;
            if ratio < small: small = ratio; prow = i;

    return prow, pcol, GetValue(tableau, prow, pcol);

def PivotOnC(tableau):
    posc = PositiveOnC(tableau);
    if posc == -1: return -1, -1,'optimal';

    pv = GetPositiveRow(tableau,posc);
    if pv[0] == -1: return pv;

    prow = pv[0]; pcol = pv[1]; pval = pv[2];
    small = GetValueB(tableau, prow)/pval;
    starting = tableau['t'].index(prow);
    ending = len(tableau['t']);
    for k in range(starting+1,ending):
        i = tableau['t'][k];
        pval = GetValue(tableau, i, pcol);
        if pval > 0.0: 
            ratio = GetValueB(tableau, i)/pval;
            if ratio < small: small = ratio; prow = i;

    return prow, pcol, GetValue(tableau, prow, pcol);

def FeasibleB( tableau ):
    pvB = PivotOnB( tableau ); print('Pivot = ',pvB);
    while pvB[0] >= 0 and pvB[1] >= 0: 
        PivotTransform( tableau , pvB );
        pvB = PivotOnB( tableau );
        print('Pivot = ',pvB);

    return pvB;

def UnboundedC( tableau ):
    pvC = PivotOnC( tableau ); print('Pivot = ',pvC);
    while pvC[0] >= 0 and pvC[1] >= 0: 
        PivotTransform( tableau , pvC );
        pvC = PivotOnC( tableau );
        print('Pivot = ',pvC);

    return pvC;

def GetBasicSolution(tableau):
    tableau['basic'] = {};
    # Independent Basic solution
    for j in tableau['x']:
        sj = 'x'+str(j);
        if tableau[sj][0] == 'x': 
            tableau['basic'][tableau[sj]] = 0.0;

    # Dependent Basic solution
    for i in tableau['t']:
        si = 't'+str(i);
        if tableau[si][0] == 'x': 
            tableau['basic'][tableau[si]] = GetValueB(tableau, i);

# Simplex for Canonical Maximum Linear Programming

def SimplexMaxLP( tableau, b, c, d, method='augmented' ):
    SetSizeOfTableau( tableau ); # Fix the size of the tableau

    # Complete tableau by attaching B, C, and d
    if method == 'augmented':
        AssembleAugLP( tableau , b, c, d );
    else:
        AssembleRedLP( tableau , b, c, d );

    PrintTableau('Canonical Maximum Tableau', tableau ) 
    original = CopyTableau( tableau );

    Bstatus = FeasibleB( tableau );
    if Bstatus[2] == 'feasible': # Make B >= 0 
        result = UnboundedC( tableau ); # Make C <= 0
        if result[2] == 'unbounded': 
            PrintTableau('Solution-> *Feasible unbounded', tableau ); 
        else: 
            PrintTableau('Optimal solution->', tableau );
            PrintBasicSolution( tableau );
    elif Bstatus[2] == 'infeasible': 
        PrintTableau('Solution-> *Infeasible', tableau );

    return original;

#end

if __name__ == '__main__':
    print(__name__);
