# -*- coding: utf-8 -*-

import sys
import SimplexMethods as sms 

#from Samples.Problem_0 import Aij, b, c, d
#from Samples.Problem_1 import Aij, b, c, d
#from Samples.Problem_2 import Aij, b, c, d
#from Samples.Problem_3 import Aij, b, c, d
#from Samples.Problem_4 import Aij, b, c, d
#from Samples.Problem_5 import Aij, b, c, d
#from Samples.Problem_6 import Aij, b, c, d
from Samples.Problem_7 import Aij, b, c, d

# Solve the Canonical Maximum LP

print('Aij = ', Aij);
print('B = ', b);
print('C = ', c);
tableau0 = sms.SimplexMaxLP( Aij, b, c, d );

#sms.PrintTableau('Original-Tableau ', tableau0 );
#sms.DelRow( tableau0, 1 );
#print('Reduced(row) = ', tableau0 );
#sms.PrintTableau( 'Reduced-Tableau(1-row)', tableau0 );
#sms.DelColumn( tableau0, 2);
#print('Recuced(col) = ', tableau0)
#sms.PrintTableau( 'Reduced-Tableau(2-col)', tableau0 );

#end
