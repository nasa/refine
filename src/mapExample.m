format long

p7 = 1.0/sqrt(2);

eVect1 = [  p7
            p7 ]
eVect2 = [ -p7
            p7 ]

eigenVector = [ eVect1 eVect2 ]

h1 = 2.0*sqrt(2)
h2 = sqrt(2)
eVal1 = 1/h1/h1
eVal2 = 1/h2/h2
sqrteVal1 = 1/h1
sqrteVal2 = 1/h2

eigenValue  = [ eVal1   0.0
                  0.0 eVal2 ]

sqrteigenValue  = [ sqrteVal1   0.0
                     0.0     sqrteVal2 ]

map = sqrteigenValue * eigenVector'

Mmatrix = map'*map
Mmatrix = eigenVector*eigenValue*eigenVector'

x1 = [ 2
       2 ]
e1=map*x1
len1=norm(e1)
len1=x1'*Mmatrix*x1
x2 = [ -1
        1 ]

e2=map*x2
len2=norm(e2)
len2=x2'*Mmatrix*x2
