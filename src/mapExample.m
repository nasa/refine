format long

p7 = 1.0/sqrt(2);
p7 = 0.5;

eVect1 = [  p7
            p7 ]
eVect2 = [ -p7
            p7 ]

eigenVector = [ eVect1 eVect2 ]

eVal1 = 0.5
eVal2 = 1.0

eigenValue  = [ eVal1   0.0
                  0.0 eVal2 ]

map = eigenValue * eigenVector'

x1 = [ 2
       2 ]
e1=map*x1
len1=norm(e1)

x2 = [ -1
        1 ]

e2=map*x2
len2=norm(e2)