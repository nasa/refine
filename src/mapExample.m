format long

p7 = 1.0/sqrt(2);

eVect1 = [  p7
            p7 ]
eVect2 = [ -p7
            p7 ]

eigenVector = [ eVect1 eVect2 ]

eVal1 = 0.5
eVal2 = 1.0

eigenValue  = [ eVal1   0.0
                  0.0 eVal2 ]

map = eigenValue * eigenVector

v1 = [ 2
       2 ]
map*v1

v2 = [ -1
        1 ]

map*v2
