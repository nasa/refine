1;

function form_test(n,t,y,rSpace,tSpace,ySpace)
r=[ n t y ]


h=[ 
   1/rSpace 0 0
   0 1/tSpace 0
   0 0 1/ySpace
]

d=h*h

d_times_r = d*r

j = h*r'

m = r*d*r'

[V,LAMDA]=eig(m)

n(1)*n(1)*(1/rSpace/rSpace) + t(1)*t(1)*(1/tSpace/tSpace) + y(1)*y(1)*(1/ySpace/ySpace)

n(1)*n(2)*(1/rSpace/rSpace) + t(1)*t(2)*(1/tSpace/tSpace) + y(1)*y(2)*(1/ySpace/ySpace)

n(1)*n(3)*(1/rSpace/rSpace) + t(1)*t(3)*(1/tSpace/tSpace) + y(1)*y(3)*(1/ySpace/ySpace)

n(2)*n(2)*(1/rSpace/rSpace) + t(2)*t(2)*(1/tSpace/tSpace) + y(2)*y(2)*(1/ySpace/ySpace)

n(2)*n(3)*(1/rSpace/rSpace) + t(2)*t(3)*(1/tSpace/tSpace) + y(2)*y(3)*(1/ySpace/ySpace)

n(3)*n(3)*(1/rSpace/rSpace) + t(3)*t(3)*(1/tSpace/tSpace) + y(3)*y(3)*(1/ySpace/ySpace)

endfunction

n=[
1
0
0
]

t=[
0
0
1
]

y=[
0
1
0
]


rSpace = 0.1
tSpace = 0.25
ySpace = 0.5

r=[ n t y ]
form_test(n,t,y,rSpace,tSpace,ySpace)

n=[
sqrt(0.5)
0
sqrt(0.5)
]

t=[
-sqrt(0.5)
0
sqrt(0.5)
]

y=[
0
1
0
]

rSpace = 0.1
tSpace = 0.25
ySpace = 0.5

form_test(n,t,y,rSpace,tSpace,ySpace)

n=[
1
0
0
]

t=[
0
1
0
]

y=[
0
0
1
]

rSpace = 1
tSpace = 2
ySpace = 4

form_test(n,t,y,rSpace,tSpace,ySpace)
