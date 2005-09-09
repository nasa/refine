1;

function row = form_row(edge)

  % P' M P

  % m11 m12 m13
  % m12 m22 m23
  % m13 m23 m33

  % x*m11+y*m12+z*m13 x*m12+y*m22+z*m23 x*m13+y*m23+z*m33

  % x(x*m11+y*m12+z*m13) + y(x*m12+y*m22+z*m23) + z(x*m13+y*m23+z*m33)

  % xxm11 + xym12 + xzm13 + yxm12 + yym22 + yzm23 + zxm13 + zym23 + zzm33

  % xx 2xy 2xz yy 2yz zz

  row = [  edge(1)*edge(1)
	 2*edge(1)*edge(2)
	 2*edge(1)*edge(3)
	   edge(2)*edge(2)
	 2*edge(2)*edge(3)
	   edge(3)*edge(3)
	 ]';

endfunction

function m = implied_metric(n1,n2,n3,n4)

  edge1 = n2-n1;
  edge2 = n3-n1;
  edge3 = n4-n1;
  edge4 = n3-n2;
  edge5 = n4-n2;
  edge6 = n4-n3;

  equ = [ form_row(edge1)
        form_row(edge2)
        form_row(edge3)
        form_row(edge4)
        form_row(edge5)
        form_row(edge6)
  ];

  m = equ \ [1;1;1;1;1;1];
endfunction

n1 = [ 0 0 0 ];
n2 = [ 1 0 0 ];
n3 = [ 0 1 0 ];
n4 = [ 0 0 1 ];

implied_metric(n1,n2,n3,n4)

n1 = [ 0.000, 0.000, 0.000 ];
n2 = [ 1.000, 0.000, 0.000 ];
n3 = [ 0.500, 0.866, 0.000 ];
n4 = [ 0.500, 0.289, 0.823 ];

implied_metric(n1,n2,n3,n4)

n1 = 0.5*[ 0.000, 0.000, 0.000 ];
n2 = 0.5*[ 1.000, 0.000, 0.000 ];
n3 = 0.5*[ 0.500, 0.866, 0.000 ];
n4 = 0.5*[ 0.500, 0.289, 0.823 ];

implied_metric(n1,n2,n3,n4)

delta = 1.0e-6;
n1 = [-1.000, 0.000, 0.000 ];
n2 = [ 1.000, 0.000, 0.000 ];
n3 = [ 0.000, delta, 0.000 ];
n4 = [ 0.000, 0.000, delta ];

implied_metric(n1,n2,n3,n4)

delta = 1.0e-6;
n1 = [ 0.000, 0.000, 0.000 ];
n2 = [ 1.000, 0.000, 0.000 ];
n3 = [ 0.000, 1.000, 0.000 ];
n4 = [ 0.000, 0.000, delta ];

implied_metric(n1,n2,n3,n4)
