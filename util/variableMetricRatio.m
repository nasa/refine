1;

function m = metric(x,s0,s1)
  m0 = (s0)^-2;
  m1 = (s1)^-2;
  m = m0 * (1-x) + m1 * x;
endfunction

function ratio = integrated_length(w,s0,s1)
  w
  x = 0;
  ratio = 0;
  step = 0;
  while (x < 1)
    step+=1;
    x0 = x;
    m0 = metric(x0,s0,s1);
    s0 = 1/sqrt(m0);
    x1 = min(1.0,x0+w*s0);
    m1 = metric(x1,s0,s1);
    s1 = 1/sqrt(m1);
    dx = (x1-x0);
    x += dx;
    m = 0.5*(m0+m1);
    ratio += sqrt(dx*m*dx);
  endwhile
  steps=step
  ratio;
endfunction

s0 = 0.001
s1 = 1.0

dx = 1;
s0_0 = sqrt(dx*metric(0,s0,s1)*dx)
s1_0 = sqrt(dx*metric(1,s0,s1)*dx)

ratio  = integrated_length(1000,s0,s1)
ratio100 = integrated_length(1,s0,s1)
ratio010 = integrated_length(.1,s0,s1)
ratio005 = integrated_length(.05,s0,s1)

fudge = ratio005 / ratio

s0 = 0.1
s1 = 1.0

dx = 1;
s0_0 = sqrt(dx*metric(0,s0,s1)*dx)
s1_0 = sqrt(dx*metric(1,s0,s1)*dx)

ratio  = integrated_length(1000,s0,s1)
ratio100 = integrated_length(1,s0,s1)
ratio010 = integrated_length(.1,s0,s1)
ratio005 = integrated_length(.05,s0,s1)

fudge = ratio005 / ratio

s0 = 1.0
s1 = 1.0

dx = 1;
s0_0 = sqrt(dx*metric(0,s0,s1)*dx)
s1_0 = sqrt(dx*metric(1,s0,s1)*dx)

ratio  = integrated_length(1000,s0,s1)
ratio100 = integrated_length(1,s0,s1)
ratio010 = integrated_length(.1,s0,s1)
ratio005 = integrated_length(.05,s0,s1)

fudge = ratio005 / ratio
