1;

function m = metric(x)
  m0 = (0.001)^-2;
  m1 = (1.0)^-2;
  m = m0 * (1-x) + m1 * x;
endfunction

dx = 1;
s0_0 = sqrt(dx*metric(0)*dx)
s1_0 = sqrt(dx*metric(1)*dx)

function integrated_length(w)
  w
  x = 0;
  ratio = 0;
  step = 0;
  while (x < 1)
    step+=1;
    x0 = x;
    m0 = metric(x0);
    s0 = 1/sqrt(m0);
    x1 = min(1.0,x0+w*s0);
    m1 = metric(x1);
    s1 = 1/sqrt(m1);
    dx = (x1-x0);
    x += dx;
    m = 0.5*(m0+m1);
    ratio += sqrt(dx*m*dx);
  endwhile
  steps=step
  ratio
endfunction

integrated_length(1000)
integrated_length(1)
integrated_length(.1)
integrated_length(.05)
