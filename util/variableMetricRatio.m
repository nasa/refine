1;

function m = metric(x)
  m0 = (0.1)^-2;
  m1 = (1.0)^-2;
  m = m0 * (1-x) + m1 * x;
endfunction

s0_0 = metric(0)
s1_0 = metric(1)
s0_5 = metric(0.5)

x = 0;
ratio = 0;
step = 0;
w = 1.0
while (x < 1)
  step+=1
  x0 = x
  m0 = metric(x0)
  s0 = 1/sqrt(m0)
  x1 = min(1.0,x0+w*s0)
  m1 = metric(x1)
  s1 = 1/sqrt(m1)
  dx = (x1-x0)
  x += dx
  ratio += 0.5*(s0+s1)*dx
endwhile
