1;

function s = matric(x)
  s0 = 0.1;
  s1 = 1.0;
  s = s0 * (1-x) + s1 * x; 
endfunction

s0_0 = matric(0)
s1_0 = matric(1)
s0_5 = matric(0.5)

x = 0;
ratio = 0;
step = 0;
w = 0.5
while (x < 1)
  step+=1
  x0 = x
  s0 = matric(x0)
  x1 = min(1.0,x0+w*s0)
  s1 = matric(x1)
  dx = (x1-x0)
  x += dx
  ratio += 0.5*(s0+s1)*dx
endwhile
