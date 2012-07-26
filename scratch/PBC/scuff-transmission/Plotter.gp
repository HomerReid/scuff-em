Eps=10.0;
Mu=1.0;
n = sqrt(Eps*Mu);
g =sqrt(Eps/Mu); # gamma
r(x) = (g*g-1)*sin(n*x) / ( (1+g*g)*sin(n*x) + 2*i*g*cos(n*x))
t(x) = 2*i*g*exp(-i*x) / ( (1+g*g)*sin(n*x) + 2*i*g*cos(n*x))
