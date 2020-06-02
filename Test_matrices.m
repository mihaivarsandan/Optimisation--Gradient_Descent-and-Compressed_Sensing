randn('state',0);
rand('state',0);
m = 100;
n = 50;
A = randn(m,n);
b = rand(m,1);
c = A'*rand(m,1); 

% Solve using linprog.
x = linprog(c, A, b);
optval = c'*x;

% Make a change of variables  x := x + s*c, so that optimal value is 1.
s = (1-optval)/(c'*c);
b = b + s*A*c;
x0 = s*c; 


% xc is the point on the central path with barrier parameter t=1.
t = 1;  
x = x0;
d = b-A*x;
f = t*c'*x-sum(log(d))


