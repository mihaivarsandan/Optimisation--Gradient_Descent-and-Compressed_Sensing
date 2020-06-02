A=A5;
b=b5;


[m n]= size(A);
m
n

 
%l2 norm minimisation


tic
x_opt = linsolve(A,b);
r= A*x_opt - b
fval=sqrt(sum(r.^2));
toc
r;
fval


histogram(r,50)
xlabel('r')
ylabel('Number of residuals')