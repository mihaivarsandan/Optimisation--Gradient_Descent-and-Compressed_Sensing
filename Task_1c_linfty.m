A_0=A5;
b_0=b5;


[m n]= size(A_0);
m
n

 
%l_infty norm minimisation

f=[zeros(n,1);1];

A = [A_0,-ones(m,1);-A_0,-ones(m,1)];
b = [b_0; -b_0];

options = optimoptions('linprog','Algorithm','dual-simplex','Display','final');

tic
[solution,fval] = linprog(f,A,b);
toc
solution;
fval
x = solution(1:n);

r = A_0*x - b_0;
histogram(r,40)
xlabel('r')
ylabel('Number of residuals')

