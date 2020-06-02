A_0=A3;
b_0=b3;


[m n]= size(A_0);
m
n

 
%l1 norm minimisation

f=[zeros(n,1);ones(m,1)];

A = [A_0,-eye(m);-A_0,-eye(m)];
b = [b_0; -b_0];

options = optimoptions('linprog','Algorithm','dual-simplex','Display','final');

tic
[solution,fval] = linprog(f,A,b);
toc
solution;
fval

x = solution(1:n);
r = A_0*x-b_0;
norm(r,1)
histogram(r,40)
xlabel('r')
ylabel('Number of residuals')