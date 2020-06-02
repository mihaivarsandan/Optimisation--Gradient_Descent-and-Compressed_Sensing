A_0=A3;
b_0=b3;


[m n]= size(A_0);
m;
n;

 
%l1 norm minimisation

f=[zeros(n,1);ones(m,1)];

A = [A_0,-eye(m);-A_0,-eye(m)];
b = [b_0; -b_0];
t=1;
f=[zeros(n,1);ones(m,1)];

A_opt =[A,-eye(m+m) ; -A,-eye(m+m)];
b_opt =[b;-b];
f_opt =[zeros(size(A,2),1);ones(size(b))];
[solution,fval] = linprog(f_opt,A_opt,b_opt);
x =solution(1:m+n);
objective_function(x,f,t,b,A)


function func = objective_function(x,c,t,b,A)
    func = t*c.'*x - sum(log(b-A*x)); 
end