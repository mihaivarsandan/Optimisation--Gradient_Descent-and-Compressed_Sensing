A_0=A3;
b_0=b3;


[m n]= size(A_0);
m
n

 
%l1 norm minimisation

f=[zeros(n,1);ones(m,1)];

A = [A_0,-eye(m);-A_0,-eye(m)];
b = [b_0; -b_0];



beta=0.5;
alpha=0.01;
t=1;
nr=0;
x=x0;
smin=0.5;
mu=2;
function_values=[]
 while norm(gradient(x,f,t,b,A))>0.0013
  s=1;
  nr = nr + 1;
  while min(b-A*(x-s*gradient(x,f,t,b,A))) < 0  
      s = beta*s; 
  end
  while objective_function(x-s*gradient(x,f,t,b,A),f,t,b,A)>objective_function(x,f,t,b,A)-alpha*s*norm(gradient(x,f,t,b,A))
    s = beta * s;
  end
   
   x = x - s*gradient(x,f,t,b,A);
   function_values=[function_values;objective_function(x,f,t,b,A)];
   if s >= smin
         t=mu*t;
       
   
 end
 
 
 nr
 objective_function(x,f,t,b,A)
 f'*x
 norm(gradient(x,f,t,b,A))
 
 
 function_values=function_values(1:(size(function_values)-1));
 f_prime_values =ones(size(function_values)) * objective_function(x,f,t,b,A);
 
 c= function_values -f_prime_values;
 k= 1:(nr-1);
 semilogy(k,c)
 xlabel('Iteration $k$','interpreter','latex')
 ylabel('$f(\mathbf{\tilde{x}}_k)-f(\mathbf{\tilde{x}}^*)$','interpreter','latex')
 
 
 
function func = objective_function(x,c,t,b,A)
    func = t*c.'*x - sum(log(b-A*x)); 
end

function delta= gradient(x,c,t,b,A)
 v = b-A*x;
 delta = t*c + A'*(1./v);
 end

