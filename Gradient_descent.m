A_0=A3;
b_0=b3;


[m n]= size(A_0);
m
n

 
%l1 norm minimisation

f=[zeros(n,1);ones(m,1)];
t=1;
A = [A_0,-eye(m);-A_0,-eye(m)];
b = [b_0; -b_0];

F(x,f,t,b,A)
norm(D(x,f,t,b,A))
x=x0;
k=0;
j=0;
step=3;

beta =0.8;
alpha = 0.8;
epsilon =0.1;

while norm(D(x,f,t,b,A))>epsilon
    s = 1;
    while step > 0
        
        x_new = x - s/step * D(x,f,t,b,A);
        
        if imag(F(x_new,f,t,b,A)) ~= 0
            break
        end
        if F(x_new,f,t,b,A) < F(x,f,t,b,A) - s/step * alpha * norm(D(x,f,t,b,A))
            j = j + 1;
            x = x_new;
            step = step -1;
            s=1;
            
        
        else
            s = beta*s;
        end
    end
    k = k+1;
    j = 0;
    step =3;
    norm(D(x,f,t,b,A));
    F(x,f,t,b,A)
    if imag(F(x_new,f,t,b,A)) ~= 0
            break
    end
    
end
norm(D(x,f,t,b,A))
F(x,f,t,b,A)


function func = F(x,c,t,b,A)
    func = t*c.'*x - sum(log(b-A*x)); 
end

function delta= D(x,c,t,b,A)
 v = b-A*x;
 delta = t*c + A'*(1./v);
 end