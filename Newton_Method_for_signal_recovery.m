%Parameters Initialisation
x_0 = x;


A=A;
b = A*x_0;
% i = 1:size(b);
% plot(i,b)
t=1;
lambda_max = norm(2*A'*b,Inf);
lambda = 0.01 * lambda_max;
u =0.01 *ones(size(x_0));

setGlobalA(A);
setGlobalb(b);
setGloballambda(lambda);
setGlobalt(t);

% END of Parameters Initialisation
F(x_0,u);
g(x_0,u);
inv(H(x_0,u));



%Hyperaramets Initialisation
alpha = 0.01;
beta = 0.5;
epsilon = 10^-4;
x_1 = zeros(size(x));
maxiter = 5;
smin=0.5;
mu = 2;

 while sqrt(g(x_1,u)' * (H(x_1,u)\g(x_1,u))) > epsilon
     delta   = -H(x_1,u)\g(x_1,u);
     delta_x = delta(1:size(x_1));
     delta_u = delta(size(x_1)+1:end);
     s=1;
     for iter=1:maxiter
        x_new = x_1+s*delta_x; 
        u_new = u+s*delta_u;
        f = [x_new-u_new;-x_new-u_new];
        if (max(f) < 0)
            if (F(x_new,u_new)-F(x_1,u) <= alpha*s*g(x_1,u)'*delta)
                break;
            end
        end
        s = beta*s;
     end    
     x_1 = x_new;
     u = u_new;
   
     if s >= smin
         t=mu*t;
         setGlobalt(t);
     end
 end
i = 1:size(x_1);
sqrt(g(x_1,u)' * (H(x_1,u)\g(x_1,u)))
sqrt(g(x_0,u)' * (H(x_0,u)\g(x_0,u)))
norm(A*x_0-b) + lambda*(norm(x_0,1))
norm(A*x_1-b) + lambda*(norm(x_1,1))


 x_LS = (A'*A)\(A'*b);
 x_Tik = (A'*A+lambda*eye(256))\(A'*b);
%tiledlayout(2,1)
%nexttile
stem(x_1,'Marker','.')
xlim([0,256])
%ylim([-0.5,max(b)+0.05])
xlabel('Sample')
ylabel('Amplitude')

% nexttile
% stem(x_LS,'Marker','.')
% xlim([0,256])
% %ylim([-0.5,max(b)+0.05])
% xlabel('Sample')
% ylabel('Amplitude')

% nexttile
% stem(b)
% title('Compressed signal')
% 
% nexttile
% stem(i,x_1)
% xlim([0,256])
% title('Uncompressed Signal')
% 

% 
% nexttile
% stem(A*x_1)
% title('Minimum Energy Signal')

 








function f = F(x,u)
global b_g 
global A_g 
global t_g 
global lambda_g
f = t_g * norm(A_g*x-b_g)^2 + t_g*lambda_g * sum(u) -  sum(log(u.^2-x.^2));
end

function grad = g(x,u)
global b_g 
global A_g 
global t_g 
global lambda_g
g_x = 2 * t_g * A_g' * A_g * x - 2 * t_g * A_g' * b_g  + 2 * x./(u.^2 - x.^2);
g_u = t_g * lambda_g * ones(size(u)) - 2 * u./(u.^2-x.^2);
grad = [g_x ; g_u];
end

function hess=H(x,u) 
global A_g 
global t_g
D_1 = diag(2*(x.^2 + u.^2)./(u.^2-x.^2).^2);
D_2 = diag(-4*(x.*u)./(u.^2-x.^2).^2);
H_xx = 2 * t_g * A_g' * A_g + D_1;
H_xu = D_2;
H_ux = H_xu;
H_uu = D_1;

hess = [H_xx, H_xu ; H_ux, H_uu];
end

function setGlobalA(val)
global A_g
A_g = val;
end

function setGlobalb(val)
global b_g
b_g = val;
end

function setGlobalt(val)
global t_g
t_g = val;
end

function setGloballambda(val)
global lambda_g
lambda_g = val;
end

