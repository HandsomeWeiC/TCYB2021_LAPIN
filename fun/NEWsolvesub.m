%求解子问题 min_x x'Dx+b'x  s.t. x'1=1,x>=0, ||x||_0<=k
% min_x x'(Av+x)+mu/2*||x-v+y/mu||^2  s.t. x'1=1,x>=0, ||x||_0<=k
function x = NEWsolvesub(D,b,d,c,lambda)

% d：x的维数
x = zeros(d,1);
v = x; y = v;%拉格朗日乘子
mu = 1e-6;
rho = 1.5;
iternum = 100;
%目标值
objv = v'*D*v+b'*v+mu/2*norm(x-v+y/mu,2);
%更新
for iter=1:iternum
    %更新x，固定v
    w = v-y/mu;

    [x,S] = GSSP1(w,lambda,c);%ICML 的方法1和 老师的方法2
    % 不加稀疏约束
%     x = EProjSimplex_new(w);
    
    %更新v，固定x
    temp = 2*D+mu*eye(d);
    
    v = temp\(mu*x+y-b); 
    %更新乘子mu y
     y = y+mu*(x-v);
    mu = rho*mu;    
    %目标值
    objv(iter+1) = v'*D*v+b'*v+mu/2*norm(x-v+y/mu,2);
    leq(iter)=norm(x-v,2)/norm(v,2);
    %判决条件
    if norm(x-v,2)/norm(v,2) < 1e-7
%          disp(norm(x-v,2)/norm(v,2));
        break;
    end
end

end