function [y2,idx,Z,A] = LAPIN(X,k,c,lambda,gamma)

%% input
itermax = 500; %最大迭代步数
tol = 10e-4; %收敛精度
% X = rand(5,50);
%X = (X-min(min(X)))/(max(max(X))-min(min(X)));
% k = 2;  %聚类个数
% c = 2;  %系数矩阵一列有且只有c个非零值  ？？？？
% lambda = 10;  %超参 
% gamma = 0.5;  % 字典选取参数

 %% 初始化
   [m,n] = size(X);
   % 字典 A
   [A,d,order] = intial_A(X,gamma,1); %  m*d
   % 系数矩阵 Z
    Z = intial_Z(d,n,2);  % d*n
%    Z = cell2mat(struct2cell(load('Z.mat')));
   
   % 矩阵 F
   G=zeros(d+n,d+n);
   G(1:d,d+1:d+n) = Z; G(d+1:d+n,1:d) = Z';
   d1 = sum(G,2); D=diag(d1); 
   Lg = D-G;  % (d+n)*(d+n)
   [V,ev0,ev]=eig1(Lg,k,0);  % 要最小的特征值和特征向量
   F0 = V;
   if sum(ev(1:k+1)) < tol
      error('The original graph has more than %d connected component', k);
   end;   %y有问题  随机初始，可能开始是不对的
   if sum(ev0) < tol  % 有k个连通分量
       Gs = sparse(G);
       [~,y] = graphconncomp(Gs);
       y1 = y(1:d)';
       y2 = y(d+1:end)';
       return;
   end
   
   % 目标函数值
   objv = norm(X-A*Z,'fro')+lambda*(trace(F0'*Lg*F0));
   
%% 主循环
   F = F0;
   for iter = 1:itermax
       %更新系数矩阵 Z
       
       Z = solveZ(X,A,F,lambda,c); %求解系数矩阵Z    ?2 c怎么确定

       % 更新字典矩阵 A
        A = X*Z'/(Z*Z'+eye(d)*0.0001);
%        A= cell2mat(struct2cell(load('toydictionary.mat')));
       
       % 更新矩阵 特征向量矩阵F
       G=zeros(d+n,d+n);
       G(1:d,d+1:d+n) = Z; G(d+1:d+n,1:d) = Z';
       d1 = sum(G,2); D = diag(d1); 
       Lg = D-G;  % (d+n)*(d+n)
      [V,ev0,ev]=eig1(Lg,k,0);  % 要最小的特征值和特征向量
       F = V;
       
      % 目标函数值
      objv(iter+1) = norm(X-A*Z,'fro')+lambda*(trace(F'*Lg*F));
       
      % 判决条件
      fn1 = sum(ev0);  % ？验证fn1 == sum(ev(1:k))
      fn2 = sum(ev(1:k+1));
      
      if fn1 < tol && fn2 > tol
          break;
      elseif abs(objv(iter+1)-objv(iter))/objv(iter) < tol
          break;
      end      
%       if abs(objv(iter+1)-objv(iter)) < tol
%           break;
%       end  
   end
    
 % 得到分类结果

Gs = sparse(G);
[clusternum,y] = graphconncomp(Gs);
% imshow(Z);
y1 = y(1:d)';
y2 = y(d+1:end)';
idx = kmeans(Z',k,'emptyaction','singleton','replicates',20,'display','off');
if clusternum ~= k
    sprintf('Can not find the correct cluster number: %d', k)
end
end

  
  
  
  
  
  
  
  
   
   
   