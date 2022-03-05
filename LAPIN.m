function [y2,idx,Z,A] = LAPIN(X,k,c,lambda,gamma)

%% input
itermax = 500; %����������
tol = 10e-4; %��������
% X = rand(5,50);
%X = (X-min(min(X)))/(max(max(X))-min(min(X)));
% k = 2;  %�������
% c = 2;  %ϵ������һ������ֻ��c������ֵ  ��������
% lambda = 10;  %���� 
% gamma = 0.5;  % �ֵ�ѡȡ����

 %% ��ʼ��
   [m,n] = size(X);
   % �ֵ� A
   [A,d,order] = intial_A(X,gamma,1); %  m*d
   % ϵ������ Z
    Z = intial_Z(d,n,2);  % d*n
%    Z = cell2mat(struct2cell(load('Z.mat')));
   
   % ���� F
   G=zeros(d+n,d+n);
   G(1:d,d+1:d+n) = Z; G(d+1:d+n,1:d) = Z';
   d1 = sum(G,2); D=diag(d1); 
   Lg = D-G;  % (d+n)*(d+n)
   [V,ev0,ev]=eig1(Lg,k,0);  % Ҫ��С������ֵ����������
   F0 = V;
   if sum(ev(1:k+1)) < tol
      error('The original graph has more than %d connected component', k);
   end;   %y������  �����ʼ�����ܿ�ʼ�ǲ��Ե�
   if sum(ev0) < tol  % ��k����ͨ����
       Gs = sparse(G);
       [~,y] = graphconncomp(Gs);
       y1 = y(1:d)';
       y2 = y(d+1:end)';
       return;
   end
   
   % Ŀ�꺯��ֵ
   objv = norm(X-A*Z,'fro')+lambda*(trace(F0'*Lg*F0));
   
%% ��ѭ��
   F = F0;
   for iter = 1:itermax
       %����ϵ������ Z
       
       Z = solveZ(X,A,F,lambda,c); %���ϵ������Z    ?2 c��ôȷ��

       % �����ֵ���� A
        A = X*Z'/(Z*Z'+eye(d)*0.0001);
%        A= cell2mat(struct2cell(load('toydictionary.mat')));
       
       % ���¾��� ������������F
       G=zeros(d+n,d+n);
       G(1:d,d+1:d+n) = Z; G(d+1:d+n,1:d) = Z';
       d1 = sum(G,2); D = diag(d1); 
       Lg = D-G;  % (d+n)*(d+n)
      [V,ev0,ev]=eig1(Lg,k,0);  % Ҫ��С������ֵ����������
       F = V;
       
      % Ŀ�꺯��ֵ
      objv(iter+1) = norm(X-A*Z,'fro')+lambda*(trace(F'*Lg*F));
       
      % �о�����
      fn1 = sum(ev0);  % ����֤fn1 == sum(ev(1:k))
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
    
 % �õ�������

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

  
  
  
  
  
  
  
  
   
   
   