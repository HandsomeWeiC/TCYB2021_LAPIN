%求解系数矩阵 Z
function Z = solveZ(X,A,F,lambda,c)
% X:m*n A:m*d Z:d*n  F:（d+n)*k  c：Z的稀疏约束
[m,n] = size(X);
[m,d] = size(A);
Da = A'*A;
U = F(1:d,:)';

V = F(d+1:d+n,:)';
H = L2_distance_1(U,V);
H = real(H);
B = lambda*H-2*A'*X;
Z = zeros(d,n);
obj=zeros(1,n);
for i=1:n
    b = B(:,i);
%     Z(:,i) = solvesub(Da,b,d,c,1);%矩阵会变奇异？
    Z(:,i) = NEWsolvesub(Da,b,d,c,1);
    obj(i)=norm(X(:,i)-A*Z(:,i),2)^2;
end
obj;

end

