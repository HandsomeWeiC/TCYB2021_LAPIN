%��������� min_x x'Dx+b'x  s.t. x'1=1,x>=0, ||x||_0<=k
% min_x x'(Av+x)+mu/2*||x-v+y/mu||^2  s.t. x'1=1,x>=0, ||x||_0<=k
function x = NEWsolvesub(D,b,d,c,lambda)

% d��x��ά��
x = zeros(d,1);
v = x; y = v;%�������ճ���
mu = 1e-6;
rho = 1.5;
iternum = 100;
%Ŀ��ֵ
objv = v'*D*v+b'*v+mu/2*norm(x-v+y/mu,2);
%����
for iter=1:iternum
    %����x���̶�v
    w = v-y/mu;

    [x,S] = GSSP1(w,lambda,c);%ICML �ķ���1�� ��ʦ�ķ���2
    % ����ϡ��Լ��
%     x = EProjSimplex_new(w);
    
    %����v���̶�x
    temp = 2*D+mu*eye(d);
    
    v = temp\(mu*x+y-b); 
    %���³���mu y
     y = y+mu*(x-v);
    mu = rho*mu;    
    %Ŀ��ֵ
    objv(iter+1) = v'*D*v+b'*v+mu/2*norm(x-v+y/mu,2);
    leq(iter)=norm(x-v,2)/norm(v,2);
    %�о�����
    if norm(x-v,2)/norm(v,2) < 1e-7
%          disp(norm(x-v,2)/norm(v,2));
        break;
    end
end

end