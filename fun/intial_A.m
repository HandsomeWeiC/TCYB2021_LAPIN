% initialization for dictionary A (two ways)!
function [A,m,order]=intial_A(X,gamma,b)

if b==2   %gamma 为整数
    [d,n]=size(X);
    %select the gamma of original points as the dictionary
    col=randperm(n);
    m=gamma;
    order=sort(col(1:m));
    A=X(:,order);   
elseif b==1
    [d,n]=size(X);
    %select the gamma of original points as the dictionary
    col=randperm(n);
    m=floor(gamma*n);%fix(0.8*n);  % 80% 是变量！！！！！！！

    order=sort(col(1:m));
%     order=cell2mat(struct2cell(load('order.mat')));
%     order=1:100;
    A=X(:,order);
end


end

%     % othogonalize X to get A , refer to the code of LRR
%     Q=orth(X');
%     A=X*Q;
%     m=d;
