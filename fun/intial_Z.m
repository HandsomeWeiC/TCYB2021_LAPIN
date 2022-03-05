% initialization for dictionary  (two ways)!
function Z=intial_Z(m,n,b)
if b==1
    % to let the element of matrix all be 0;
    Z=zeros(m,n);   
else
    %to let the sum of each column in Z be 1, satisfy Z^T*1=1  
    Z=rand(m,n);
    s=sum(Z);
    for i=1:n
        Z(:,i)=Z(:,i)/s(i);
    end
end
end