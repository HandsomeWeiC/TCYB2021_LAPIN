function [beta_star, S_star] = GSSP1(w, lambda, k)
%min||x-w||^2 s.t.x>=0,x'1=lambda,||x||_0<=k
[~, ind] = sort(w, 'descend');
S_star = ind(1:k);

beta_star = zeros(size(w));
w1 = w(S_star);
rho=1;
 for j=1:k
     if w1(j)-(sum(w1(1:j))-lambda)/j > 0
         rho = j;
     end
 end

tau = (1/rho)*(lambda- sum(w1(1:rho)));
beta_star(S_star) = max(w(S_star) + tau.*ones(size(w(S_star))), 0);
end