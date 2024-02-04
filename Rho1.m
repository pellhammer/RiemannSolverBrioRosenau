function Y = Rho1(U,a)

u = U(1,:);
v = U(2,:);

%first Riemann invariant Rho1
Y = (lambda1(u,v,a) - 4*u).^((a-3)/(2*(a-2))) .* (lambda1(u,v,a) - 2*a*u).^((a-1)/(2*(a-2)));

%second Riemann invariant
% rho2 = (lambda2(u,v,a) - 4*u).^((a-3)/(2*(a-2))) * (lambda2(u,v,a) - 2*a*u).^((a-1)/(2*(a-2)));




end