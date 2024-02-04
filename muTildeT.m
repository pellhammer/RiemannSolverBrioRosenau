function Y = muTildeT(uv,a,c)
% find x s.d. c=sigma(x)


uL = uv(1);
vL = uv(2);


% hugoniot
nu = @(uv,mu) 2*((-uv(2)*mu.^2-uv(1)*mu*(a-1)+uv(2))./((a-2)*mu.^2+1));
u =  @(uv,mu) uv(1) + mu.*nu(uv,mu);
v =  @(uv,mu) uv(2) - nu(uv,mu);
s =  @(uv,mu) -2*uv(2)*mu + 2*u(uv,mu);
%coeffs = [a*vL, (-2*c + a*c + a*uL) , -vL,-uL+c];
coeffs = [2*a*vL, - 2*c + a*c + 2*a*uL, -2*vL, c - 2*uL];


rts = sort(roots(coeffs));
%rts(check)=[];
rts(imag(rts)~=0)=[];

Y = rts;

end