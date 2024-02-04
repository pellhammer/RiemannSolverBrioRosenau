function Y = muE(uv,a)

uL = uv(1,:);
vL = uv(2,:);

format long



lambda1 = @(u,v) ((a+1)*u-sqrt((a-1).^2.*u.^2+4.*v.^2));
lambda = lambda1(uL,vL);

rts = muTildeT([uL,vL],a,lambda);
Y = rts(2);


end