function Y = muA(uv,a)

uL = uv(1,:);
vL = uv(2,:);

Y = -(1/(4*vL))*(2*(a-1)*uL-sqrt(4*(a-1)^2*uL^2+16*vL^2));

end