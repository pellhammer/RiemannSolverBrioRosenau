function Y = muF(uv,a)

uL = uv(1,:);
vL = uv(2,:);


rts = sort(roots([-2*(a-2)*vL*a, 0, (4-8*a)*vL,-8*(a-1)*uL , 2*vL ]));
Y = rts(2);


end