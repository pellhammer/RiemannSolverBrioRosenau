function Y = implicitHugoniotFun(Ul,U,a)

ul =Ul(1);
vl =Ul(2);
u =U(1);
v =U(2);




% hugoniotfunction with s eliminated 
% it is positive below the whole hug curve, positive inside the loop

Y = 2*u^2*v + 2*ul^2*vl + v*vl^2 + v^2*vl - v^3 - vl^3 - a*u^2*v + a*u^2*vl + a*ul^2*v - a*ul^2*vl - 2*u*ul*v - 2*u*ul*vl;



end