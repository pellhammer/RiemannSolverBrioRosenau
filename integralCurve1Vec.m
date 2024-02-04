function Y = integralCurve1Vec(Ul,Ur,a,shift,tau)


ul = Ul(1);
vl = Ul(2);
tau = tau-shift;

Y=NaN(2,length(tau));

targetFun = @(muS) Rho1([ul;vl],a) - Rho1([muS;0],a);
uAxintersec = fzero(targetFun,1);

for kk =1:length(tau)

    Y(:,kk) = interalCurve1(Ul,Ur,a,uAxintersec,tau(kk));

end

end