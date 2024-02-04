function Y = integralCurve2Vec(Ul,Ur,a,shift,tau)

ul = Ul(1);
vl = Ul(2);

tau=tau-shift;


Y=NaN(2,length(tau));
    
targetFun = @(muS) Rho2([ul;vl],a) - Rho2([muS;0],a);
uAxintersec = fzero(targetFun,0);


for kk =1:length(tau)
   
    Y(:,kk) = integralCurve2(Ul,Ur,a,uAxintersec,tau(kk));

end




end