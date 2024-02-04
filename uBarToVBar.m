function Y = uBarToVBar(a,mixStartPt,Axintersec,VvalBar)

    targetFun = @(uS) Rho1(mixStartPt,a) - Rho1([uS;VvalBar],a);

    tolerance = 1e-12;

    if targetFun(mixStartPt(1))*targetFun(Axintersec) < 0

       Y = fzero(targetFun,[mixStartPt(1),Axintersec]);

    elseif abs(targetFun(mixStartPt(1))) < tolerance
       Y = mixStartPt(1);
    elseif abs(targetFun(Axintersec)) < tolerance
       Y = Axintersec;
    end




end