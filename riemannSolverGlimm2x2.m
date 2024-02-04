function [Y1,Y2] = riemannSolverGlimm2x2(Ul,Ur,a,evalPT)

Yh = NaN(size(Ul));

    for kk = 1:length(Ul(1,:))
        Ul(:,kk);
        Ur(:,kk);
        Yh(:,kk) = evaluateRiemannSol2x2(Ul(:,kk),Ur(:,kk),a,evalPT(kk));

    end
    

    Y1 = Yh(1,:);
    Y2 = Yh(2,:);
end