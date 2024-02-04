function [Y1,Y2] = riemannSolverGodunov2x2(Ul,Ur,a)
    shift=0;
    


Yh = NaN(size(Ul));

for kk = 1:length(Ul(1,:))


    Ulrun = Ul(:,kk);
    Urrun = Ur(:,kk);

    [~,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersection(Ulrun,Urrun,a);



    if 0<=speed1(1)
        Yh(:,kk) = Ulrun;

    elseif length(speed1) == 1 && (speed1<0  &&  0 <=speed2(1)) % first wave is shock wave
        Yh(:,kk) = statesWave1(:,2);

    elseif length(speed1) == 2 && (speed1(1)<0  &&  0 <=speed1(2)) % first wave is rarefaction

        if strcmp(type1,'R1')
            Yh(:,kk) = integralCurve1Vec(statesWave1(:,1),statesWave1(:,2),a,shift,0);
        elseif strcmp(type1,'RS1')
            Yh(:,kk) = integralCurve1Vec(statesWave1(:,1),statesWave1(:,2),a,shift,0);
        end


    elseif speed1(end)<0  &&  0 <=speed2(1)
    
        Yh(:,kk) = statesWave2(:,1);

    elseif length(speed2) == 1 && speed2<0 % second wave is shock wave

        Yh(:,kk) = Urrun;

    elseif length(speed2) == 2 && speed2(1)<0  &&  0 <=speed2(2) % second wave is rarefaction
        if strcmp(type2,'R2')

            Yh(:,kk) = integralCurve2Vec(statesWave2(:,1),statesWave2(:,2),a,shift,0);
        elseif strcmp(type2,'SR2')
            Yh(:,kk) = integralCurve2Vec(statesWave2(:,2),statesWave2(:,3),a,shift,0);
        end

    elseif speed2(end)<=0
        
        Yh(:,kk) = Urrun;

    end

    Y1 = Yh(1,:);
    Y2 = Yh(2,:);
    
end




















end