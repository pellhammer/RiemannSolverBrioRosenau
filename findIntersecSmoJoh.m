function [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2]  = findIntersecSmoJoh(Ul,Ur,a)

tolerance = 1e-14;

ul = Ul(1);
vl = Ul(2);
ur = Ur(1);
vr = Ur(2);




Rho1Val_l = Rho1([ul;vl],a);
Rho1Val_r = Rho1([ur;vr],a);

Rho2Val_l = Rho2([ul;vl],a);
Rho2Val_r = Rho2([ur;vr],a);



%% type 1 smoller johnson

if vl<0 && vr<=0 % Smoller Johnson Data
        
    if  (Rho1Val_r <= Rho1Val_l +tolerance)  &&  (Rho2Val_r <= Rho2Val_l + tolerance) % true if ur is to the right of the 1 rarefaction and 2 rarefaction
        
        type1 = 'R1';
        type2 = 'R2';
        
             % calculate the intersection of the rarefaction curve with the u axis   
             targetFun = @(muS) Rho1Val_l - Rho1([muS;0],a);
             uAxintersec = fzero(targetFun,1);

             targetFun = @(v) Rho2([uBarToVBar(a,[ul;vl],uAxintersec,v);v],a) - Rho2Val_r;

            if targetFun(vl)*targetFun(0)<0
                vintersec = fzero(targetFun,[vl,0]);
            elseif abs(targetFun(vl))<10e-10
                vintersec = vl;
            elseif abs(targetFun(0))<10e-10
                vintersec = 0;
            else
                Ul
                Ur
                targetFun(vl)
                targetFun(0)
               error( 'Weird')
            end
            intersecPoint = [uBarToVBar(a,[ul;vl],uAxintersec,vintersec);vintersec];
            speed1 = [lambda1(ul,vl,a),lambda1(intersecPoint(1),intersecPoint(2),a)];
            speed2 = [lambda2(intersecPoint(1),intersecPoint(2),a),lambda2(ur,vr,a)];
            statesWave1 = [Ul,intersecPoint];
            statesWave2 = [intersecPoint,Ur];

    elseif ur<=ul  && (implicitHugoniotFun([ul;vl],[ur;vr],a) <=0) % ul left of ur and HugoniotFun negative, ur is in region left 
        type1 = 'S1';
        type2 = 'S2';

        %parametrisiere kurve durch (ul,vl) durch mu.
        %find first a initial gues for the left intervall boundary by
        %finding the intersection with the rarefaction curve of (ur,vr)
        
        muAA = muA([ul;vl],a); 
        targetFun = @(muS) Rho2([U([ul;vl],a,muS);V([ul;vl],a,muS)],a) - Rho2Val_r;
        
        muintersec = fzero(targetFun,muAA-1);

        intersecParPoint = [U([ul;vl],a,muintersec);V([ul;vl],a,muintersec)];
        if intersecParPoint(1)>ul+10e-8
            ul
            vl
            ur
            vr
            error('Weird,not correctly computed')
        end
        
        %find actual intersection point
        targetFun = @(muS) implicitHugoniotFun([ur;vr],[U([ul;vl],a,muS);V([ul;vl],a,muS)],a);

        if targetFun(muintersec)*targetFun(muAA)<0
                muIntActual = fzero(targetFun,[muintersec,muAA]);
        elseif abs(targetFun(muintersec))<10e-10
                muIntActual  = muintersec;
        elseif abs(targetFun(muAA))<10e-10
               muIntActual  = muAA;
        else
                targetFun(muintersec)
                targetFun(muAA)
                Ul
                Ur
                error('Weird');
        end

        intersecPoint = [U([ul;vl],a,muIntActual);V([ul;vl],a,muIntActual)];
        speed1 = lambda1((1/2)*(ul+intersecPoint(1)),(1/2)*(vl+intersecPoint(2)),a);
        speed2 = lambda2((1/2)*(ur+intersecPoint(1)),(1/2)*(vr+intersecPoint(2)),a);
        statesWave1 = [Ul,intersecPoint];
        statesWave2 = [intersecPoint,Ur];
    elseif vr<=vl &&  (Rho2Val_r >= Rho2Val_l-tolerance) && (implicitHugoniotFun([ul;vl],[ur;vr],a) >=-tolerance) % vr<vl and pos hugFun, ur in the lower area

       
        type1 = 'S1';
        type2 = 'R2';
        
        %parametrisiere kurve durch (ul,vl) durch mu.
        %finding the intersection with the rarefaction curve of (ur,vr)
        
        muAA = muA([ul;vl],a); 
        targetFun = @(muS) Rho2([U([ul;vl],a,muS);V([ul;vl],a,muS)],a) - Rho2Val_r;
        
        muintersec = fzero(targetFun,muAA-0.1);

%         intersecPoint = [U([ul;vl],a,muintersec);V([ul;vl],a,muintersec)];
        kk = 0.1;
        while targetFun(muAA-kk)*targetFun(muAA)>=0 %loop only to find a starting intervall for fzero

            kk = kk+1;

        end

        muintersec = fzero(targetFun,muAA-kk);
        intersecPoint = [U([ul;vl],a,muintersec);V([ul;vl],a,muintersec)];

        speed1 = lambda1((1/2)*(ul+intersecPoint(1)),(1/2)*(vl+intersecPoint(2)),a);
        speed2 = [lambda2(intersecPoint(1),intersecPoint(2),a),lambda2(ur,vr,a)];
        statesWave1 = [Ul,intersecPoint];
        statesWave2 = [intersecPoint,Ur];



    elseif vr>=vl && (Rho1Val_r >= Rho1Val_l-tolerance) && (implicitHugoniotFun([ul;vl],[ur;vr],a) >=-tolerance) % inside the loop

        type1 = 'R1';
        type2 = 'S2';
        
        % calculate the intersection of the rarefaction curve with the u axis   
        targetFun = @(muS) Rho1Val_l - Rho1([muS;0],a);
        uAxintersec = fzero(targetFun,1);
        targetFun = @(v) implicitHugoniotFun([ur;vr],[uBarToVBar(a,[ul;vl],uAxintersec,v);v],a);


        if targetFun(vl)*targetFun(vr)<0
                vintersec = fzero(targetFun,[vl,vr]);
        elseif abs(targetFun(vl))<10e-10
                vintersec = vl;
        elseif abs(targetFun(vr))<10e-10
                vintersec = vr;
        else
                'Weird'
        end
        intersecPoint = [uBarToVBar(a,[ul;vl],uAxintersec,vintersec);vintersec];
        speed1 = [lambda1(ul,vl,a),lambda1(intersecPoint(1),intersecPoint(2),a)];
        speed2 = lambda2((1/2)*(ur+intersecPoint(1)),(1/2)*(vr+intersecPoint(2)),a);
        statesWave1 = [Ul,intersecPoint];
        statesWave2 = [intersecPoint,Ur];
    end

end































end