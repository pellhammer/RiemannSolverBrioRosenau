function [intersecPoint,intersecInUpperPlane,speedConsec,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersectionInUpperHalfplane(Ul,Ur,a)
% find intersection in upper halfplane. Uses only values vr>0 vl<0
tolerance = 10e-12;
    
ul = Ul(1);
vl = Ul(2);
ur = Ur(1);
vr = Ur(2);


RhoVal2r = Rho2([ur;vr],a);

    

if vl <= -sqrt(a)*ul % left value inside of A_1 aor A_2
    

     % calculate the intersection of the rarefaction curve with the u axis   
    targetFun = @(muS) Rho1([ul;vl],a)-Rho1([muS;0],a);
    uAxintersec = fzero(targetFun,1);

    % calsulate alpha_f 
    intersecWithSqA = @(u) Rho1([u;sqrt(a)*u],a) - Rho1([ul;vl],a);
    alpha_fu = fzero(intersecWithSqA,[0,uAxintersec]);
    alpha_f = [alpha_fu;-sqrt(a)*alpha_fu];
    beta_f = -alpha_f;
    
    Rho2_betaf = Rho2(beta_f,a);

    if (RhoVal2r >= Rho2_betaf+tolerance) % (ur,vr) is to the left of the rarefaction curve trough alphaf
        
        intersecPoint = NaN;
        speed1 = NaN;
        speed2 = NaN;
        intersecInUpperPlane = false;
        speedConsec = false;
        type1 = 'noType';
        type2 = 'noType';
        statesWave1 = NaN;
        statesWave2 = NaN;
    elseif abs(RhoVal2r - Rho2_betaf)<tolerance % (ur,vr) is on rarefaction curve trough alphaf

        if vr>beta_f(2) % true if vr is over betaf
            type1 = 'RS1';
            type2 = 'R2';
            intersecPoint = beta_f; 
            intersecInUpperPlane = true; 
            speedConsec = true;
            ptOnRare = alpha_f;
            speed1 = [lambda1(ul,vl,a),lambda1(ptOnRare(1),ptOnRare(2),a)];
            speed2 = [lambda2(intersecPoint(1),intersecPoint(2),a),lambda2(ur,vr,a)];
    
            statesWave1 = [Ul,ptOnRare,intersecPoint];
            statesWave2 = [intersecPoint,Ur];  
        else
    
            intersecPoint = beta_f;
            intersecInUpperPlane = true;
            speedConsec = false;
            type1 = 'noType';
            type2 = 'noType';
            statesWave1 = NaN;
            statesWave2 = NaN; 
            speed1 = NaN;
            speed2 = NaN;
        end


    elseif (RhoVal2r <= Rho2_betaf-tolerance) % (ur,vr) is above to the right of the rarefaction through beta_f

        
        %first find intersec with rarefaction and mixedCurve

        targetFun = @(muS) Rho2(mixedCurveByBarV([ul,vl],a,alpha_f,uAxintersec,muS),a) - RhoVal2r;

        VvalBar=fzero(targetFun,[alpha_f(2),0]);%computes \bar{v} according to intersection
        intersecRarPoint = mixedCurveByBarV([ul,vl],a,alpha_f,uAxintersec,VvalBar);

        if intersecRarPoint(2)<=vr
            intersecPoint = intersecRarPoint;
            intersecInUpperPlane = true;
            speedConsec = true;
            type1 = 'RS1';
            type2 = 'R2';
            
            uValue = uBarToVBar(a,alpha_f,uAxintersec,VvalBar);

            ptOnRare = [uValue;VvalBar];
            speed1 = [lambda1(ul,vl,a),lambda1(ptOnRare(1),ptOnRare(2),a)];
            speed2 = [lambda2(intersecRarPoint(1),intersecRarPoint(2),a),lambda2(ur,vr,a)];
    
            statesWave1 = [Ul,ptOnRare,intersecPoint];
            statesWave2 = [intersecPoint,Ur]; 

        else
            % findIntersec of shock Curve and mixedCurve
            targetFun = @(muS) implicitHugoniotFun([ur,vr],mixedCurveByBarV([ul,vl],a,alpha_f,uAxintersec,muS),a);
            if targetFun(VvalBar)*targetFun(0)<0
                sol2 = fzero(targetFun,[VvalBar,0]);

            elseif abs(targetFun(VvalBar))<=10e-10 
                sol2 = VvalBar;
            else
               Ul
               Ur
                error('Weird1')
                targetFun(VvalBar)
                targetFun(0)
            end
            intersecPoint = mixedCurveByBarV([ul,vl],a,alpha_f,uAxintersec,sol2);
            intersecInUpperPlane = true;

            %find the value of \bar{u}
            uValue = uBarToVBar(a,alpha_f,uAxintersec,VvalBar);

            ptOnRare = [uValue;sol2];
            speed1 = lambda1(ptOnRare(1),ptOnRare(2),a);        
            speed2 = lambda2((1/2)*(intersecPoint(1) + ur), (1/2)*(intersecPoint(2)+vr),a);
            if speed1 > speed2 + tolerance
                speedConsec = false;
                type1 = 'noType';
                type2 = 'noType';
                statesWave1 = NaN;
                statesWave2 = NaN;     
            else
                speedConsec = true;
                type1 = 'RS1';
                type2 = 'S2';
                statesWave1 = [Ul,ptOnRare,intersecPoint];
                statesWave2 = [intersecPoint,Ur]; 
                speed1 = [lambda1(ul,vl,a),speed1];

            end
        end

    end


%% The case where (ul,vl) is in A_3

elseif vr>-sqrt(a)*ul % The case where (ul,vl) is in A_3

    % calculate the intersection of the rarefaction curve with the u axis   
    targetFun = @(muS) Rho1([ul;vl],a)-Rho1([muS;0],a);
    uAxintersec = fzero(targetFun,0);


    %find stat wher the forward wave curve is tangential to the 2
    %eigenvector
    muFF = muF([ul;vl],a);
    endOfWaveCurve = [U([ul;vl],a,muFF);V([ul;vl],a,muFF)];

    %find first stat wher the forward wave consists of 1shock
    muEE = muE([ul;vl],a);
    beginnOf1shock = [U([ul;vl],a,muEE);V([ul;vl],a,muEE)];

    Rho2_endofWave = Rho2(endOfWaveCurve,a);
    Rho2_beginnOf1shock = Rho2(beginnOf1shock,a);

    if (RhoVal2r >= Rho2_endofWave+tolerance) % (ur,vr) is to the left of the rarefaction curve trough endOfWaveCurve 

        intersecPoint = NaN;
        intersecInUpperPlane = false;
        speedConsec = false;
        type1 = 'noType';
        type2 = 'noType';
        statesWave1 = NaN;
        statesWave2 = NaN;
        speed1 = NaN;
        speed2 = NaN;

    elseif abs(RhoVal2r - Rho2_endofWave)<tolerance

        if vr>endOfWaveCurve(2) % true if vr is over betaf
            type1 = 'S1';
            type2 = 'R2';
            intersecPoint = endOfWaveCurve; 
            intersecInUpperPlane = true; 
            speedConsec = true;
            speed1 = S([ul;vl],a,muFF);
            speed2 = [lambda2(endOfWaveCurve(1),endOfWaveCurve(2),a),lambda2(ur,vr,a)];
        else
    
            intersecPoint = endOfWaveCurve;
            intersecInUpperPlane = true;
            speedConsec = false;
            type1 = 'noType';
            type2 = 'noType';
            statesWave1 = NaN;
            statesWave2 = NaN;
            speed1 = NaN;
            speed2 = NaN;
        end


    elseif (RhoVal2r < Rho2_endofWave)  && (RhoVal2r >= Rho2_beginnOf1shock) % (ur,vr) lies between rarefaction curve trough endOFWAveCurve and beginn1shockCurve
            %hence, the solution contains of a 1-shock in the upper
            %halfplane followed by a rarefaction or a shock 
            
            %first find intersec with rarefaction and the shockCurve
            targetFun = @(muS) Rho2([U([ul;vl],a,muS);V([ul;vl],a,muS)],a) - RhoVal2r;
            muOfInt = fzero(targetFun,[muEE,muFF]);%computes parameter mu of Intersection
            intersecRarPoint = [U([ul;vl],a,muOfInt);V([ul;vl],a,muOfInt)];
            

            if intersecRarPoint(2)<=(vr+tolerance)
                intersecPoint = intersecRarPoint;
                intersecInUpperPlane = true;
                speedConsec = true;
                type1 = 'S1';
                type2 = 'R2';
                speed1 = S([ul;vl],a,muOfInt);
                speed2 = [lambda2(intersecPoint(1),intersecPoint(2),a),lambda2(ur,vr,a)];
                statesWave1 = [Ul,intersecPoint];
                statesWave2 = [intersecPoint,Ur];               

            else
                
                % findIntersec of shock Curve and mixedCurve. Note that
                % this intersection point can also lie in the region to the
                % very right that is handled in the nex section. We have to
                % deal with this. We use another if clase
                targetFun = @(muS) implicitHugoniotFun([ur;vr],[U([ul;vl],a,muS);V([ul;vl],a,muS)],a);
%                 'bapedi'
%                 targetFun(muEE)
%                 targetFun(muOfInt)

                if abs(targetFun(muEE))<tolerance
                    sol2 = muEE;
                    intersecPoint = beginnOf1shock;
                    intersecInUpperPlane = true;
                    speed1 = S([ul;vl],a,sol2);
                    speed2 = lambda2((1/2)*(intersecPoint(1) + ur), (1/2)*(intersecPoint(2)+vr),a);
                    speedConsec = true;
                    type1 = 'S1';
                    type2 = 'S2';
                    statesWave1 = [Ul,intersecPoint];
                    statesWave2 = [intersecPoint,Ur];    
                elseif abs(targetFun(muOfInt))<tolerance
                    sol2 = muOfInt;
                    intersecPoint = intersecRarPoint;
                    intersecInUpperPlane = true;
                    speed1 = S([ul;vl],a,sol2);
                    speed2 = lambda2((1/2)*(intersecPoint(1) + ur), (1/2)*(intersecPoint(2)+vr),a);
                    speedConsec = true;
                    type1 = 'S1';
                    type2 = 'S2';
                    statesWave1 = [Ul,intersecPoint];
                    statesWave2 = [intersecPoint,Ur];    
                elseif targetFun(muEE)*targetFun(muOfInt)<0 % if true, there is a change of signs for the fzero function

                    sol2 = fzero(targetFun,[muEE,muOfInt]);
                    intersecPoint = [U([ul;vl],a,sol2);V([ul;vl],a,sol2)];
                    intersecInUpperPlane = true;

                    speed1 = S([ul;vl],a,sol2);
                    speed2 = lambda2((1/2)*(intersecPoint(1) + ur), (1/2)*(intersecPoint(2)+vr),a);
                    if speed1 > speed2 + tolerance
                        speedConsec = false;
                        type1 = 'noType';
                        type2 = 'noType';
                        statesWave1 = NaN;
                        statesWave2 = NaN;
                    else
                        speedConsec = true;
                        type1 = 'S1';
                        type2 = 'S2';
                        statesWave1 = [Ul,intersecPoint];
                        statesWave2 = [intersecPoint,Ur];                         
                    end

                 
                else
                    %we have to find the intersection with the mixed curve and the 2 shock curve
                     % findIntersec of shock Curve and mixedCurve
                    targetFun = @(muS) implicitHugoniotFun([ur;vr],mixedCurveByBarV([ul,vl],a,[ul;vl],uAxintersec,muS),a);
                    
                    if targetFun(0)*targetFun(vl)<0
                        sol2 = fzero(targetFun,[vl, 0]); 
                    elseif abs(targetFun(vl))<tolerance
                        sol2 = vl;
                    elseif abs(targetFun(0))<tolerance
                        sol2 = 0;
                    else
%                         targetFun(vl)
%                         targetFun(0)
%                         Ul
%                         Ur
                        'Weird2'
                    end

                        intersecPoint = mixedCurveByBarV([ul,vl],a,[ul;vl],uAxintersec,sol2);
                    intersecInUpperPlane = true;

                    %find the value of \bar{u}
                    uValue = uBarToVBar(a,[ul;vl],uAxintersec,sol2);

                    ptOnRare = [uValue;sol2];
                    speed1 = lambda1(ptOnRare(1),ptOnRare(2),a);        
                    speed2 = lambda2((1/2)*(intersecPoint(1) + ur), (1/2)*(intersecPoint(2)+vr),a);
                    if speed1 > speed2 + tolerance
                        speedConsec = false;
                        type1 = 'noType';
                        type2 = 'noType';
                        statesWave1 = NaN;
                        statesWave2 = NaN;    
                    else
                        speedConsec = true;
                        type1 = 'RS1';
                        type2 = 'S2';
                        statesWave1 = [Ul,ptOnRare,intersecPoint];
                        statesWave2 = [intersecPoint,Ur]; 
                        speed1 = [lambda1(ul,vl,a),speed1];
                    end
                end
            end


    elseif (RhoVal2r < Rho2_beginnOf1shock) % (ur,vr) lies to the right of the rarefaction curve through te point where the 1-shock starts
            
        %first find intersec with rarefaction and mixedCurve
        targetFun = @(muS) Rho2(mixedCurveByBarV([ul,vl],a,[ul;vl],uAxintersec,muS),a) - RhoVal2r;
        VvalBar=fzero(targetFun,[vl,0]);%computes \bar{v} according to intersection
        intersecRarPoint = mixedCurveByBarV([ul,vl],a,[ul;vl],uAxintersec,VvalBar);

        if intersecRarPoint(2)<=vr
            intersecPoint = intersecRarPoint;
            intersecInUpperPlane = true;
            speedConsec = true;
            type1 = 'RS1';
            type2 = 'R2';

          %uncomment if the speeds are needes as an output
            %find the point on the rarefaction 
            %find the value of \bar{u}
            
            uValue = uBarToVBar(a,[ul;vl],uAxintersec,VvalBar);

            ptOnRare = [uValue;VvalBar];
            speed1 = [lambda1(ul,vl,a),lambda1(ptOnRare(1),ptOnRare(2),a)];
            speed2 = [lambda2(intersecRarPoint(1),intersecRarPoint(2),a),lambda2(ur,vr,a)];
            statesWave1 = [Ul,ptOnRare,intersecPoint];
            statesWave2 = [intersecPoint,Ur]; 


        else
            % findIntersec of shock Curve and mixedCurve
            targetFun = @(muS) implicitHugoniotFun([ur,vr],mixedCurveByBarV([ul,vl],a,[ul;vl],uAxintersec,muS),a);

            if targetFun(0)*targetFun(VvalBar)<0
                sol2 = fzero(targetFun,[VvalBar,0]);
            elseif abs(targetFun(0))<10e-10
                sol2 = 0;
            elseif abs(targetFun(VvalBar))<10e-10
                sol2 = VvalBar;
            else
                'Weird'
            end
            intersecPoint = mixedCurveByBarV([ul,vl],a,[ul;vl],uAxintersec,sol2);
            intersecInUpperPlane = true;

            %find the value of \bar{u}
            uValue = uBarToVBar(a,[ul;vl],uAxintersec,sol2);

            ptOnRare = [uValue;sol2];
            speed1 = lambda1(ptOnRare(1),ptOnRare(2),a);        
            speed2 = lambda2((1/2)*(intersecPoint(1) + ur), (1/2)*(intersecPoint(2)+vr),a);
            if speed1 > speed2 + tolerance
                speedConsec = false;
                type1 = 'noType';
                type2 = 'noType';
                statesWave1 = NaN;
                statesWave2 = NaN;   
            else
                speedConsec = true;
                type1 = 'RS1';
                type2 = 'S2';
                statesWave1 = [Ul,ptOnRare,intersecPoint];
                statesWave2 = [intersecPoint,Ur]; 
                speed1 = [lambda1(ul,vl,a),speed1];
            end
        end

    end


end


end