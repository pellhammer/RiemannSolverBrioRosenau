function Y = interalCurve1(Ul,Ur,a,uAxintersec,tau)
% This function paramterizes the rarefactioncurve of Ul by its eigenvalue
% lambda1 

ul = Ul(1);
vl = Ul(2);
ur = Ur(1);
vr = Ur(2);

tolerance = 1e-14;



    if abs(vl) > tolerance || (abs(vl) < tolerance && ul>0)
        % function lambda1=tau paramterizen by u 
        if vl<=0 && vr<=0
            VV = @(u) -((tau - 2*u).*(tau - 2*a*u)).^(1/2)./2;
        elseif vl>=0 && vr>=0 
            VV = @(u) ((tau - 2*u).*(tau - 2*a*u)).^(1/2)./2;
        else
            Ul
            Ur
        end

        % calculate the intersection of the rarefaction curve with the u
        % axis---später in functioncall direkt übergeben
%         targetFun = @(muS) Rho1([ul;vl],a)-Rho1([muS;0],a);
%         uAxintersec = fzero(targetFun,0);
        Lambda1Max = lambda1(uAxintersec,0,a);

        if tau < Lambda1Max
            % Value where VV hits  the u axis:
            KK = (ul>=0) *(tau/2) + (ul<0) *(tau/(2*a));
    
            targetFun = @(uS) Rho1(Ul,a) - Rho1([uS;VV(uS)],a);
            
            if targetFun(KK)*targetFun(uAxintersec)<0
                sol = fzero(targetFun,[KK,uAxintersec]);
            elseif abs(targetFun(KK))<tolerance
                sol = KK;
            elseif abs(targetFun(uAxintersec))<tolerance
                sol = uAxintersec;
            else
                error('Something here went wrong');
            end
    
            Y = [sol;VV(sol)];
        
        else
            Y = [uAxintersec;0];
        end


    elseif abs(vl) < tolerance

        if ul<0

            Y = [tau/(2*a);0];

        end
    end

end



















