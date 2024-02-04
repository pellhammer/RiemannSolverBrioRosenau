function Y = mixedCurveByBarV(Ul,a,barStartPt,uAxintersec,vSS)
% this function calculates the mixed curve parametrized by the u component
% on the rarefaction curve denotetn \bar{u}
%
% barStarPt :  this is the first point on the rarefaction curve in the
%   lower halfplane that is associated to a mixed curve point. For (ul,vl) in
%   the regions A1 and A2 this is alpha_f. For  (ul,vl) in
%   the region A3 is is (ul,vl). Is an input to the function such that it is
%   not computed several times
% uAxintersec :  this is the u value where the rarefaction curve in the
%   lower halfplane hits the u axis
% 
%
%
%

% vSS
% barStartPt

ul = Ul(1);
vl = Ul(2);


% %calculate alpha, is maybe absolute if alphf gets to be a infut to the
% %function
%  intersecWithSqA = @(u) Rho1([u;sqrt(a)*u],a) - Rho1([ul;vl],a);
%             alpha_fu = fzero(intersecWithSqA,-vl*(1/sqrt(a)));
%             alpha_f = [alpha_fu;-sqrt(a)*alpha_fu];
% 
%     % calculate the intersection of the rarefaction curve with the u axis   
%     targetFun = @(muS) Rho1(alpha_f,a)-Rho1([muS;0],a);
%     uAxintersecParam = fzero(targetFun,0);

    %find the value of \bar{u}
    targetFun = @(uS) Rho1([ul;vl],a) - Rho1([uS;vSS],a);

    % vValue = fzero(targetFun,[alpha_f(2),0]);
    uValue = fzero(targetFun,-1);
    ptOnRare = [uValue;vSS];

if vSS/uValue<=-sqrt(a)+10e-8

    Y = -barStartPt;%stimmt nur for barstarpt=alpha_f, fall sollte aber nie auftreten

elseif  vSS/uValue>=-10e-7

    Y = [uAxintersec;0];
else

     Y = [U(ptOnRare,a,muE(ptOnRare,a));V(ptOnRare,a,muE(ptOnRare,a))];
end
end

