function [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersectionInLowerHalfplane(Ul,Ur,a)
% find intersection in upper halfplane. Uses only values vr>0 vl<0

%just like in the upper plane,just for negative values, use symmetry

ul = -Ur(1);
vl = -Ur(2);
ur = -Ul(1);
vr = -Ul(2);


[intersecPoint,~,~,type1_h,type2_h,statesWave1_h,statesWave2_h,speed1_h,speed2_h] = findIntersectionInUpperHalfplane([ul;vl],[ur;vr],a);


%set values in order that is works with symmetrie
intersecPoint = -intersecPoint;



statesWave2 = -fliplr(statesWave1_h);
statesWave1 = -fliplr(statesWave2_h);

speed2 = -fliplr(speed1_h);
speed1 = -fliplr(speed2_h);


if strcmp(type1_h,'S1')
    type2 = 'S2';
elseif strcmp(type1_h,'RS1')
    type2 = 'SR2';
end
if strcmp(type2_h,'S2')
    type1 = 'S1';
elseif strcmp(type2_h,'R2')
    type1 = 'R1';
else
    [ul;vl]
    [ur;vr]
    type2_h
    type1_h
    error('??')
end

end