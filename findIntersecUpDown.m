function [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersecUpDown(Ul,Ur,a)

%find the intersection in the case vl<0 vr>0


[intersecPoint,intersecInUpperPlane,speedConsec,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersectionInUpperHalfplane(Ul,Ur,a);

 if (intersecInUpperPlane == false) || (speedConsec == false)

    [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2]= findIntersectionInLowerHalfplane(Ul,Ur,a);
 end
%type1
%type2
end

