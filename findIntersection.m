function [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersection(Ul,Ur,a)


ul = Ul(1);
vl = Ul(2);
ur = Ur(1);
vr = Ur(2);

% set initial values to expolid symmetrie
tolerance = 1e-10;


if abs(ul-ur)<tolerance && abs(vl-vr)<tolerance

    intersecPoint = Ul;
    type1 = 'S1';
    type2 = 'S2';
    statesWave1 = [Ul,Ul];
    statesWave2 = [Ul,Ul];
    speed1 = lambda1((1/2)*(ul+ur),(1/2)*(vl+vr),a);
    speed2 = lambda2((1/2)*(ul+ur),(1/2)*(vl+vr),a);


elseif vl<-tolerance && vr>tolerance % updownClassic

    [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersecUpDown([ul;vl],[ur;vr],a);

elseif vr<-tolerance && vl>tolerance % updownMirror

    [intersecPoint_h,type1,type2,statesWave1_h,statesWave2_h,speed1,speed2] = findIntersecUpDown([ul;-vl],[ur;-vr],a);
    intersecPoint = [intersecPoint_h(1);-intersecPoint_h(2)];
    statesWave1 = [statesWave1_h(1,:);-statesWave1_h(2,:)];
    statesWave2 = [statesWave2_h(1,:);-statesWave2_h(2,:)];

elseif vl<-tolerance && vr<-tolerance % SmollerJohnsonClassic


    [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersecSmoJoh([ul;vl],[ur;vr],a);
%     smojo2=intersecPoint

elseif vl>tolerance && vr>tolerance % SmollerJohnsonClassicMirror


    [intersecPoint_h,type1,type2,statesWave1_h,statesWave2_h,speed1,speed2] = findIntersecSmoJoh([ul;-vl],[ur;-vr],a);
    intersecPoint = [intersecPoint_h(1);-intersecPoint_h(2)];
    statesWave1 = [statesWave1_h(1,:);-statesWave1_h(2,:)];
    statesWave2 = [statesWave2_h(1,:);-statesWave2_h(2,:)];


elseif abs(vl)<tolerance && vr>tolerance %vl =0 classic
    vl = 0;
    [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersecvLNull([ul;vl],[ur;vr],a);

elseif abs(vl)<tolerance && vr<-tolerance %vl =0 classic mirror
    vl = 0;
    [intersecPoint_h,type1,type2,statesWave1_h,statesWave2_h,speed1,speed2] = findIntersecvLNull([ul;-vl],[ur;-vr],a);
    intersecPoint = [intersecPoint_h(1);-intersecPoint_h(2)];
    statesWave1 = [statesWave1_h(1,:);-statesWave1_h(2,:)];
    statesWave2 = [statesWave2_h(1,:);-statesWave2_h(2,:)];


elseif abs(vr)<tolerance && vl<-tolerance

    vr = 0;
    [intersecPoint_h,type1_h,type2_h,statesWave1_h,statesWave2_h,speed1_h,speed2_h] = findIntersecvLNull(-[ur;vr],-[ul;vl],a);


    intersecPoint = -intersecPoint_h;
    statesWave2 = -fliplr(statesWave1_h);
    statesWave1 = -fliplr(statesWave2_h);

    speed2 = -fliplr(speed1_h);
    speed1 = -fliplr(speed2_h);

    if strcmp(type1_h,'S1')
      type2 = 'S2';
    elseif strcmp(type1_h,'R1')
        type2 = 'R2';
    end
    if strcmp(type2_h,'S2')
        type1 = 'S1';
    elseif strcmp(type2_h,'R2')
        type1 = 'R1';
    end

    

elseif abs(vr)<tolerance && vl>tolerance

    vr = 0;
    [intersecPoint_h,type1_h,type2_h,statesWave1_h,statesWave2_h,speed1_h,speed2_h] = findIntersecvLNull(-[ur;-vr],-[ul;-vl],a); %Note that we have to exploit both symmetries!!!!!!

    
    %first invert the first symmetrie
    intersecPoint_h2 = -intersecPoint_h;


    statesWave2_h2 = -fliplr(statesWave1_h);
    statesWave1_h2 = -fliplr(statesWave2_h);

    speed2 = -fliplr(speed1_h);
    speed1 = -fliplr(speed2_h);
    

    if strcmp(type1_h,'S1')
      type2 = 'S2';
    elseif strcmp(type1_h,'R1')
        type2 = 'R2';
    end
    if strcmp(type2_h,'S2')
        type1 = 'S1';
    elseif strcmp(type2_h,'R2')
        type1 = 'R1';
    end

    %undo second symmetry
    intersecPoint = [intersecPoint_h2(1);-intersecPoint_h2(2)];
    statesWave1 = [statesWave1_h2(1,:);-statesWave1_h2(2,:)];
    statesWave2 = [statesWave2_h2(1,:);-statesWave2_h2(2,:)];




elseif abs(vr)<tolerance && abs(vl)<tolerance
     vr = 0;
     vl = 0;


    [intersecPoint,type1,type2,statesWave1,statesWave2,speed1,speed2] = findIntersecvLvRNull([ul;vl],[ur;vr],a);



end

end


