% Parameter file: Choose between `Case` 1-7

% Main contributors: Anna-Simone Frank, Kamila Larripa, Hwayeon Ruy and Susanna
% Roeblitz

function params=parameters(Case)
 
% Baseline parameters (together with Case =1):

    params.pT=0.023;    %1/day; proliferation of tumor cells
    params.K=0.5;	%max relative amount of macrophages compared to tumor cells
    params.f=0.5;     %relative required tumor size to supports M2 polarization
    params.dT=0.01;       %1/day; macrophage independent tumor cell death rate
    params.r1=10;     %1/cells; contribution of M1 to the killing of tumor cells
    params.r2=10;     %1/cells; contribution of M2 to the proliferation of tumor cells
    params.p0=0.7;     %1/day; proliferation rate of M0
    params.g=1;    %1;        %contribution of M1 to recruitment of M0
    params.d0=0.1;       %1/day; M0 natural death rate
    params.alpha01=1;   %1/day; differentiation of M0 to M1
    params.alpha02=0.1;   %1/day; differentiation of M0 to M2
    params.alpha10=0.0001;%1/day; differentiation of M1 to M0
    params.alpha20=0.0001;%1/day; differentiation of M2 to M0
    params.b=2;     %contribution of tumor to M2 polarization
    params.d1=0.1;       %1/day; M1 natural death rate
    params.alpha1m=0.001;  %1/day; polarization of M1 to Mm
    params.d2=0.1;       %1/day; M2 natural death rate
    params.alpham2=0.01; %1/day; polarization of Mm to M2
    params.dm=0.1;       %1/day; Mm natural death rate

    %assumption: from M2 to M1 only over M0, takes long; M1 to M2 can go
    %over Mm or M0
    params.alpha12=0;
    params.alpha21=0;
    params.alpham1=0;
    params.alpha2m=0;

    if Case==1  %low
        params.alpham2=0.01;
        params.alpha02=0.1;
        params.alpha1m=0.001;
        params.alpha01=1;
    end

    % Specific cases to distuinguish between monot and bi-stability:

    if Case==6  %high
        params.alpham2=0.01;
        params.alpha02=1;
        params.alpha1m=0.001;
        params.alpha01=1;
    end

    if Case==3  %original bistable (medium/low)
        params.alpham2=0.01;
        params.alpha02=1;
        params.alpha1m=0.1;
        params.alpha01=10;
        
    end

    if Case==2  %bistable
        params.alpham2=0.01;
        params.alpha02=0.1;
        params.alpha1m=0.1; %increased
        params.alpha01=1;
    end

    if Case==4  %bistable 
        params.alpham2=0.01;
        params.alpha02=1;
        params.alpha1m=0.001;
        params.alpha01=1;
        params.alpha20=0.1; %changed
    end

    if Case==5  %bistable
        params.alpham2=0.01;
        params.alpha02=0.5;
        params.alpha1m=0.001;
        params.alpha01=1;
    end

    if Case==7  %high
        params.alpham2=0.01;
        params.alpha02=0.1;
        params.alpha1m=0.5; %increased
        params.alpha01=1;
    end


end

