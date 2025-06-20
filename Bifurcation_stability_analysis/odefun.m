% ODE system of tumor-macrophage interactions in the TME for bifurcation analysis

% Main contributor: Anna-Simone Frank
%-------------------------------------------------------------------------

function [dy,J]=odefun(t,y,b,params,Case,bif_param)
    T=y(1);     %tumor cell
    M0=y(2);    %M0 macrophage
    M1=y(3);    %M1 macrophage
    M2=y(4);    %M2 macophage
    Mm=y(5);    %M12 mixed phenotype macrophage

    
 % Handle bifurcation parameter selection dynamically
    switch bif_param
        case 'b'
            params.b = b;  
        case 'r1'
            params.r1 = b;  
        case 'r2'
            params.r2 = b;  
        case 'f'
            params.f = b;  
        case 'K'
            params.K = b;  
        case 'p0'
            params.p0 = b;
        case 'pT'
            params.pT = b; 
        case 'dT'
            params.dT = b; 
        case 'd2'
            params.d2 = b;
        case 'alpha21'
            params.alpha21 = b;  
        otherwise
            disp('No such bifurcation parameter!');
    end

%==========================================================================
% ODE functions
dy=zeros(5,1);
    %T
    dy(1)=params.pT*T*(1-T)*(1-params.r1*M1+params.r2*M2)-params.dT*T;
    % M0
    dy(2)=params.p0*T*(1-(M0+M1+M2+Mm)/params.K)*(1+params.g*M1)+...
        +params.alpha10*M1+params.alpha20*M2-...
        params.d0*M0-params.alpha01*M0-...
        params.alpha02*M0*(1+params.b*T/(T+params.f));
    %M1
    dy(3)=-params.d1*M1-params.alpha10*M1+...
        params.alpha01*M0+params.alpha21*M2+params.alpham1*Mm-...
        (params.alpha12*M1+params.alpha1m*M1)*(1+params.b*T/(T+params.f));
    %M2
    dy(4)=-params.d2*M2-params.alpha20*M2-params.alpha21*M2-params.alpha2m*M2+...
        (params.alpha02*M0+params.alpha12*M1+params.alpham2*Mm)*(1+params.b*T/(T+params.f));

    %Mm
    dy(5)=-params.dm*Mm-params.alpham1*Mm+params.alpha2m*M2+...
        (params.alpha1m*M1-params.alpham2*Mm)*(1+params.b*T/(T+params.f));


 dydt = [dy(1); dy(2); dy(3); dy(4); dy(5)];
end
