%function [dy,J]=odefun(t,y,b,params,Case,bif_param)
%Input;

function [dy,J]=odefun(t,y,b,params,Case,bif_param)
%function [dy,J]=odefun(t,y,f,params,Case)
    T=y(1);     %tumor cell
    M0=y(2);    %M0 macrophage
    M1=y(3);    %M1 macrophage
    M2=y(4);    %M2 macophage
    Mm=y(5);    %M12 mixed phenotype macrophage

    %bif_param
 % Handle bifurcation parameter selection dynamically
    switch bif_param
        case 'b'
            params.b = b;  % !!! Update the 'b' parameter
        case 'r1'
            params.r1 = b;  % !!! Update the 'r1' parameter
        case 'r2'
            params.r2 = b;  % !!! Update the 'r2' parameter
        case 'f'
            params.f = b;  % Update the 'f' parameter
        case 'K'
            params.K = b;  % Update the 'K' parameter
        case 'p0'
            params.p0 = b;
        case 'pT'
            params.pT = b; %!!!
        case 'dT'
            params.dT = b; % !!!
        case 'd2'
            params.d2 = b;%!!!
        case 'alpha21'
            params.alpha21 = b; % 
        otherwise
            disp('No such bifurcation parameter!');
    end


  % params.alpha01=alpha01;
    %params.alpha02=alpha02;
    % params.alpha12=alpha12;
      % params.alpha21=alpha21;
        %params.alpha1m=alpha1m;
        %params.alpham2=alpham2;
      %  params.alpham1=alpham1;
      %  params.alpha2m=alpha2m;
     % params.p0=p0;
%==========================================================================
% New functions
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
%==========================================================================
% Old functions
% dy=zeros(5,1);
%     %T
%     dy(1)=params.pT*T*(1-T)*(params.h-params.r1*M1+params.r2*M2)-params.dT*T;
%     %dy(1)=pT*T*(1-T)*(M2^n/(M2^n+M1^n)-0.5)-dT*T;
%     %dy(1)=pT*M2^n/(M2^n+M1^n)-dT*T;
%     %dy(1)=pT*T*(1-T)*(h+r2*M2)-dT*T*(h+r1*M1);
%     %M0
%     %dy(2)=p0*T*(1+g*M1)-d0*M0-alpha01*M0-...
%     %dy(2)=p0*T*M0*(1-(M0+M1+M2+Mm)/f)*(1+g*M1)-d0*M0-alpha01*M0-...
%     %dy(2)=p0*T*M0*(1-M0/f)*(1+g*M1)-d0*M0-alpha01*M0-...
%     %dy(2)=p0*T*(1+g*M1)-d0*M0-alpha01*M0-...
%     dy(2)=params.p0*T*(1-(M0+M1+M2+Mm)/params.K)*(1+params.g*M1)+...
%         +params.alpha10*M1+params.alpha20*M2-...
%         params.d0*M0-params.alpha01*M0-...
%         params.alpha02*M0*(1+params.b*T/(T+params.f));
%     %M1
%     dy(3)=-params.d1*M1-params.alpha10*M1+...
%         params.alpha01*M0+params.alpha21*M2+params.alpham1*Mm-...
%         (params.alpha12*M1+params.alpha1m*M1)*(1+params.b*T/(T+params.f));
%     %M2
%     dy(4)=-params.d2*M2-params.alpha20*M2-params.alpha21*M2-params.alpha2m*M2+...
%         (params.alpha02*M0+params.alpha12*M1+params.alpham2*Mm)*(1+params.b*T/(T+params.f));
% 
%     %Mm
%     dy(5)=-params.dm*Mm-params.alpham1*Mm+params.alpha2m*M2+...
%         (params.alpha1m*M1-params.alpham2*Mm)*(1+params.b*T/(T+params.f));
%  dydt = [dy(1); dy(2); dy(3); dy(4); dy(5)];