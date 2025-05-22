% Code simply calculates the jacobian matrix of the tumor-macrophage model
% Jacobian matrix is needed to determine eigenvalues and the satbility of steady states

% Main contributor: Anna-Simone Frank 

%--------------------------------------------------------------------------

syms T M0 M1 M2 Mm pT r1 r2 dT p0 f g d0 alpha01 alpha02 b alpha10 alpha20 d1 alpha alpha21 alpha12 alpham1 alpha1m ...
    d2 alpha2m alpham2 dm K

 %New function
  f1=pT*T*(1-T)*(1-r1*M1+r2*M2)-dT*T;
    % M0
    f2=p0*T*(1-(M0+M1+M2+Mm)/K)*(1+g*M1)+alpha10*M1+alpha20*M2-d0*M0-alpha01*M0-alpha02*M0*(1+b*T/(T+f));
    %M1
    f3=-d1*M1-alpha10*M1+...
        alpha01*M0+alpha21*M2+alpham1*Mm-...
        (alpha12*M1+alpha1m*M1)*(1+b*T/(T+f));
    %M2
    f4=-d2*M2-alpha20*M2-alpha21*M2-alpha2m*M2+...
        (alpha02*M0+alpha12*M1+alpham2*Mm)*(1+b*T/(T+f));

    %Mm
    f5=-dm*Mm-alpham1*Mm+alpha2m*M2+...
        (alpha1m*M1-alpham2*Mm)*(1+b*T/(T+f));

    jacobian([f1,f2,f3,f4,f5],[T,M0,M1,M2,Mm])