syms T M0 M1 M2 Mm pT r1 r2 dT p0 f g d0 alpha01 alpha02 b alpha10 alpha20 d1 alpha alpha21 alpha12 alpham1 alpha1m ...
    d2 alpha2m alpham2 dm K
%parameters;

%Old functions
 % f1=pT*T*(1-T)*(0.1-r1*M1+r2*M2)-dT*T;
 % f2=p0*T*(1-(M0+M1+M2+Mm)/f)*(1+g*M1)-d0*M0-alpha01*M0-...
 %        alpha02*M0*(1+b*T/(T+f))+alpha10*M1+alpha20*M2;
 % 
 %  f3=-d1*M1+alpha01*M0-alpha10*M1+alpha21*M2-alpha12*M1*(1+b*T/(T+f))+...
 %        alpham1*Mm-alpha1m*M1*(1+b*T/(T+f));
 % 
 %  f4=-d2*M2-alpha20*M2-alpha21*M2-alpha2m*M2+alpha02*M0*(1+b*T/(T+f))+...
 %        +alpha12*M1*(1+b*T/(T+f))+alpham2*Mm*(1+b*T/(T+f));
 % 
 %  f5=-dm*Mm-alpham1*Mm-alpham2*Mm*(1+b*T/(T+f))+...
 %        alpha1m*M1*(1+b*T/(T+f))+alpha2m*M2;
 % 

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