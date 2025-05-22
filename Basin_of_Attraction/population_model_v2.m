% Function solves the ODE system, and plots time series trajectories
% It also calculates steasy states and their stability

% Main contributor: Susanna Roblitz

function population_model()

    
    disp('Which parameter case do you want to run?')
    Case = input('Choose between 1 (low), 2 (bistable medium/low), 3 (high): ');
    params=parameters(Case) ;

    %define initial conditions:

    %high initial tumor state T0
    y0_high = [0.9, 0.1,0.1,0.1,0.1]; 
    
   

    %low T0
    y0_low = [0.09, 0.1,0.1,0.1,0.1]; 
    

    tspan=[0,500];

    y0=y0_high;  %change here to y0_high/y0_low if needed  

    [t,y]=ode23s(@(t,y) rhs(t,y,params),tspan,y0);

    % Figure initation: 
    fig=figure(1);
    ax=axes;
    %---
    yyaxis right
    hT=plot(t,y(:,1),'k-','LineWidth',2);
    ylim([0,1])
    ylabel('tumor size')
    xlabel('time')
    hold off
   %---
    yyaxis left
    h0=plot(t,y(:,2),'c','LineWidth',2);
    hold on
    h1=plot(t,y(:,3),'r--','LineWidth',2);
    h2=plot(t,y(:,4),'b:','LineWidth',2);
    hM=plot(t,y(:,5),'m-.','LineWidth',2);
    ylim([0,0.4])
    ylabel('cell populations')
    hold off
    xlim(tspan)
    legend([hT,h0,h1,h2,hM],'Tumor','M0','M1','M2','Mm')
    fontsize(fig,18,'points')
    
    % Set the color of each axis to black
    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(2).Color = [0 0 0];


    % Stability analysis:

    %finding steady states as roots of the rhs
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
    %---
    x0 = y0_low;
    if Case==1
        x0 = y0_low;
    end
    if Case==3
        x0 = y0_high;
    end
    if Case==2
        x0 = y0_high;
    end
    % Solving for steady states and calculating stability for Case =2:
    x=fsolve(@(y) rhs(0,y,params),x0)
    J=jacob(0,x,params);
    la=eig(J)
    if Case==2
        x0(1)=1-x0(1);
        x=fsolve(@(y) rhs(0,y,params),x0)
        J=jacob(0,x,params);
        la=eig(J)
    end

end

function [dy,J]=rhs(t,y,params)


    T=y(1);     %tumor cell
    M0=y(2);    %M0 macrophage
    M1=y(3);    %M1 macrophage
    M2=y(4);    %M2 macophage
    Mm=y(5);    %M12 mixed phenotype macrophage

    

    J=jacob(t,y,params);

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

end

function Jac=jacob(~,y,params)

    T=y(1);     %tumor cell
    M0=y(2);    %M0 macrophage
    M1=y(3);    %M1 macrophage
    M2=y(4);    %M2 macophage
    Mm=y(5);    %M12 mixed phenotype macrophage

    Jac=zeros(5,5);
    Jac(1,1)=params.pT*(1-2*T)*(1-params.r1*M1+params.r2*M2)-params.dT;
    Jac(1,2)=0;
    Jac(1,3)=-params.pT*T*(1-T)*params.r1;
    Jac(1,4)=params.pT*T*(1-T)*params.r2;
    Jac(1,5)=0;
    Jac(2,1)=params.p0*(1-(M0+M1+M2+Mm)/params.K)*(1+params.g*M1)-params.alpha02*M0*params.b*params.f/(T+params.f)^2;
    Jac(2,2)=-params.p0*T*(1+params.g*M1)/params.K-params.d0-params.alpha01-params.alpha02*(1+params.b*T/(T+params.f));
    Jac(2,3)=-params.p0*T*(1+params.g*M1)/params.K+params.p0*T*(1-(M0+M1+M2+Mm)/params.K)*params.g+params.alpha10;
    Jac(2,4)=-params.p0*T*(1+params.g*M1)/params.K+params.alpha20;
    Jac(2,5)=-params.p0*T*(1+params.g*M1)/params.K;
    Jac(3,1)=-params.alpha12*M1*params.b*params.f/(T+params.f)^2-params.alpha1m*M1*params.b*params.f/(T+params.f)^2;
    Jac(3,2)=params.alpha01;
    Jac(3,3)=-params.d1-params.alpha10-params.alpha12*(1+params.b*T/(T+params.f))-params.alpha1m*(1+params.b*T/(T+params.f));
    Jac(3,4)=params.alpha21;
    Jac(3,5)=params.alpham1;
    Jac(4,1)=params.alpha02*M0*params.b*params.f/(T+params.f)^2+params.alpha12*M1*params.b*params.f/(T+params.f)^2+params.alpham2*Mm*params.b*params.f/(T+params.f)^2;
    Jac(4,2)=params.alpha02*(1+params.b*T/(T+params.f));
    Jac(4,3)=params.alpha12*(1+params.b*T/(T+params.f));
    Jac(4,4)=-params.d2-params.alpha20-params.alpha21-params.alpha2m;
    Jac(4,5)=params.alpham2*(1+params.b*T/(T+params.f));
    Jac(5,1)=-params.alpham2*Mm*params.b*params.f/(T+params.f)^2+params.alpha1m*M1*params.b*params.f/(T+params.f)^2;
    Jac(5,2)=0;
    Jac(5,3)=params.alpha1m*(1+params.b*T/(T+params.f));
    Jac(5,4)=params.alpha2m;
    Jac(5,5)=-params.dm-params.alpham1-params.alpham2*(1+params.b*T/(T+params.f));

end