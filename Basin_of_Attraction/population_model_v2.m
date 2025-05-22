% function [t,y1,y2]=population_model_v2
% Case=2
% params=parameters(Case);
%     %fields = fieldnames(params)
%  %for  i=1:1:numel(fields)
% 
% %     fields = fieldnames(params.b)
% %  for  i=1:1:numel(fields)
% % params.b=params.b.(fields{i})
% 
% 
%    % y0=[0.9,0.001,0.001,0.001,0.001];
% y0=[0.9,1,1,1,1];
% 
%     %y0=[0.003,0.0003,0.001,0.001,0];
%     %y0=[0.0009    0.0001    0.0057    0.0006    0.0001];
%     %y0=[0.0009    0.0001    0.0001    0.06    0.0001];
% 
%     tspan=[0,1000];
% 
% 
%     [t,y1]=ode23s(@(t,y) rhs(t,y,params),tspan,y0);
% 
%    % fields = fieldnames(params.b)
%     % save steady state of tumor
%     %Tss=zeros(numel(fields),1);
%   % Tss=zeros(7,1);
%   %for  i=1:1:numel(fields)
%      % params.b.(fields{i})
%        %params.b=params.b.b2
%          %[t,y]=ode23s(@(t,y) rhs(t,y,params),tspan,y0);
%     % Tss= [, y(end,1)];
% 
%  %end
% 
% 
%     figure(1)
%     plot(t,y1(:,1),'b-*','LineWidth', 2,'MarkerSize',3)
%     xlabel('Time t');
%     ylabel('Cell population dynamics');
%     legend('tumor size','Location','NorthWest')
%     hold on
%     plot(t,y1(:,2),'k','LineWidth', 2,'MarkerSize',3)
%     hold on
%     plot(t,y1(:,3),'r-','LineWidth', 2,'MarkerSize',3)
%     plot(t,y1(:,4),'b-.','LineWidth', 2,'MarkerSize',3)
%     plot(t,y1(:,5),'m-.','LineWidth', 2,'MarkerSize',3)
%     legend('T','M0','M1','M2','Mm')
% hold all;
%    ylim([0 1]);
% 
% 
%     %figure(3)
%     %plot(t,y(:,2)+y(:,3)+y(:,4)+y(:,5))
%     %legend('sum of all macrophages')
% 
%     z0 = y0; %Making use of same inital conditions
%     z0(1)=1-y0(1);
%     %z0=[1-y0(1),0.1,0.1,0.1,0.01];
%     [t,y2]=ode23s(@(t,y) rhs(t,y,params),tspan,z0);
%    % figure(1)
%     %plot(t,y2(:,1),'r-*')
%      %xlabel('Time t');
%     %ylabel('Cell population dynamics');
%     %legend('tumor size','Location','NorthWest')
% 
%    % hold on;
% figure(2)
%   plot(t,y2(:,1),'b-*','LineWidth', 2,'MarkerSize',3)
%   hold on
% xlabel('Time t');
%     ylabel('Cell population dynamics');
%     legend('tumor size','Location','NorthWest')
%     plot(t,y2(:,2),'k','LineWidth', 2,'MarkerSize',3)
%     hold on
%     plot(t,y2(:,3),'r-','LineWidth', 2,'MarkerSize',3)
%     plot(t,y2(:,4),'b-.','LineWidth', 2,'MarkerSize',3)
%     plot(t,y2(:,5),'m-.','LineWidth', 2,'MarkerSize',3)
%     legend('T2','M0','M1','M2','Mm')
%     hold on
%   ylim([0 1])
% %     figure(4)
% %     for f=0.1:0.1:1
% %         f
% %         [t,y]=ode23s(@(t,y) rhs(t,y,f),tspan,y0);
% %         ss=sum(y(end,2:5));
% %         plot(f,ss,'*')
% %         hold on
% %     end
% 
%     %finding steady states as roots of the rhs
%     options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
%     %x0=[1,0.1,0.1,0.1,0.01];
%     if Case==1
%         x0=[0.01,0.006,0.06,0.006,0.0001];
%     end
%     if Case==3
%         x0=[0.9,0.006,0.06,0.006,0.0001];
%     end
%     if Case==2
%         x0=[0.9,0.01,0.01,0.01,0.01];% 0.9,0.006,0.06,0.006,0.0001
%     end
%      if Case==7
%         x0=[0.9,0.006,0.06,0.006,0.0001];
%     end
%     x=fsolve(@(y) rhs(0,y,params),x0) %Equilibrium for high T IC
%     J=jacob(0,x,params);
%     la=eig(J)
%     if Case==2
%         x0(1)=1-x0(1);
%         x2=fsolve(@(y) rhs(0,y,params),x0) % Equilibrium for low T IC
%         J=jacob(0,x,params);
%         la=eig(J)
%     end
% 
% end
% 
% function [dy,J]=rhs(t,y,params)
% 
% 
%     T=y(1);     %tumor cell
%     M0=y(2);    %M0 macrophage
%     M1=y(3);    %M1 macrophage
%     M2=y(4);    %M2 macophage
%     Mm=y(5);    %M12 mixed phenotype macrophage
% 
% 
% 
%     J=jacob(t,y,params);
% 
% % New functions
% dy=zeros(5,1);
%     %T
%     dy(1)=params.pT*T*(1-T)*(1-params.r1*M1+params.r2*M2)-params.dT*T;
%     % M0
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
% 
% 
%  %dydt = [dy(1); dy(2); dy(3); dy(4); dy(5)];
% 
%     % dy=zeros(5,1);
%     % %T
%     % dy(1)=params.pT*T*(1-T)*(params.h-params.r1*M1+params.r2*M2)-params.dT*T;
%     % %dy(1)=pT*T*(1-T)*(M2^n/(M2^n+M1^n)-0.5)-dT*T;
%     % %dy(1)=pT*M2^n/(M2^n+M1^n)-dT*T;
%     % %dy(1)=pT*T*(1-T)*(h+r2*M2)-dT*T*(h+r1*M1);
%     % %M0
%     % %dy(2)=p0*T*(1+g*M1)-d0*M0-alpha01*M0-...
%     % %dy(2)=p0*T*M0*(1-(M0+M1+M2+Mm)/f)*(1+g*M1)-d0*M0-alpha01*M0-...
%     % %dy(2)=p0*T*M0*(1-M0/f)*(1+g*M1)-d0*M0-alpha01*M0-...
%     % %dy(2)=p0*T*(1+g*M1)-d0*M0-alpha01*M0-...
%     % dy(2)=params.p0*T*(1-(M0+M1+M2+Mm)/params.K)*(1+params.g*M1)+...
%     %     +params.alpha10*M1+params.alpha20*M2-...
%     %     params.d0*M0-params.alpha01*M0-...
%     %     params.alpha02*M0*(1+params.b*T/(T+params.f));
%     % %M1
%     % dy(3)=-params.d1*M1-params.alpha10*M1+...
%     %     params.alpha01*M0+params.alpha21*M2+params.alpham1*Mm-...
%     %     (params.alpha12*M1+params.alpha1m*M1)*(1+params.b*T/(T+params.f));
%     % %M2
%     % dy(4)=-params.d2*M2-params.alpha20*M2-params.alpha21*M2-params.alpha2m*M2+...
%     %     (params.alpha02*M0+params.alpha12*M1+params.alpham2*Mm)*(1+params.b*T/(T+params.f));
%     % 
%     % %Mm
%     % dy(5)=-params.dm*Mm-params.alpham1*Mm+params.alpha2m*M2+...
%     %     (params.alpha1m*M1-params.alpham2*Mm)*(1+params.b*T/(T+params.f));
% 
% end
% 
% function Jac=jacob(t,y,params)
% 
%     T=y(1);     %tumor cell
%     M0=y(2);    %M0 macrophage
%     M1=y(3);    %M1 macrophage
%     M2=y(4);    %M2 macophage
%     Mm=y(5);    %M12 mixed phenotype macrophage
% 
%     %parameters
% 
%     % New functions
%       Jac=zeros(5,5);
%      Jac(1,1)=params.pT*(1-2*T)*(1-params.r1*M1+params.r2*M2)-params.dT;
%      Jac(1,2)=0;
%      Jac(1,3)=-params.pT*T*(1-T)*params.r1;
%      Jac(1,4)=params.pT*T*(1-T)*params.r2;
%      Jac(1,5)=0;
%      Jac(2,1)=params.p0*(1-(M0+M1+M2+Mm)/params.K)*(1+params.g*M1)-params.alpha02*M0*params.b*params.f/(T+params.f)^2;
%      Jac(2,2)=-params.p0*T*(1+params.g*M1)/params.K-params.d0-params.alpha01-params.alpha02*(1+params.b*T/(T+params.f));
%      Jac(2,3)=-params.p0*T*(1+params.g*M1)/params.K+params.p0*T*(1-(M0+M1+M2+Mm)/params.K)*params.g+params.alpha10;
%      Jac(2,4)=-params.p0*T*(1+params.g*M1)/params.K+params.alpha20;
%      Jac(2,5)=-params.p0*T*(1+params.g*M1)/params.K;
%      Jac(3,1)=-params.alpha12*M1*params.b*params.f/(T+params.f)^2-params.alpha1m*M1*params.b*params.f/(T+params.f)^2;
%      Jac(3,2)=params.alpha01;
%      Jac(3,3)=-params.d1-params.alpha10-params.alpha12*(1+params.b*T/(T+params.f))-params.alpha1m*(1+params.b*T/(T+params.f));
%      Jac(3,4)=params.alpha21;
%      Jac(3,5)=params.alpham1;
%      Jac(4,1)=params.alpha02*M0*params.b*params.f/(T+params.f)^2+params.alpha12*M1*params.b*params.f/(T+params.f)^2+params.alpham2*Mm*params.b*params.f/(T+params.f)^2;
%      Jac(4,2)=params.alpha02*(1+params.b*T/(T+params.f));
%      Jac(4,3)=params.alpha12*(1+params.b*T/(T+params.f));
%      Jac(4,4)=-params.d2-params.alpha20-params.alpha21-params.alpha2m;
%      Jac(4,5)=params.alpham2*(1+params.b*T/(T+params.f));
%      Jac(5,1)=-params.alpham2*Mm*params.b*params.f/(T+params.f)^2+params.alpha1m*M1*params.b*params.f/(T+params.f)^2;
%      Jac(5,2)=0;
%      Jac(5,3)=params.alpha1m*(1+params.b*T/(T+params.f));
%      Jac(5,4)=params.alpha2m;
%      Jac(5,5)=-params.dm-params.alpham1-params.alpham2*(1+params.b*T/(T+params.f));
% 
%      % Old functions
%     % Jac=zeros(5,5);
%     % Jac(1,1)=params.pT*(1-2*T)*(params.h-params.r1*M1+params.r2*M2)-params.dT;
%     % Jac(1,2)=0;
%     % Jac(1,3)=-params.pT*T*(1-T)*params.r1;
%     % Jac(1,4)=params.pT*T*(1-T)*params.r2;
%     % Jac(1,5)=0;
%     % Jac(2,1)=params.p0*(1-(M0+M1+M2+Mm)/params.K)*(1+params.g*M1)-params.alpha02*M0*params.b*params.f/(T+params.f)^2;
%     % Jac(2,2)=-params.p0*T*(1+params.g*M1)/params.K-params.d0-params.alpha01-params.alpha02*(1+params.b*T/(T+params.f));
%     % Jac(2,3)=-params.p0*T*(1+params.g*M1)/params.K+params.p0*T*(1-(M0+M1+M2+Mm)/params.K)*params.g+params.alpha10;
%     % Jac(2,4)=-params.p0*T*(1+params.g*M1)/params.K+params.alpha20;
%     % Jac(2,5)=-params.p0*T*(1+params.g*M1)/params.K;
%     % Jac(3,1)=-params.alpha12*M1*params.b*params.f/(T+params.f)^2-params.alpha1m*M1*params.b*params.f/(T+params.f)^2;
%     % Jac(3,2)=params.alpha01;
%     % Jac(3,3)=-params.d1-params.alpha10-params.alpha12*(1+params.b*T/(T+params.f))-params.alpha1m*(1+params.b*T/(T+params.f));
%     % Jac(3,4)=params.alpha21;
%     % Jac(3,5)=params.alpham1;
%     % Jac(4,1)=params.alpha02*M0*params.b*params.f/(T+params.f)^2+params.alpha12*M1*params.b*params.f/(T+params.f)^2+params.alpham2*Mm*params.b*params.f/(T+params.f)^2;
%     % Jac(4,2)=params.alpha02*(1+params.b*T/(T+params.f));
%     % Jac(4,3)=params.alpha12*(1+params.b*T/(T+params.f));
%     % Jac(4,4)=-params.d2-params.alpha20-params.alpha21-params.alpha2m;
%     % Jac(4,5)=params.alpham2*(1+params.b*T/(T+params.f));
%     % Jac(5,1)=-params.alpham2*Mm*params.b*params.f/(T+params.f)^2+params.alpha1m*M1*params.b*params.f/(T+params.f)^2;
%     % Jac(5,2)=0;
%     % Jac(5,3)=params.alpha1m*(1+params.b*T/(T+params.f));
%     % Jac(5,4)=params.alpha2m;
%     % Jac(5,5)=-params.dm-params.alpham1-params.alpham2*(1+params.b*T/(T+params.f));
% 
% end
function population_model()

    
    disp('Which parameter case do you want to run?')
    Case = input('Choose between 1 (low), 2 (bistable medium/low), 3 (high): ');
    params=parameters(Case) ;

    %define initial conditions:

    %high initial tumor state T0
    y0_high = [0.9, 0.1,0.1,0.1,0.1]; 
    %y0=[0.9,0.0003,0.001,0.001,0]; %original
   

    %low T0
    y0_low = [0.09, 0.1,0.1,0.1,0.1]; 
    %y0=[0.1,0.0003,0.001,0.001,0]; %original

    tspan=[0,500];

    y0=y0_high;  %change here to y0_high/y0_low if needed  

    [t,y]=ode23s(@(t,y) rhs(t,y,params),tspan,y0);
    fig=figure(1);
    ax=axes;
    yyaxis right
    hT=plot(t,y(:,1),'k-','LineWidth',2);
    ylim([0,1])
    ylabel('tumor size')
    xlabel('time')
    %legend('T')
    hold off
    %figure(2)
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
    %set(fig,'DefaultLineLineWidth',3)
    % Set the color of each axis to black
    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(2).Color = [0 0 0];

%     figure(3)
%     plot(t,y(:,2)+y(:,3)+y(:,4)+y(:,5))
%     legend('sum of all macrophages')

    %z0=[1-y0(1),0.1,0.1,0.1,0.01];
%     z0=[1-y0(1),0,0.1,0.1,0.01];
%     [t,y]=ode23s(@(t,y) rhs(t,y,params),tspan,z0);
%     figure(1)
%     plot(t,y(:,1),'r')
    


%     figure(4)
%     for f=0.1:0.1:1
%         f
%         [t,y]=ode23s(@(t,y) rhs(t,y,f),tspan,y0);
%         ss=sum(y(end,2:5));
%         plot(f,ss,'*')
%         hold on
%     end

    %finding steady states as roots of the rhs
    options = optimoptions('fsolve','SpecifyObjectiveGradient',true);
    %x0=[1,0.1,0.1,0.1,0.01];
    x0 = y0_low;
    if Case==1
        %x0=[0.01,0.006,0.06,0.006,0.0001];
        x0 = y0_low;
    end
    if Case==3
        %x0=[0.9,0.006,0.06,0.006,0.0001];
        x0 = y0_high;
    end
    if Case==2
        %x0=[0.9,0.006,0.06,0.006,0.0001];
        x0 = y0_high;
    end
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
    %dy(1)=params.pT*T*(1-T)*(params.h-params.r1*M1+params.r2*M2)-params.dT*T;
    %dy(1)=pT*T*(1-T)*(M2^n/(M2^n+M1^n)-0.5)-dT*T;
    %dy(1)=pT*M2^n/(M2^n+M1^n)-dT*T;
    %dy(1)=pT*T*(1-T)*(h+r2*M2)-dT*T*(h+r1*M1);
    %M0
    %dy(2)=p0*T*(1+g*M1)-d0*M0-alpha01*M0-...
    %dy(2)=p0*T*M0*(1-(M0+M1+M2+Mm)/f)*(1+g*M1)-d0*M0-alpha01*M0-...
    %dy(2)=p0*T*M0*(1-M0/f)*(1+g*M1)-d0*M0-alpha01*M0-...
    %dy(2)=p0*T*(1+g*M1)-d0*M0-alpha01*M0-...
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

    %parameters

    Jac=zeros(5,5);
    %Jac(1,1)=params.pT*(1-2*T)*(params.h-params.r1*M1+params.r2*M2)-params.dT;
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