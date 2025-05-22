%% To plot Basin of Attraction for a given set of equations using ode45 solver
%% Bistable Case

clc; clear all
warning('off') 
%==========================================================================
%load general information
myFolder = pwd;

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
%==========================================================================
% Load important information
Case_basin = 1
%=========================================================================
% The roots of the given governing equations per parameter cases
% Makes only sense for bistable cases
%=========================================================================
%Case 3:
%r1 = [ 0.0121 ;   0.0063 ;   0.0627 ;   0.0067;    0.0006] ; %for bi-stability case
% r2 = [ 0.5492  ;  0.0037 ;   0.1211 ;   0.1175 ;   0.2057] ; %From
% population_model_v2
%---------------------------------------------
% Case 2: Check that SS are correct, population_Model_v2; Jacodian_dots_.m
% file calculates the jacobian matrix
 r1 = [0.5511    0.0343    0.1126    0.1095    0.1914]; 

 r2 = [0.0629    0.0200    0.0901    0.0365    0.0982];
%---------------------------------------------
%Case 4:
%r1 = [0.6286    0.0212    0.2075    0.2246    0.0036];
%r2 = [0.0575    0.0141    0.1387    0.0848    0.0015];
%---------------------------------------------
%Case 5:
% r1 = [0.0575    0.0141    0.1387    0.0849    0.0015];
%   r2 = [0.6293    0.0212    0.2075    0.2247    0.0036];
%---------------------------------------------

switch Case_basin
    case {1}
    y1 = linspace(0,2,50) ;
    y2 = linspace(0,2,50) ;
    y3 = 0.001;
    y4 = 0.001;
    y5 = 0.001;
    ylabelMessage = 'Macrophage (M0)';
    xlabelMassage = 'Tumor state';
    figureMessage = 'BoA_M0';
    case {2}
    y1 = linspace(0,2,50) ;
    y2 = 0.001;
    y3 = linspace(0,2,50);
    y4 = 0.001;
    y5 = 0.001;
    ylabelMessage = 'Macrophage (M1)';
    xlabelMassage = 'Tumor state';
    figureMessage = 'BoA_M1';
    case {3}
    y1 = linspace(0,2,50) ;
    y2 = 0.001;
    y3 = 0.001;
    y4 = linspace(0,2,50);
    y5 = 0.001;
    ylabelMessage = 'Macrophage (M2)';
    xlabelMassage = 'Tumor state';
    figureMessage = 'BoA_M2';
    case {4}
    y1 = linspace(0,2,50) ;
    y2 = 0.001;
    y3 = 0.001;
    y4 = 0.001;
    y5 = linspace(0,2,50);
    ylabelMessage = 'Macrophage (Mm)';
    xlabelMassage = 'Tumor state';
    figureMessage = 'BoA_Mm';
end
%==========================================================================

switch Case_basin
     %-----------------------------------------------------------------------
     %-----------------------------------------------------------------------
    case {1}
    disp('Running case 1: M0')
  Yr1 = [] ; Yr2 = [] ; 
% Time Span
tmax=80;
tspan=[0,tmax];

tic 
% varying inital condition for M1 (y3) and Tumor (y1), therewith
% investigating if  M1 subpopualtion levels contribute to tumor
% bistability.
for i = 1:length(y1)
     for j = 1:length(y2)
          y0 = [y1(i);y2(j);y3;y4;y5] ; %initial condition between 0 and 2 both for x1 and x2
          % Solve the system of Equations using ode45 solver
          
          %[T,Y]=ode45(@rhs, tspan,y0);
          [T,Y]=ode23s(@rhs,tspan,y0);
          FP=transpose(Y(end,:)) ; %modified 2 D: Y(end,1:2) or Y(end,1:2:3 

          % Locating the initial conditions according to error
          if norm(FP-r1)<1e-8 
             Yr1 = [y0 Yr1]  ;
          elseif norm(FP-r2)<1e-8
             Yr2 = [y0 Yr2] ;
          else
              x=[norm(FP-r1),norm(FP-r2)];
              I=find(x==min(x));
              if I==1
                  Yr1 = [y0 Yr1];
              elseif I==2
                  Yr2 = [y0 Yr2];
              end
          end
          
     end 
end
toc
%---------------------------------------------------------------------------
warning('on') % Remove the warning off constraint 
    % Initialize figure
    figure;
    %set(gcf,'color','w') 
    hold on
    plot(Yr1(1,:),Yr1(2,:),'.','color','r') ;
    plot(Yr2(1,:),Yr2(2,:),'.','color','b') ;
    plot(r1(1), r1(2),'k*','MarkerSize',20) ;
    plot(r2(1), r2(2),'kd','MarkerSize',10) ;

    legend('High Tumor (IC)','Low Tumor (IC)', 'fontsize', 12) ;
    xlabel('Tumor state', 'fontsize', 16) ;
    xlim([0 1])
    ylabel(ylabelMessage, 'fontsize', 16) ;
    legend;
    ylim([0 1])
    grid off;
    % Create a dynamic figure name based on these parameters
    figure_name = figureMessage; 
    folder_name = 'Figures'
    % Define the file names with full path
    base_file_name = fullfile(folder_name, figure_name);

    % Save the figure in EPS, PDF, and JPEG formats
    saveas(gcf, [base_file_name, '.eps'], 'epsc');  % EPS format
    saveas(gcf, [base_file_name, '.pdf'], 'pdf');  % PDF format
    saveas(gcf, [base_file_name, '.jpg'], 'jpg');  % JPEG format

    % Optionally, you can close the figure after saving
    close(gcf);
   %Add other 3 cases 5clean
   %-----------------------------------------------------------------------
    %-----------------------------------------------------------------------
     case {2}
    disp('Running case 2: M1')
  Yr1 = [] ; Yr2 = [] ; 
% Time Span
tmax=80;
tspan=[0,tmax];

tic 
% varying inital condition for M1 (y3) and Tumor (y1), therewith
% investigating if  M1 subpopualtion levels contribute to tumor
% bistability.
for i = 1:length(y1)
     for j = 1:length(y3)
          y0 = [y1(i);y2;y3(j);y4;y5] ; %initial condition between 0 and 2 both for x1 and x2
          % Solve the system of Equations using ode45 solver
          
          %[T,Y]=ode45(@rhs, tspan,y0);
          [T,Y]=ode23s(@rhs,tspan,y0);
          FP=transpose(Y(end,:))  %modified 2 D: Y(end,1:2) or Y(end,1:2:3 

          % Locating the initial conditions according to error
          if norm(FP-r1)<1e-8 
             Yr1 = [y0 Yr1]  ;
          elseif norm(FP-r2)<1e-8
             Yr2 = [y0 Yr2] ;
          else
              x=[norm(FP-r1),norm(FP-r2)];
              I=find(x==min(x));
              if I==1
                  Yr1 = [y0 Yr1];
              elseif I==2
                  Yr2 = [y0 Yr2];
              end
          end
          
     end 
end
toc
%---------------------------------------------------------------------------
warning('on') % Remove the warning off constraint 
    % Initialize figure
    figure;
    %set(gcf,'color','w') 
    hold on
    plot(Yr1(1,:),Yr1(3,:),'.','color','r') ;
    plot(Yr2(1,:),Yr2(3,:),'.','color','b') ;
    plot(r1(1), r1(3),'k*','MarkerSize',20)  ;
    plot(r2(1), r2(3),'kd','MarkerSize',10)  ;

    legend('High Tumor (IC)','Low Tumor (IC)', 'fontsize', 12) ;
    xlabel('Tumor state', 'fontsize', 16) ;
     xlim([0 1])
    ylabel(ylabelMessage, 'fontsize', 16) ;
    legend;
    ylim([0 1])
    grid off;

    % Create a dynamic figure name based on these parameters
    figure_name = figureMessage; 
    folder_name = 'Figures'
    % Define the file names with full path
    base_file_name = fullfile(folder_name, figure_name);

    % Save the figure in EPS, PDF, and JPEG formats
    saveas(gcf, [base_file_name, '.eps'], 'epsc');  % EPS format
    saveas(gcf, [base_file_name, '.pdf'], 'pdf');  % PDF format
    saveas(gcf, [base_file_name, '.jpg'], 'jpg');  % JPEG format

    % Optionally, you can close the figure after saving
    close(gcf);
   %-----------------------------------------------------------------------
   %-----------------------------------------------------------------------
    case{3}
        disp('Running case 3: M2')
       % Initialize the required matrices
Yr1 = [] ; Yr2 = [] ; 
% Time Span
tmax=80;
tspan=[0,tmax];

tic 
% varying inital condition for M1 (y3) and Tumor (y1), therewith
% investigating if  M1 subpopualtion levels contribute to tumor
% bistability.
for i = 1:length(y1)
     for j = 1:length(y4)
          y0 = [y1(i);y2;y3;y4(j);y5] ; %initial condition between 0 and 2 both for x1 and x2
          % Solve the system of Equations using ode45 solver
          
          %[T,Y]=ode45(@rhs, tspan,y0);
          [T,Y]=ode23s(@rhs,tspan,y0);
          FP=transpose(Y(end,:)); %modified 2 D: Y(end,1:2) or Y(end,1:2:3 

          % Locating the initial conditions according to error
          if norm(FP-r1)<1e-8 
             Yr1 = [y0 Yr1]  ;
          elseif norm(FP-r2)<1e-8
             Yr2 = [y0 Yr2] ;
          else
              x=[norm(FP-r1),norm(FP-r2)];
              I=find(x==min(x));
              if I==1
                  Yr1 = [y0 Yr1];
              elseif I==2
                  Yr2 = [y0 Yr2];
              end
          end
          
     end 
end
toc
warning('on') % Remove the warning off constraint 
    % Initialize figure
    figure;
    %set(gcf,'color','w') 
    hold on
    plot(Yr1(1,:),Yr1(4,:),'.','color','r') ;
    plot(Yr2(1,:),Yr2(4,:),'.','color','b') ;
    plot(r1(1), r1(4),'k*','MarkerSize',20)   ;
    plot(r2(1), r2(4),'kd','MarkerSize',10)    ;

    legend('High Tumor (IC)','Low Tumor (IC)', 'fontsize', 12) ;
    xlabel('Tumor state', 'fontsize', 16) ;
    xlim([0 1])
    ylabel(ylabelMessage, 'fontsize', 16) ;
    legend;
    ylim([0 1])
    grid off;
    % Create a dynamic figure name based on these parameters
    figure_name = figureMessage; 
    folder_name = 'Figures'
    % Define the file names with full path
    base_file_name = fullfile(folder_name, figure_name);

    % Save the figure in EPS, PDF, and JPEG formats
    saveas(gcf, [base_file_name, '.eps'], 'epsc');  % EPS format
    saveas(gcf, [base_file_name, '.pdf'], 'pdf');  % PDF format
    saveas(gcf, [base_file_name, '.jpg'], 'jpg');  % JPEG format

    % Optionally, you can close the figure after saving
    close(gcf);
    %-----------------------------------------------------------------------
   %-----------------------------------------------------------------------
    case {4}
        disp('Running case 4: Mm')
       % Initialize the required matrices
Yr1 = [] ; Yr2 = [] ; 
% Time Span
tmax=80;
tspan=[0,tmax];

tic 
% varying inital condition for M1 (y3) and Tumor (y1), therewith
% investigating if  M1 subpopualtion levels contribute to tumor
% bistability.
for i = 1:length(y1)
     for j = 1:length(y5)
          y0 = [y1(i);y2;y3;y4;y5(j)] ; %initial condition between 0 and 2 both for x1 and x2
          % Solve the system of Equations using ode45 solver
          
          %[T,Y]=ode45(@rhs, tspan,y0);
          [T,Y]=ode23s(@rhs,tspan,y0);
          FP=transpose(Y(end,:)) ; %modified 2 D: Y(end,1:2) or Y(end,1:2:3 

          % Locating the initial conditions according to error
          if norm(FP-r1)<1e-8 
             Yr1 = [y0 Yr1]  ;
          elseif norm(FP-r2)<1e-8
             Yr2 = [y0 Yr2] ;
          else
              x=[norm(FP-r1),norm(FP-r2)];
              I=find(x==min(x));
              if I==1
                  Yr1 = [y0 Yr1];
              elseif I==2
                  Yr2 = [y0 Yr2];
              end
          end
          
     end 
end
toc
warning('on') % Remove the warning off constraint 
    % Initialize figure
    figure;
    %set(gcf,'color','w') 
    hold on
    plot(Yr1(1,:),Yr1(5,:),'.','color','r') ;
    plot(Yr2(1,:),Yr2(5,:),'.','color','b') ;
    plot(r1(1), r1(5),'k*','MarkerSize',20)  ;
    plot(r2(1), r2(5),'kd','MarkerSize',10) ;

    legend('High Tumor (IC)','Low Tumor (IC)', 'fontsize', 12) ;
    xlabel('Tumor state', 'fontsize', 16) ;
    xlim([0 1])
    ylabel(ylabelMessage, 'fontsize', 16) ;
    legend;
    ylim([0 1])
    grid off;
    % Create a dynamic figure name based on these parameters
    figure_name = figureMessage; 
    folder_name = 'Figures'
    % Define the file names with full path
    base_file_name = fullfile(folder_name, figure_name);

    % Save the figure in EPS, PDF, and JPEG formats
    saveas(gcf, [base_file_name, '.eps'], 'epsc');  % EPS format
    saveas(gcf, [base_file_name, '.pdf'], 'pdf');  % PDF format
    saveas(gcf, [base_file_name, '.jpg'], 'jpg');  % JPEG format

    % Optionally, you can close the figure after saving
    close(gcf);
    end





%==========================================================================
function [dy,J]=rhs(t,y)

%parameters
%parameters_update;
 Case = 2 ;%input('Choose between 1 (low), 2 (bistable medium/low), 3 (high): ');
    params=parameters(Case) ;

    T=y(1);     %tumor cellf
    M0=y(2);    %M0 macrophagef
    M1=y(3);    %M1 macrophagef
    M2=y(4);    %M2 macophagef
    Mm=y(5);    %M12 mixed phenotype macrophagef

    
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

