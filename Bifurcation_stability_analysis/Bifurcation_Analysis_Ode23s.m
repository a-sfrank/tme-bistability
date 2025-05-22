close all; clear all;
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
Case = 1
params=parameters(Case);
Tumor_Case= 1 % Can be 1 or 2, it specifies high or low inital condition
bif_param='r1'

% Case = 3;
% 
% if Case == 3
%     fprintf('Running simulation for Case = %d: .\n', Case);
% elseif Case == 1
%         fprintf('Running simulation for Case = %d: .\n', Case);
% 
% end

%==========================================================================
% Initial Conditions
y0(1)=0.9;
y0(2)=0.01;
y0(3)=0.01;
y0(4)=0.01;
y0(5)=0.01;

switch Tumor_Case
    case {1}
        xlabelMessage = ' - High Tumor (IC)'
        figureMessage = 'High_IC'
        T0=y0(1); % high inital tumor
    case {2}
         xlabelMessage = ' - Low Tumor (IC)'
         figureMessage = 'Low_IC'
       T0=0.09; %low inital tumor
 otherwise
                disp('No such case')
end

M0_0=y0(2); M1_0=y0(3); M2_0=y0(4);Mm_0=y0(5);

y0=[T0, M0_0,M1_0, M2_0,Mm_0];

tspan = [0 500];
%==========================================================================
%% Bif. parameters
switch bif_param
    case 'b'
        param_range = 0.0:0.1:4.0; %beta range values
    case 'r1'
        param_range = 0.5:0.1:20.5; % r1 values range
    case 'r2'
        param_range = 0.5:0.1:20.5; % r2 values range
    case 'f'
        param_range = 0.1:0.2:2.4; % f values range
    case 'K'
        param_range = 0.4:0.1:1.1; % K range values
    case 'p0'
        param_range = 0.0:0.02:1.0; % p0 range values
    case 'pT'
         param_range = 0.0:0.1:0.5; % pT range values
    case 'dT'
        param_range = 0.0:0.01:0.1; % dT range values
    case 'd2'
        param_range = 0.0:0.01:0.2; % d2 range values
    case 'alpha21' 
        param_range = 0.0:0.02:0.5; % alpha21 range values
    otherwise
        disp('No such case')
end

% Alphas are boring for bifurcation analysis, so they are not used
%param_range = 0.1:0.1:1.1; % alpha01 range values
%param_range = 0.0:0.02:0.5; % alpha02 range values
%param_range = 0.0:0.02:0.5; % alpha12 range values
%param_range = 0.0:0.02:0.5; % alpha21 range values
%param_range = 0.001:0.001:0.5; % alpha1m range values
%param_range = 0.0:0.01:0.5; % alpham2 range values
%param_range = 0.0:0.02:0.5; % alpham1 range values
%param_range = 0.0:0.02:0.5; % alpha2m range values




%==========================================================================
%Preallocate arrays for bifurcation points
 bifurcation_values = zeros(length(param_range), 5);

%Loop over parameter range and solve the ODE
for i = 1:length(param_range)
   param =param_range(i)

    % Options for ode23s solver
    %options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
y0
    % Solve the ODE
    try

[t, y] = ode23s(@(t, y) odefun(t, y,param, params,Case,bif_param), [0 500], y0) %, options);


        % Extract the final point as bifurcation value
        bifurcation_values(i, :) = y(end, :);
    catch ME
        fprintf('Error at parameter value %f: %s\n', param, ME.message);
        % Handle error, perhaps use NaN to indicate failure
        bifurcation_values(i, :) = NaN;
    end
end


% Create a new folder for saving figures
folder_name = 'Figures';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

%==========================================================================
switch bif_param
    case 'b'
        bif_param_l = '$\beta$'; 
    case 'r1'
        bif_param_l = '$r_1$'; 
    case 'r2'
        bif_param_l = '$r_2$'; 
    case 'p0'
        bif_param_l = '$p_0$'; 
    case 'pT'
        bif_param_l = '$p_T$';
    case 'dT'
        bif_param_l = '$d_T$';
    case 'd2'
        bif_param_l = '$d_2$';
    case 'alpha21'
        bif_param_l = '$\alpha_{21}$'; % LaTeX format for alpha21   
    otherwise
        bif_param_l = bif_param;
end


% Plotting the bifurcation diagram
% figure;
% plot(param_range, bifurcation_values(:, 1), 'b', 'DisplayName', 'T','LineWidth', 3,'MarkerSize',4);
% hold on;
% plot(param_range, bifurcation_values(:, 2), 'r', 'DisplayName', 'M0','LineWidth', 3,'MarkerSize',4);
% plot(param_range, bifurcation_values(:, 3), 'g', 'DisplayName', 'M1','LineWidth', 3,'MarkerSize',4);
% plot(param_range, bifurcation_values(:, 4), 'c', 'DisplayName', 'M2','LineWidth', 3,'MarkerSize',4);
% plot(param_range, bifurcation_values(:, 5), 'm', 'DisplayName', 'Mm','LineWidth', 3,'MarkerSize',4);
% xlabel(['Bifurcation for ', bif_param_l,  xlabelMessage],'FontSize',16, 'Interpreter', 'latex');
% ylabel('Steady State Values','FontSize',16);
% title('');
% legend;
% ylim([0 1])
% grid on;

% Susannas plot
% figure=figure(1);
%     ax=axes;
%     yyaxis right
%     hT=plot(param_range,bifurcation_values(:,1),'k-','LineWidth', 3,'MarkerSize',4)
%     ylim([0,1])
%     ylabel('Steady State Values -- Tumor size')
%     xlabel(['Bifurcation for ', bif_param_l,  xlabelMessage],'FontSize',16, 'Interpreter', 'latex');
%    % xlabel('time')
%     %legend('T')
%     hold off
%     %figure(2)
%     yyaxis left
%     %h0=plot(t,y(:,2),'c','LineWidth',2);
%     h0=plot(param_range, bifurcation_values(:, 2), 'c', 'LineWidth', 3,'MarkerSize',4);
%     hold on
%     %h1=plot(t,y(:,3),'r--','LineWidth',2);
%     h1=plot(param_range, bifurcation_values(:, 3), 'r--', 'LineWidth', 3,'MarkerSize',4);
%     %h2=plot(t,y(:,4),'b:','LineWidth',2);
%     h2=plot(param_range, bifurcation_values(:, 4), 'b:', 'LineWidth', 3,'MarkerSize',4);
%     %hM=plot(t,y(:,5),'m-.','LineWidth',2);
%     hM=plot(param_range, bifurcation_values(:, 5), 'm-.', 'LineWidth', 3,'MarkerSize',4);
%     ylim([0,0.4])
%     ylabel('Steady State Values -- Cell populations')
%     hold off
%     xlim(tspan)
%     legend([hT,h0,h1,h2,hM],'Tumor','M0','M1','M2','Mm')
%     %fontsize(figure,16,'points')
%     %set(fig,'DefaultLineLineWidth',3)
%     % Set the color of each axis to black
%     ax.YAxis(1).Color = [0 0 0];
%     ax.YAxis(2).Color = [0 0 0];
%     grid on;
%-----------------------------------------------------------------------
figure;
ax = axes;

% ---- Plot Tumor (Right Y-axis) ----
yyaxis right
hT = plot(param_range, bifurcation_values(:,1), 'k-', 'LineWidth', 2);%, 'MarkerSize', 4);
ylabel('Tumor size steady state')
ylim([0, 1]);
hold off
% ---- Plot Immune Cells (Left Y-axis) ----
yyaxis left
h0 = plot(param_range, bifurcation_values(:,2), 'c-', 'LineWidth', 2);%, 'MarkerSize', 4);
hold on;
h1 = plot(param_range, bifurcation_values(:,3), 'r--', 'LineWidth', 2); %, 'MarkerSize', 4);
h2 = plot(param_range, bifurcation_values(:,4), 'b:', 'LineWidth', 2); %, 'MarkerSize', 4);
hM = plot(param_range, bifurcation_values(:,5), 'm-.', 'LineWidth', 2); %, 'MarkerSize', 4);
ylabel('Cell populations steady state')
ylim([0, 0.4]);
hold off
% ---- X-Axis and Labels ----
xlabel(['Bifurcation for ', bif_param_l, xlabelMessage], 'FontSize', 16, 'Interpreter', 'latex');
legend([hT, h0, h1, h2, hM], {'Tumor','M0','M1','M2','Mm'});
grid on;
xlim([min(param_range), max(param_range)]);
%xlim(tspan)

% Set colors of both Y-axes
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

% Optional font size
fontsize(gcf, 16, 'points');


% % Colorblind friendly
% % Create figure and axes
% figHandle = figure;
% ax = axes('Parent', figHandle);
% 
% % ---- Tumor on Right Y-axis ----
% yyaxis(ax, 'right');
% hT = plot(ax, param_range, bifurcation_values(:,1), '-', ...
%     'LineWidth', 2.5, 'Color', [0 0 0]); % Black
% ylabel(ax, 'Tumor Size (T)', 'Color', [0 0 0], 'FontSize', 14);
% ylim(ax, [0 1]);
% ax.YAxis(2).Color = [0 0 0]; % Color of the right Y-axis
% 
% % ---- Cell populations on Left Y-axis ----
% yyaxis(ax, 'left');
% hold(ax, 'on');
% h0 = plot(ax, param_range, bifurcation_values(:,2), '-', ...
%     'Color', [0.3 0.75 0.93], 'LineWidth', 2);  % Light Blue (M0)
% h1 = plot(ax, param_range, bifurcation_values(:,3), '--', ...
%     'Color', [0.85 0.33 0.1], 'LineWidth', 2);  % Reddish (M1)
% h2 = plot(ax, param_range, bifurcation_values(:,4), ':', ...
%     'Color', [0 0.45 0.74], 'LineWidth', 2);    % Blue (M2)
% hM = plot(ax, param_range, bifurcation_values(:,5), '-.', ...
%     'Color', [0.49 0.18 0.56], 'LineWidth', 2); % Purple (Mm)
% 
% ylabel(ax, 'Cell Populations (M0, M1, M2, Mm)', 'Color', [0 0 0], 'FontSize', 14);
% ylim(ax, [0 0.4]);
% ax.YAxis(1).Color = [0 0 0]; % Color of the left Y-axis
% %
% 
% % ---- X-axis setup ----
% xlabel(['', bif_param_l, ' ', xlabelMessage], ...
%     'Interpreter', 'latex', 'FontSize', 14);
% 
% %xlabel(['Bifurcation parameter: ', bif_param_l, ' ', xlabelMessage], ...
%   %  'Interpreter', 'latex', 'FontSize', 14);
% xlim(ax, [min(param_range), max(param_range)]);
% grid(ax, 'on');
% 
% % ---- Legend ----
% legend([hT, h0, h1, h2, hM], ...
%     {'Tumor (T)', 'M0', 'M1', 'M2', 'Mm'}, ...
%     'Location', 'best', 'FontSize', 12);
% 
% % Optional: Title
% %title(['Bifurcation Diagram for ', bif_param], 'Interpreter', 'none', 'FontSize', 16);
% 
% % Improve visual appeal
% set(gca, 'FontSize', 12);  % Axis tick font size
% box on;


%==========================================================================
% Create a dynamic figure name based on these parameters
figure_name = ['Bifurcation_', bif_param, '_', figureMessage , '_', 'case','_ ',num2str(Case)]; 


% Define the file names with full path
base_file_name = fullfile(folder_name, figure_name);

% Save the figure in EPS, PDF, and JPEG formats
saveas(gcf, [base_file_name, '.eps'], 'epsc');  % EPS format
saveas(gcf, [base_file_name, '.pdf'], 'pdf');  % PDF format
saveas(gcf, [base_file_name, '.jpg'], 'jpg');  % JPEG format

% Optionally, you can close the figure after saving
close(gcf);










