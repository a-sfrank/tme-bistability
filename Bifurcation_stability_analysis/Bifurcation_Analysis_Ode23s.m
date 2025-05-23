% This code plots bifrucation analysis for sensitive parameters

% Main contributor: Anna-Simone Frank

%==========================================================================
% DO NOT CHANGE
close all; clear all;

%load general information
myFolder = pwd;

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
%==========================================================================
% ADD/CHANGE IMPORTANT INFORMATION:
%---
% Choose between bi-or monostable cases: 1-7:

Case = 2

% Loading of prameters
params=parameters(Case);

%---
% Choose for bistable cases between high or low inital tumor state:

Tumor_Case= 2 % Can be 1 (high) or 2 (low)

%---
% Specify the bifurcation parameter e.g., r1,b,dT,f,K,alpha21,p0,d2,pT etc.

bif_param='p0'

%==========================================================================
% Default Initial Conditions
y0(1)=0.9;
y0(2)=0.01;
y0(3)=0.01;
y0(4)=0.01;
y0(5)=0.01;
%==========================================================================
% DO NOT CHANGE BELOW
%==========================================================================
% Set inital condition for run:
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
%% Bifurcation parameters
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


%==========================================================================
%Preallocate arrays for bifurcation points
 bifurcation_values = zeros(length(param_range), 5);

%Loop over parameter range and solve the ODE
for i = 1:length(param_range)
   param = param_range(i)

y0
    % Solve the ODE
    try

[t, y] = ode23s(@(t, y) odefun(t, y,param, params,Case,bif_param), [0 500], y0);


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
        bif_param_l = '$\mathbf{\beta}$'; 
    case 'r1'
        bif_param_l = '$\mathbf{r_1}$'; 
    case 'r2'
        bif_param_l = '$\mathbf{r_2}$'; 
    case 'p0'
        bif_param_l = '$\mathbf{p_0}$'; 
    case 'pT'
        bif_param_l = '$\mathbf{p_T}$';
    case 'dT'
        bif_param_l = '$\mathbf{d_T}$';
    case 'd2'
        bif_param_l = '$\mathbf{d_2}$';
    case 'alpha21'
        bif_param_l = '$\mathbf{\alpha_{21}}$';  
    otherwise
        bif_param_l = bif_param;
end


%-----------------------------------------------------------------------
figure;
ax = axes;

fontSize = 18;  % Common font size
fontName = 'Arial';  % Common font name

% Apply font settings to axes
set(ax, 'FontSize', fontSize, 'FontName', fontName);

% ---- Plot Tumor (Right Y-axis) ----
yyaxis right
hT = plot(param_range, bifurcation_values(:,1), 'k-', 'LineWidth', 2);
ylabel('tumor size steady state', 'FontSize', fontSize, 'FontName', fontName)
ylim([0, 1]);
hold off

% ---- Plot Immune Cells (Left Y-axis) ----
yyaxis left
h0 = plot(param_range, bifurcation_values(:,2), 'c-', 'LineWidth', 2);
hold on;
h1 = plot(param_range, bifurcation_values(:,3), 'r--', 'LineWidth', 2);
h2 = plot(param_range, bifurcation_values(:,4), 'b:', 'LineWidth', 2);
hM = plot(param_range, bifurcation_values(:,5), 'm-.', 'LineWidth', 2);
ylabel('cell population steady state', 'FontSize', fontSize, 'FontName', fontName)
ylim([0, 0.4]);
hold off

% ---- X-Axis and Labels ----
%xlabel(['Bifurcation for ', bif_param_l, xlabelMessage], ...
 %      'FontSize', fontSize, 'FontName', fontName, 'Interpreter', 'latex');
xlabel(bif_param_l, ...
       'FontSize', 22, 'FontName', fontName , 'Interpreter', 'latex');
legend([hT, h0, h1, h2, hM], {'Tumor','M0','M1','M2','Mm'});
grid on;
xlim([min(param_range), max(param_range)]);



% Set colors of both Y-axes
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];

% Optional font size
%fontsize(gcf, 18, 'points');


%==========================================================================
% Create a dynamic figure name based on these parameters
figure_name = ['Bifurcation_', bif_param, '_', figureMessage , '_', 'case','_ ',num2str(Case)]; 


% Define the file names with full path
base_file_name = fullfile(folder_name, figure_name);

% Save the figure in EPS, PDF, and JPEG formats
saveas(gcf, [base_file_name, '.eps'], 'epsc');  % EPS format
saveas(gcf, [base_file_name, '.pdf'], 'pdf');  % PDF format
saveas(gcf, [base_file_name, '.jpg'], 'jpg');  % JPEG format
saveas(gcf, [base_file_name, '.fig'], 'fig');  % FIG format

% Optionally, you can close the figure after saving
close(gcf);










