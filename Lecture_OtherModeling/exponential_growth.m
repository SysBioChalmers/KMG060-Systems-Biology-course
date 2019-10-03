function [t,x] = exponential_growth(cell_init, growth_rate)
% Author: Daniel Cook
% Date: 10/01/2018
% Copyrighted under Creative Commons Share Alike
% To run, use: [t,x] = exponential_growth(1,.2);

% Set plotting and printing (1=show results, 2=suppress results)
shouldPlot = 1; colorChoice = 'k';

%% Section 1: Set parameter values
% Set time (in hours)
tStart = 0;
tEnd = 48; % Number of hours to run simulation

% Set number of cells to start
x0 = cell_init;

% Define model parameters
k = growth_rate;

%% Section 2: Run model
% Call ODE solver
[t,x] = ode45(@(t,x)odefun(t,x,k,x0), [tStart:tEnd], x0);

%% Section 3: Plot results
% Plot results
if shouldPlot == 1
    % Plot cell growth
    figure(1); hold on; plot(t,x,'-','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (hours)'); ylabel('Number of cells')
%     legend('x1')
end

end

function dxdt = odefun(t,x,k,x0)
% This function is the growth equation

% Set parameters
growth_rate = k;

% Set up differential equations
dxdt = growth_rate*x;
end