function [t,x] = exponential_growth_pred_prey(cell_init, growth_rate ,death_rate)
% Author: Daniel Cook
% Date: 10/01/2018
% Copyrighted under Creative Commons Share Alike
% To run, use: [t,x] = exponential_growth_pred_prey([10,10],[.2,.005],[.01,.1]);

% Set plotting and printing (1=show results, 2=suppress results)
shouldPlot = 1; colorChoice = 'k';

%% Section 1: Set parameter values
% Set time (in days)
tStart = 0;
tEnd = 100; % Number of hours to run simulation

% Set number of cells to start
x0(1) = cell_init(1); % cell population 1
x0(2) = cell_init(2); % cell population 2

% Define model parameters
k = [growth_rate,death_rate];

%% Section 2: Run model
% Call ODE solver
timeStep = 0.1; % Days
[t,x] = ode45(@(t,x)odefun(t,x,k,x0), [tStart:tEnd], x0);

%% Section 3: Plot results
% Plot results
if shouldPlot == 1
    % Plot cell growth
    figure(1); hold on; plot(t,x(:,1),'-','color',colorChoice,'linewidth',2);
    figure(1);hold on; plot(t,x(:,2),'--','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (hours)'); ylabel('Level [A.U.]')
    legend('x1:Prey','x2:Predator')
end

end

function dxdt = odefun(t,x,k,x0)
% This function is the growth equation

% Set parameters
growth_rate = k(1:2);
death_rate = k(3:4);

% Set up differential equations
dxdt(1) = growth_rate(1)*x(1) - death_rate(1)*x(1)*x(2); % Cell population (1)
dxdt(2) = growth_rate(2)*x(2)*x(1) - death_rate(2)*x(2); % Cell population (2)


dxdt = dxdt'; % This needs to be a column vector
end