function [t,x] = exponential_growth_pred_prey(cell_init, growth_rate ,death_rate)
%Simulate exponential growth with predator-prey dynamics
%
% Input:
%
%   cell_init       vector specifying the initial number of the predator
%                   and prey cells, respectively (e.g., [10, 10])
%
%   growth_rate     vector specifying the growth rate (1/h) for the
%                   predator and prey cells, respectively (e.g., [0.005, 0.2])
%
%   death_rate      vector specifying the death rate (1/h) for the
%                   predator and prey cells, respectively (e.g., [0.1, 0.01])
%
% Output:
%
%   t       vector of time points (h)
%
%   x       Nx2 matrix of the number of predator and prey cells at each
%           time point
%
% Usage:
%
%   [t,x] = exponential_growth_pred_prey(cell_init, growth_rate, death_rate);
%
%
%
% Author: Daniel Cook, 2018-10-01
% Updated: Jonathan Robinson, 2019-10-03
% Copyrighted under Creative Commons Share Alike


%% Section 1: Set parameter values 

% Set plotting and printing (true=show results, false=suppress results)
shouldPlot = true;

% Set start and end times (in hours), as well as the step size
tStart = 0;
tEnd = 100;

% Set number of cells to start with
x0(1) = cell_init(1);  % predator
x0(2) = cell_init(2);  % prey

% Define model parameters
k = growth_rate;
d = death_rate;


%% Section 2: Run model 

% Call ODE solver
[t,x] = ode45(@(t,x) growthFunc(t,x,k,d,x0), [tStart, tEnd], x0);


%% Section 3: Plot results 
% Plot results
if shouldPlot
    plot(t,x,'LineWidth',2);
    set(gca,'fontsize',14);
    xlabel('Time (hours)');
    ylabel('Level [A.U.]');
    legend('Predator','Prey');
end

end


%% Secion 4: Growth function 
function dxdt = growthFunc(t,x,k,d,x0)
% This function is the growth equation

% Set parameters
growth_rate = k;
death_rate = d;

% Set up differential equations
dxdt(1) = growth_rate(1)*x(1)*x(2) - death_rate(1)*x(1); % Predator population
dxdt(2) = growth_rate(2)*x(2) - death_rate(2)*x(2)*x(1); % Prey population

dxdt = dxdt'; % This needs to be a column vector
end



