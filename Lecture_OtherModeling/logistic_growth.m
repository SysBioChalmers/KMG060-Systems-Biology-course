function [t,x] = logistic_growth(cell_init, growth_rate)
%Simulate logistic cell growth
%
% Input:
%
%   cell_init       initial number of cells (e.g., 1)
%
%   growth_rate     growth rate (1/h) (e.g., 0.2)
%
% Output:
%
%   t       vector of time points (h)
%
%   x       vector of number of cells at each time point
%
% Usage:
%
%   [t,x] = logistic_growth(cell_init, growth_rate);
%
%
% Author: Daniel Cook, 2018-10-01
% Updated: Jonathan Robinson, 2019-10-03
% Copyrighted under Creative Commons Share Alike 


%% Section 1: Set parameter values

%Set plotting and printing (true=show results, false=suppress results)
shouldPlot = true;

% Set start and end times (in hours), as well as the step size
tStart = 0;
tEnd = 48;
tStep = 0.1;

% Set number of cells to start with
x0 = cell_init;

% Define model parameters
k = growth_rate;


%% Section 2: Run simulation

% Call ODE solver
[t,x] = ode45(@(t,x) growthFunc(t,x,k,x0), tStart:tStep:tEnd, x0);


%% Section 3: Plot results
% Plot results
if shouldPlot
    plot(t,x,'LineWidth',2);
    set(gca,'fontsize',14);
    xlabel('Time (hours)');
    ylabel('Number of cells');
end

end


%% Secion 4: Growth function

function dxdt = growthFunc(t,x,k,x0)
% This function is the growth equation

% Set parameters
growth_rate = k;
carrying_capacity = 10;

% Set up differential equations
dxdt = growth_rate * x * (1 - x/carrying_capacity); % cell population
end


