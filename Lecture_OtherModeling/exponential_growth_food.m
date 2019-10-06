function [t,x] = exponential_growth_food(cell_init, food_init, growth_rate)
%Simulate exponential growth with nutrient depletion
%
% Input:
%
%   cell_init       initial number of cells (e.g., 1)
%
%   food_init       initial amount of food (e.g., 10)
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
%   [t,x] = exponential_growth_food(cell_init, food_init, growth_rate);
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
tEnd = 48;

% Set number of cells and amount of food to start with
x0(1) = cell_init;
x0(2) = food_init;

% Define model parameters
k = growth_rate;


%% Section 2: Run simulation 

% Call ODE solver
[t,x] = ode45(@(t,x) growthFunc(t,x,k), [tStart, tEnd], x0);


%% Section 3: Plot results 

% Plot results
if shouldPlot
    plot(t,x,'LineWidth',2);
    set(gca,'fontsize',14);
    xlabel('Time (hours)');
    ylabel('Level [A.U.]');
    legend('cells','food');
end

end


%% Secion 4: Growth function 

function dxdt = growthFunc(t,x,k)
% This function is the growth equation

% Set parameters
growth_rate = k;

% Set up differential equations
dxdt(1) = growth_rate*x(1)*(1-x(1)/x(2)); % Cell population
% dxdt(2) = -0.1*x(1); % No feeding
dxdt(2) = 1 - 0.1*x(1); % Continuous feeding

dxdt = dxdt'; % This needs to be a column vector
end


