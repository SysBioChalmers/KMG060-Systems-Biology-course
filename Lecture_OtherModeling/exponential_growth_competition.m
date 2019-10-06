function [t,x] = exponential_growth_competition(cell_init, food_init, growth_rate)
%Simulate exponential growth with competition and nutrient depletion
%
% Input:
%
%   cell_init       vector specifying the initial number of each of the two
%                   cell types (e.g., [1, 1])
%
%   food_init       initial amount of food
%
%   growth_rate     vector specifying the growth rate (1/h) for each of the
%                   two cell types (e.g., [0.1, 0.2])
%
% Output:
%
%   t       vector of time points (h)
%
%   x       Nx2 matrix of number of each of the two cell types at each time
%           point
%
% Usage:
%
%   [t,x] = exponential_growth_competition(cell_init, food_init, growth_rate);
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
tEnd = 100;

% Set number of cells and amount of food to start with
x0(1) = cell_init(1);
x0(2) = cell_init(2);
x0(3) = food_init;

% Define model parameters
k = growth_rate;


%% Section 2: Run simulation

% Call ODE solver
[t,x] = ode45(@(t,x) growthFunc(t,x,k,x0), [tStart, tEnd], x0);


%% Section 3: Plot results
% Plot results
if shouldPlot
    plot(t,x,'LineWidth',2);
    set(gca,'fontsize',14);
    xlabel('Time (hours)');
    ylabel('Level [A.U.]');
    legend('cell 1','cell 2','food');
end

end


%% Secion 4: Growth function
function dxdt = growthFunc(t,x,k,x0)
% This function is the growth equation

% Set parameters
growth_rate = k;

% cells do not inhibit each other
dxdt(1) = growth_rate(1)*x(1)*(1-(x(1)+x(2))/x(3)); % Cell population (1)
dxdt(2) = growth_rate(2)*x(2)*(1-(x(1)+x(2))/x(3)); % Cell population (2)
dxdt(3) = 1 - 0.1*(x(1)+x(2)); % Continuous glucose feeding

% cell 2 inhibits cell 1 growth
% dxdt(1) = (growth_rate(1)/x(2))*x(1)*(1-(x(1)+x(2))/x(3)); % Cell population (1)
% dxdt(2) = growth_rate(2)*x(2)*(1-(x(1)+x(2))/x(3)); % Cell population (2)
% dxdt(3) = 1 - 0.1*(x(1)+x(2)); % Continuous glucose feeding

% Mutual inhibition
% dxdt(1) = (growth_rate(1)/x(2))*x(1)*(1-(x(1)+x(2))/x(3)); % Cell population (1)
% dxdt(2) = (growth_rate(2)/x(1))*x(2)*(1-(x(1)+x(2))/x(3)); % Cell population (2)
% dxdt(3) = 1 - 0.1*(x(1)+x(2)); % Continuous glucose feeding

dxdt = dxdt'; % This needs to be a column vector
end

