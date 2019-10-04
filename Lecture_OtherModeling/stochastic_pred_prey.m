function [t,x] = stochastic_pred_prey(cell_init, growth_rate ,death_rate)
%Simulate stochastic predator-prey dynamics
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
%   [t,x] = stochastic_pred_prey(cell_init, growth_rate, death_rate);
%
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

% Set number of cells to start with
x = cell_init;


%% Section 2: Run model
t_current = tStart;
t(1) = t_current;
i = 1;
while t_current <= tEnd
    % Determine which reaction(s) will occur
    % Set probabilities
    p1 = growth_rate(1) * x(i,2) * x(i,1); % Predator birth
    p2 = death_rate(1) * x(i,1); % Predator death
    p3 = growth_rate(2) * x(i,2); % Prey birth
    p4 = death_rate(2) * x(i,2) * x(i,1); % Prey death
    pTot = p1 + p2 + p3 + p4; % Total of rates
    
    % Generate random number (uniform random variable)
    reaction_rn = rand;
    
    % Select reaction
    if reaction_rn < p1/pTot  % Predator birth
        x(i+1,1) = x(i,1) + 1;
        x(i+1,2) = x(i,2);
    elseif reaction_rn < (p1+p2)/pTot  % Predator death
        x(i+1,1) = x(i,1) - 1;
        x(i+1,2) = x(i,2);
    elseif reaction_rn < (p1+p2+p3)/pTot  % Prey birth
        x(i+1,1) = x(i,1);
        x(i+1,2) = x(i,2) + 1;
    else                                % Prey death
        x(i+1,1) = x(i,1);
        x(i+1,2) = x(i,2) - 1;
    end
    
    % How long will a reaction take (exponential random variable)
    tau = exprnd(1/pTot);

    % step forward in time
    t_current = t_current+tau;
    t(i+1,1) = t_current;
    i = i+1;
end


%% Section 3: Plot results
% Plot results
if shouldPlot
    plot(t,x,'LineWidth',2);
    set(gca,'fontsize',14);
    xlabel('Time (hours)');
    ylabel('Level [A.U.]');
    legend('Predator','Prey');
end



