function [t,x] = stochastic_pred_prey(cell_init, growth_rate ,death_rate)
% Author: Daniel Cook
% Date: 10/01/2018
% Copyrighted under Creative Commons Share Alike
% To run, use: [t,x] = stochastic_pred_prey([10,10],[.2,.06],[.01,.1]);

% Set plotting and printing (1=show results, 2=suppress results)
shouldPlot = 0; colorChoice = 'k';

%% Section 1: Set parameter values
% Set time (in days)
tStart = 0;
tEnd = 500; % Number of hours to run simulation

% Set number of cells to start
x(1,1) = cell_init(1); % cell population 1
x(2,1) = cell_init(2); % cell population 2

% Define model parameters
k = [growth_rate,death_rate];

%% Section 2: Run model
t_step = tStart;
t(1) = t_step;
i = 1;
while t_step <= tEnd
    % Which reaction will occur
    % Set probabilities
    p1 = growth_rate(1)*x(1,i); % Prey birth
    p2 = death_rate(1)*x(1,i)*x(2,i); % Prey death
    p3 = growth_rate(2)*x(1,i)*x(2,i); % Predator birth
    p4 = death_rate(2)*x(2,i); % Predator death
    pTot = p1 + p2 + p3 + p4; % Total of rates
    
    % Generate random number (uniform random variable)
    reaction_rn = rand;
    
    % Select reaction
    if reaction_rn < p1/pTot
        x(1,i+1) = x(1,i) + 1;
        x(2,i+1) = x(2,i);
    elseif reaction_rn >= p1/pTot && reaction_rn < p2/pTot
        x(1,i+1) = x(1,i) - 1;
        x(2,i+1) = x(2,i);
    elseif reaction_rn >= p2/pTot && reaction_rn < p3/pTot
        x(1,i+1) = x(1,i);
        x(2,i+1) = x(2,i) + 1;
    else
        x(1,i+1) = x(1,i);
        x(2,i+1) = x(2,i) - 1;
    end
    
    % How long will a reaction take (exponential random variable)
    tau = exprnd(pTot);
%     tau = exprnd(1/pTot);

    t_step = t_step+tau;
    t(i+1) = t_step;
    i = i+1;
end

%% Section 3: Plot results
% Plot results
if shouldPlot == 1
    % Plot cell growth
    figure(1); hold on; plot(t,x(1,:),'-','color',colorChoice,'linewidth',2);
    figure(1);hold on; plot(t,x(2,:),'--','color',colorChoice,'linewidth',2);
    set(gca,'fontsize',18,'linewidth',2); box off
    xlabel('Time (hours)'); ylabel('Level [A.U.]')
    legend('x1:Prey','x2:Predator')
end

end