%% Stoichiometry exercise v1.0

%% The complete S-matrix
S=[1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 1, -1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 2, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 1, 0, 0, 0, 0, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, -1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0;
-1, 0, -1, 2, 0, 0, 0, 0, 0, 1, 0, 0, -2, 5, -1, 0, 0, -25;
0, 0, 0, 2, 0, 0, 0, 0, 0, 5, 0, -1, 0, -2, 0, 0, 0, 0;
0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, -2;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0;
0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 1, 0, 0, 1, 0, 0, -1, 0;
0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, -4.3;
0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.4;
-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

%% Question 6
% Define Sin and Sout
Sin  = S(1:14,:);
Sout = S(15:end, :);

% Define which fluxes are measured and which will be calculated
measFluxes = [1, 12, 16, 18];
calFluxes = [2:11, 13:15, 17];

% Generate the Sin,meas and Sin,cal matrices
Sinmeas = Sin(:,measFluxes);
Sincal = Sin(:,calFluxes);

% Provide the measured flux values, and determine the unknown fluxes
vmeas = [1; 0; -5; 0.05];
vcal = (- Sincal^-1) * (Sinmeas * vmeas);
% Print the fluxes of calculated fluxes
vcal

% Rearrange the fluxes
allFluxes(measFluxes) = vmeas;
allFluxes(calFluxes) = vcal;
% Give the n-th flux (here, n=13, giving us the rate of v13):
allFluxes(13)

% To make it a bit simpler, we've prepared a function called 'solveMFA',
% where we can give the following input: Sin matrix, the index number of
% those fluxes that are measured and the values of those measurements. The
% function performs the MFA and rearranges the solution to get an ordered
% list of fluxes.

solveMFA(Sin, [1, 12, 16, 18], [1, 0, -5, 0.05]);

% Note, you should be able to do the matrix calculations by hand (with
% small models of less than ~10 mtabolites and reactions), so make sure you
% understand the code above, before using the function!

%% Question 7
% Make a new Sin matrix, where some of the metabolites have changed
% position. This can be done by defining the new order, for instance
% [2:end, 1], where the first entry is now placed last. Do such an
% arrangement of the metabolites in the S-matrix
Sin_q7 = Sin(); % Change this to rearrange the metabolites, and then rerun
                % the analysis from question 6

%% Question 8
solveMFA(Sin, [1, 13, 17, 18], [0.275, 0, 1, 0.05]);

%% Question 10
solveMFA(Sin, [1, 13, 16, 18], [2.116, 0, 0, 0.05]);

%% Question 11
solveMFA(Sin, [1, 15, 12, 18], [10, 2.0597, 0, 0]);
