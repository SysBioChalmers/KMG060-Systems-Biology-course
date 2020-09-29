%% Stoichiometry exercise v1.1

%% The complete S-matrix
S=[1,  -1,   X,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   1,   X,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   X,  -1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   X,   1,   0,   0,   0,   0,  -1,  -1,  -1,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   X,   0,   1,  -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   X,   0,   0,  -1,   1,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   1,  -1,  -1,   0,   0,   0,   0,   0;
   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X,   X;
   0,   0,   X,   1,   0,   0,   0,   0,   0,   5,   0,  -1,   0,  -2,   0,   0,   0,   0;
   0,   0,   X,   0,   2,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   0,   0,   0,  -2;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,  -1,   0,   0;
   0,   0,   X,   0,   1,   0,   0,   0,   0,   3,   1,   0,   0,   1,   0,   0,  -1,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,  -4.3;
   0,   0,   X,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -0.4;
  -1,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0;
   0,   0,   X,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1];
   
%% Question 4
% Define Sin and Sout
Sin  = S(1:14, :); % Rows are specified first, then columns. Take row 1
% to 14 (1:14) and all columns (:).
Sout = S(15:end, :); % Take rows 15 until the last (15:end) and all columns (:).

% Define which fluxes are measured and which will be calculated
measFluxes = [12, 16, 17, 18];
calFluxes = [1:11, 13:15];

% Generate the Sin,meas and Sin,cal matrices
Sinmeas = Sin(:,measFluxes);
Sincal = Sin(:,calFluxes);

% Provide the measured flux values, and determine the unknown fluxes
vmeas = [0; -5; 11; 0.05];
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

solveMFA(Sin, [12, 16, 17, 18], [0, -5, 11, 0.05]);

% Note, you should be able to do the matrix calculations by hand using
% elimination in the mass balances (with small models of less than ~10
% metabolites and reactions), so make sure you understand the code above,
% before using the function!

%% Question 5
% Make a new S matrix, where some of the metabolites have changed
% position. This can be done by defining the new order, for instance
% [2:end, 1], where the first entry is now placed last. Do such an
% arrangement of the metabolites in the S-matrix.

% Here is an example, but you can make any changes to the order of columns
% and rows.
Snew = S([1,3,4,2,5:end], [2:end, 1]); % Here the first four rows are
% shuffled, and the first row is placed last. You can test whether this
% S-matrix still givse the same solution by trying to solve question 4a
% again with this new S-matrix and compare the result.
