% Stoichiometry exercise v1.1

% The complete S-matrix
S=[1,  -1,   0,   0,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0
   0,   1,  -1,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   2,  -1,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   0,   1,   0,   0,   0,   0,  -1,  -1,  -1,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   0,   0,   1,  -2,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   0,   0,   0,  -1,   1,  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1,  -1,   0,   0,   0,   0,   0;
  -1,   0,  -1,   2,   0,   0,   0,   0,   0,   1,   0,   0,  -2,   5,  -1,   0,   0,  -25;
   0,   0,   0,   1,   0,   0,   0,   0,   0,   5,   0,  -1,   0,  -2,   0,   0,   0,   0;
   0,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,  -1,   0,   0,   0,   0,  -2;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -1,   0,  -1,   0,   0;
   0,   0,   0,   0,   1,   0,   0,   0,   0,   3,   1,   0,   0,   1,   0,   0,  -1,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,  -4.3;
   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,  -0.4;
  -1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0;
   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1];
   
% Question 4
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

% Question 5
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

% In the example rearrangement done above, Sout is still rows 15:end.
SinNew  = Snew(1:14, :);
SoutNew = Snew(15:end, :);
% But the columns all shifted one position (v12 has index 11, etc.).
solveMFA(SinNew, [11, 15, 16, 17], [0, -5, 11, 0.05]);

% Question 6
% Examples that error:
solveMFA(Sin, [1, 2, 3, 4], [10, 4, -3, 2]); % Zero determinant
solveMFA(Sin, [1, 2, 3], [10, 4, -3]); % Less then degrees of freedom

% Question 8
solveMFA(Sin, [1, 13, 17, 18], [0.675, 0, 3, 0.05]);

% Question 10
solveMFA(Sin, [1, 13, 16, 18], [3.75, 0, 0, 0.05]);
%solveMFA(Sin, [1, 13, 16, 18], [1.7737, 0, 0, 0.05]);

% Question 11
solveMFA(Sin, [1, 12, 15, 18], [1, 0, 6.1915, 0]);
solveMFA(Sin, [1, 12, 15, 18], [5, 0, 6.1915, 0]);
% See how much glucose is required for just ATP maintenance (no PHB, v13)
solveMFA(Sin, [13, 12, 15, 18], [0, 0, 6.1915, 0]);

% Question 13
solveMFA(Sin, [1, 12, 15, 18], [4, 0, 6.1915, 0.3]);
