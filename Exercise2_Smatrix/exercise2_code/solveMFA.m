function result=solveMFA(S_in, measFlux, measValue, print)
% solveMFA
%   Solves metabolic flux analysis
%
%   S_in        S-matrix of internal metabolites
%   measFlux    indices of the measured fluxes
%   measValue   measured flux values
%   print       logical, whether output should be printed to command window
%               as strings (opt, default true)
%
%   result      matrix with flux indices and their value
%
%   Usage: result=solveMFA(S_in, measFlux, measValue, print)

if nargin<4
    print=true;
end

dims = size(measValue);
if dims(2)>dims(1)
    measValue=transpose(measValue);
end

fluxNo = 1:length(S_in);
calFlux = fluxNo;
calFlux(measFlux) = [];

df = size(S_in);
df = df(2) - df(1);

if length(measFlux)<df
    error('The number of measured fluxes is less than the degrees of freedom. Unable to determine a unique solution.')
end

Sinmeas = S_in(:,measFlux);
Sincal = S_in(:,calFlux);

if det(Sincal)==0
    error('The determinant of S_in,cal is zero. Unable to determine a unique solution.')
end

vcal = (- Sincal^-1) * (Sinmeas * measValue); % Calculate the unknown fluxes

output(measFlux) = measValue;
output(calFlux) = vcal;

if print==true
    for i=1:numel(output)
    disp(['v' num2str(i) ': ' num2str(output(i))])
    end
end