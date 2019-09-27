function result=solveMFA(Sin, measFlux, measValue, print)
% solveMFA
%   Solves metabolic flux analysis
%
%   Sin         S-matrix of internal metabolites
%   measFlux    indices of the measured fluxes
%   measValue   measured flux values
%   print       logical, whether output should be printed to command window
%               as strings (opt, default true)
%
%   result      matrix with flux indices and their value
%
%   Usage: result=solveMFA(Sin, measFlux, measValue, print)
%
%   Eduard Kerkhoven, 2019-09-26
%

if nargin<4
    print=true;
end

dims = size(measValue);
if dims(2)>dims(1)
    measValue=transpose(measValue);
end

fluxNo = 1:length(Sin);
calFlux = fluxNo;
calFlux(measFlux) = [];

Sinmeas = Sin(:,measFlux);
Sincal = Sin(:,calFlux);

vcal = (- Sincal^-1) * (Sinmeas * measValue); % Calculate the unknown fluxes

output(measFlux) = measValue;
output(calFlux) = vcal;

if print==true
    for i=1:numel(output)
    disp(['v' num2str(i) ': ' num2str(output(i))])
    end
end