
%% Syntax
% # [output]=functionName(input)

%% Description
% Input arguments:
% # variable name  comment
% Output arguments:
% # variable name   comment

%%  Examples
% Example

%% Related functions %Related functions

%%
function [zRes] = ideal_membrane_resonator_impedance(parametersCox, frequencies, z0, c0, specificImpedance)


f = frequencies;

d = parametersCox(1);
m = parametersCox(2);
r = parametersCox(3);
c = c0;

if specificImpedance == true
    zRes = (r / z0) + (1i / z0) * (2 * pi * f * m - (z0 * cot((2 * pi * f)/c*d)));
elseif specificImpedance == false
    zRes = r + 1i * (2 * pi * f * m - (z0 * cot((2 * pi * f)/c*d)));
end


end