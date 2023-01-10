
%% function name
% Function to generate the characteristic Impedance using Delany-Bazley-Miki model

%% Syntax
% # [Zs] = mikimodelimpedance(sigmaH, frange, c0,specificImpedance,z0)

%% Description
% Description
% sigmaH = matrix with air flow resisitivities and thicknesses for different surfaces (rows)
% fRange = frequency vector (Hz)
% specificImpedance = boolean to switch between specific (true) or wall impedance (false)
% z0 = characteristic impedance of air
% Output arguments
% Zs = 


%% Example(s)
% Related functions
function [Zs] = miki_model_impedance(sigmaH, fRange, c0,specificImpedance,z0)

omega = 2 * pi * fRange;

nSurfaces = size(sigmaH, 1);
nFrequencies = length(fRange);


Zs = ones(nSurfaces, nFrequencies);

for i = 1:1:nSurfaces

    thisSigma = sigmaH(i, 1);
    thisH = sigmaH(i, 2);


    zc = z0 * (1 + 5.5 .* (10^3 * (fRange / thisSigma)).^(-0.632) - 1i * 8.43 * (10^3 * (fRange / thisSigma)).^-(0.632));
    k = (omega / c0) .* (1 + 7.81 * (10^3 * (fRange / thisSigma)).^(-0.618) - 1i * 11.41 * (10^3 * (fRange / thisSigma)).^(-0.618));
    
    if specificImpedance==false
    Zs(i, :) = -1i * zc .* cot(k*thisH);
    else
    Zs(i, :) = -1i * zc .* cot(k*thisH)/z0;        
    end


end

end