% Description
% Calculates the surface area weighted scattering and absorption
% coefficients.

function [coefficientAverage]=calculate_average_coefficients(surfaceAreas,coefficients)




surfaceAreaNormalized=surfaceAreas/sum(surfaceAreas,1);
coefficientAverage=nansum(bsxfun(@times,surfaceAreaNormalized,coefficients));


end