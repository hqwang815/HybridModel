function [interpolatedQ] = create_q_values(interpolatedDirectivity, DirectivityMeasurementPressure)
% This function generates Q factors from the interpolated pressure values

interpolatedQ{size(interpolatedDirectivity, 1), size(interpolatedDirectivity, 2)} = [];
meanPressure = mean(DirectivityMeasurementPressure, 1);
[mSources, nReflectionLevels] = size(interpolatedDirectivity);

for iSource = 1:1:mSources
    for jReflectionlevel = 1:1:nReflectionLevels
        interpolatedQ{iSource, jReflectionlevel} = interpolatedDirectivity{iSource, jReflectionlevel} ./ meanPressure;
    end
end

end