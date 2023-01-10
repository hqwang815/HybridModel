function [interpolatedQHighResolution] = increase_resolution(interpolatedQ, samplePointsOld, samplePointsNew)


[mSources, nReflectionLevels] = size(interpolatedQ);
interpolatedQHighResolution{mSources, nReflectionLevels} = [];

for iSource = 1:1:mSources
    for jReflectionLevel = 1:1:nReflectionLevels

        if jReflectionLevel == nReflectionLevels
            thisQValue = interpolatedQ{iSource, jReflectionLevel};
            interpolatedQHighResolution{iSource, jReflectionLevel} = interp1(samplePointsOld, thisQValue, samplePointsNew, 'pchip', NaN);
        else
            thisQValue = interpolatedQ{iSource, jReflectionLevel}';
            interpolatedQHighResolution{iSource, jReflectionLevel} = interp1(samplePointsOld, thisQValue, samplePointsNew, 'pchip', NaN)';
        end
    end
end
end
