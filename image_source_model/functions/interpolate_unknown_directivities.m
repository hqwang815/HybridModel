function [interpolatedDirectivity] = interpolate_unknown_directivities(orthodromicDistances, directivityMeasurementPressure, interpolationPowerParameter, nInterpolationClosestPoints)


% orthodromicDistancesSources, sourceDirectivityMeasurementPressure, interpolationPowerParameter, nInterpolationClosestPoints


[nTransducers, nReflectionLevels] = size(orthodromicDistances);

interpolatedDirectivity{nTransducers, nReflectionLevels} = [];

for iTransducer = 1:1:nTransducers

    for jReflectionLevel = 1:1:nReflectionLevels
        theseDistances = orthodromicDistances{iTransducer, jReflectionLevel};

        %size of the matrix containing the orthodromic distances
        nImageSourcesInLevel = size(theseDistances, 1);
        %size of the matrix containing the pressures for the 8bands and known
        %measurement points.
        nFrequencyBands = size(directivityMeasurementPressure, 2);
        %%%INVERSE DISTANCE INTERPOLATION%%%%%%%%%%

        for kImageSourceInLevel = 1:1:nImageSourcesInLevel
            %sort the distances for a given unknown point on length and the index they
            %belong to
            [theseDistancesSorted, theseDistancesSortedIndex] = sort(theseDistances(kImageSourceInLevel, :));

            %put exponent on the distances
            theseDistancesSortedExponentiated = theseDistancesSorted.^interpolationPowerParameter;
            %calculate the inverse distances for the amount of points to be used in
            %the interpolation 'n'.
            invertedDistances = 1 ./ theseDistancesSortedExponentiated(1:nInterpolationClosestPoints);

            for lFrequencyBand = 1:1:nFrequencyBands
                %take the first 4 closests known measurement values to mpoint i
                thesePressures = directivityMeasurementPressure(theseDistancesSortedIndex(1:nInterpolationClosestPoints), lFrequencyBand)';
                %multiply the inverse distances of the measurements with the corresponding
                %measurement values
                interpolatedDirectivity{iTransducer, jReflectionLevel}(kImageSourceInLevel, lFrequencyBand) = ...
                    sum(invertedDistances.*thesePressures) / sum(invertedDistances);
                
            end
        end
    end
end
end
