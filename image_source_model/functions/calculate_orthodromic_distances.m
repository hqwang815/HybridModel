function [orthodromicDistances] = calculate_orthodromic_distances(anglesOfincidenceSource, positionsDirectivityMeasurement)
% This function is used to determine the orthodromic distances between
% points on a unit sphere. as an input it uses two sets of points an as an
% output it gives the distances of each point of set A to all points of set
% B

% angleOfIncidenceSource, rotatedSourceDirectivityMeasurementCoordinates


[~, nReflectionLevels] = size(anglesOfincidenceSource);

mPointSets = size(positionsDirectivityMeasurement,1);

% preallocate orthodromic distance cell matrix
orthodromicDistances{mPointSets, nReflectionLevels} = [];

for iPointSet = 1:1:mPointSets

        knownPositions = positionsDirectivityMeasurement{iPointSet, 1};
        nKnownPositions=size(knownPositions, 1);
    
    
    for jReflectionlevel = 1:1:nReflectionLevels
        % Take point to be interpolated
        theseAnglesOfIncidence = anglesOfincidenceSource{iPointSet, jReflectionlevel};
        
        
        nPoints=size(theseAnglesOfIncidence,1);
        

        for kPoint = 1:1:nPoints
            
            thisAngleOfIncidence=theseAnglesOfIncidence(kPoint, :);
           
            for lKnownPosition = 1:1:nKnownPositions
                
                
                thisKnownPosition=knownPositions(lKnownPosition, :);
                
                crossProductOfVector = norm(cross(thisAngleOfIncidence, thisKnownPosition));
                dotProductOfVector = dot(thisAngleOfIncidence, thisKnownPosition);
                
                if theseAnglesOfIncidence(kPoint, :) == knownPositions(lKnownPosition, :)
                    orthodromicDistances{iPointSet, jReflectionlevel}(kPoint, lKnownPosition) = 0.001;
                else
                    orthodromicDistances{iPointSet, jReflectionlevel}(kPoint, lKnownPosition) = atan2(crossProductOfVector, dotProductOfVector);
                end
                
            end
            
%              plotorthodromicdistances(knownPositions, thisAngleOfIncidence,orthodromicDistances{i, j}(kPoint,:));
            
            
            
        end
    end
end


end