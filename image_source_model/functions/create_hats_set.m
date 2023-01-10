% This function is used to create two spherical coordinate sets (Left ear and Right ear) for the
% directivity of a Head And Torso Simulator based on symmetry using the spherical coordinates of
% the left ear measurements (left ear at the origin with HATS facing in positive y-direction).


function [receiverDirectivityMeasurementCoordinatesHats,earPos] = create_hats_set(receiverDirectivityMeasurementCoordinates, flipAlong,centerPos)

if strcmpi(flipAlong, 'x')
    posLeft = [0, 1, 0];
    posRight = [0, -1, 0];
    vFlip = [1, -1, 1];
elseif strcmpi(flipAlong, 'y')
    posLeft = [1, 0, 0];
    posRight = [-1, 0, 0];
    vFlip = [-1, 1, 1];
elseif strcmpi(flipAlong, 'z')
    posLeft = [0, 0, 1];
    posRight = [0, 0, -1];
    vFlip = [1, 1, -1];
end

receiverDirectivityMeasurementCoordinatesHats{1,1} = receiverDirectivityMeasurementCoordinates;
%convert the left ear coordinates into the right ear coordinates (flip along y-axis)
receiverDirectivityMeasurementCoordinatesHats{2,1} = receiverDirectivityMeasurementCoordinates .* vFlip;


diameterEartoEar=0.215;
earPos{1,1} = posLeft*diameterEartoEar*(1/2);
earPos{2,1}  = posRight*diameterEartoEar*(1/2);

end