function [imageSourceReceiverDistances] = calculate_source_receiver_distance(effectiveImageSourcesPositions, sourcePositions, thisReceiverPosition)
% effectiveImageSourcesPositions, sourcePositions, thisReceiverPosition
% Function that calculates the distance from (Image) source to receiver.


[m, n] = size(effectiveImageSourcesPositions);
imageSourceReceiverDistances{m, n} = [];

for a = 1:1:m

    for b = 1:1:n + 1
        % For direct sound
        if b == n + 1
            thisVector = abs(thisReceiverPosition-sourcePositions(a, :));
            imageSourceReceiverDistances{a, b}(1, 1) = norm(thisVector);
        else
            % For reflected sound
            for c = 1:1:size(effectiveImageSourcesPositions{a, b}, 1)
                thisVector = abs(thisReceiverPosition-effectiveImageSourcesPositions{a, b}(c, :));
                imageSourceReceiverDistances{a, b}(c, 1) = norm(thisVector);
            end
        end
    end
end
end
