function [impulseResponsePerImageSource, totalImpulseResponse, totalDirectFieldImpulseResponse] = generate_transfer_function(imageSourceReceiverDistances, ...
    frequencies, speedOfSound, totalReflectionPerImageSource, interpolatedQHighResolutionSource, ...
    interpolatedQHighResolutionReceiver, directivitySourceIsOn, directivityReceiverIsOn, phaseReverseIsOn,...
    totalScatteringPerImageSource, scatteringIsOn)
%This function performs the frequency domain construction of the impulse
%response
[nSources, nLevels] = size(imageSourceReceiverDistances);
% p = length(frequencies);

%
impulseResponsePerImageSource = cell(nSources, nLevels);


% generateFrequencyDomainImpulseResponse(imageSourceReceiverDistances, frequencies, ...
%     speedOfSound, totalReflectionPerImageSource, interpolatedQHighResolutionSource, ...
%     interpolatedQHighResolutionReceiver, directivitySourceIsOn, directivityReceiverIsOn, phaseReverseIsOn);


%source number
for iSource = 1:1:nSources
    %reflection level and direct path
    for jLevel = 1:1:nLevels
        %number of Image sources in reflection level
        nImageSources = size(imageSourceReceiverDistances{iSource, jLevel}, 1);


        for kImageSource = 1:1:nImageSources
            %the last collumn in the cell array with the source
            %receiver distances contains the distances for the direct
            %for the direct field no wall absorption is occuring
            if jLevel ~= nLevels

                if scatteringIsOn == true

                    thisTotalReflection = totalReflectionPerImageSource{iSource, jLevel}(kImageSource, :);
                    thisTotalScattering = totalScatteringPerImageSource{iSource, jLevel}(kImageSource, :);

                    thisTotalReflection = thisTotalScattering .* thisTotalReflection;


                elseif scatteringIsOn == false
                    thisTotalReflection = totalReflectionPerImageSource{iSource, jLevel}(kImageSource, :);
                end

            else
                thisTotalReflection = ones(1, length(frequencies));
            end
            %time delay is (image) source receiver distance divided by
            %the speed of sound
            t0 = imageSourceReceiverDistances{iSource, jLevel}(kImageSource) / speedOfSound;
            %apply source q if qsonoff is 1.
            if directivitySourceIsOn == 1
                qSource = interpolatedQHighResolutionSource{iSource, jLevel}(kImageSource, :);
            elseif directivitySourceIsOn == 0
                qSource = ones(1, length(frequencies));
            end
            %apply receiver q if gronoff is 1.
            if directivityReceiverIsOn == 1
                qReceiver = interpolatedQHighResolutionReceiver{iSource, jLevel}(kImageSource, :);
            elseif directivityReceiverIsOn == 0
                qReceiver = ones(1, length(frequencies));
            end
            %take (image) source receiver distance from matrix
            r = imageSourceReceiverDistances{iSource, jLevel}(kImageSource);
            %Phase part
            %Principle value of the phase shift
            phi = wrapToPi(frequencies*2*pi*t0);

            %Amplitude part
            %Calculate amplitude using directionality of source /
            %receiver, total reflection coefficient for image source
            %and distance attenuation
            A = (qSource .* qReceiver) .* thisTotalReflection / (4 * pi * r);

            %Option to enable phase inversion upon reflecting in
            %boundary (simple 180 degree inversion)
            if phaseReverseIsOn == 1

                impulseResponsePerImageSource{iSource, jLevel}(:, kImageSource) = A .* exp(1i.*phi) * exp(1i*jLevel*pi);

            else
                impulseResponsePerImageSource{iSource, jLevel}(:, kImageSource) = A .* exp(1i.*phi);


            end


        end

    end

end


% %%%%%%%%%%%%%%%%%%%%

%Summation of all (image) sources.
[nSources, nLevels] = size(impulseResponsePerImageSource);

for lSource = 1:1:nSources

    for mLevel = 1:1:nLevels
        impulseResponsePerLevel{lSource, 1}(:, mLevel) = sum(impulseResponsePerImageSource{lSource, mLevel}, 2);


    end
    impulseResponsePerSource(:, lSource) = sum(impulseResponsePerLevel{lSource, 1}, 2);

    ImpulseResponseDirectFieldPerSource(:, lSource) = impulseResponsePerImageSource{lSource, end};

    totalImpulseResponse = sum(impulseResponsePerSource, 2);
    totalDirectFieldImpulseResponse = sum(ImpulseResponseDirectFieldPerSource, 2);
    %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%
end

end