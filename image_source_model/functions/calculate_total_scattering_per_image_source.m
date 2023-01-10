function [totalScatteringPerImageSource] = calculate_total_scattering_per_image_source(reflectionSequences, ...
    scatteringCoefficients)
% determine the size of the reflection sequence array. the reflection
% sequence gives the sequence of wall numbers through which a certain image
% source reflected
% totalReflectionPerImageSource
% reflectionSequences, interPolatedReflectionCoefficients
% reflectionSequences, interpolatedReflectionCoefficients, anglesOfIncidenceWall,specificOrCharacteristicImpedance,angleDependance

% (reflectionSequences,...
%     interpolatedReflectionCoefficients, anglesOfIncidenceWall,specificImpedance,angleDependance,z0)

[nSources, nLevels] = size(reflectionSequences);
% preallocate memory for the absorption coefficients array
% refarr{m,n}=[];
totalScatteringPerImageSource{nSources, nLevels} = [];

%m?=source number n?=reflection level
for iSource = 1:1:nSources

    for jLevel = 1:1:nLevels
        [nImageSources, nReflections] = size(reflectionSequences{iSource, jLevel});

        for kImageSource = 1:1:nImageSources
            %             refarr{i,j}(k,:)=ones(size(refcof,2),1);
            totalScatteringPerImageSource{iSource, jLevel}(kImageSource, :) = ones(size(scatteringCoefficients, 2), 1);

            for lReflection = 1:1:nReflections
                %                 refarr{i,j}(k,:)=refarr{i,j}(k,:).*refcof(refseq{i,j}(k,l),:);

                
      
                
                thisS = scatteringCoefficients(reflectionSequences{iSource, jLevel}(kImageSource, lReflection), :);
                R = (1-thisS);
      

                % Multiply by corresponding R
                totalScatteringPerImageSource{iSource, jLevel}(kImageSource, :) = ...
                    totalScatteringPerImageSource{iSource, jLevel}(kImageSource, :) .* R;


            end
        end
    end
end

end