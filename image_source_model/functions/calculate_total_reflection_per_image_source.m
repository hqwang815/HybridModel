function [totalReflectionPerImageSource] = calculate_total_reflection_per_image_source(reflectionSequences, ...
    impedanceValues, cosineAngleOfIncidence, specificImpedance, angledep, z0)
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
totalReflectionPerImageSource{nSources, nLevels} = [];

%m?=source number n?=reflection level
for iSource = 1:1:nSources

    for jLevel = 1:1:nLevels
        [nImageSources, nReflections] = size(reflectionSequences{iSource, jLevel});

        for kImageSource = 1:1:nImageSources
            %             refarr{i,j}(k,:)=ones(size(refcof,2),1);
            totalReflectionPerImageSource{iSource, jLevel}(kImageSource, :) = ones(size(impedanceValues, 2), 1);

            for lReflection = 1:1:nReflections
                %                 refarr{i,j}(k,:)=refarr{i,j}(k,:).*refcof(refseq{i,j}(k,l),:);

                
                % Complex coefficient

                % Cosine of the Angle of Incidence at the surface
                thisCosineAngleOfIncidence = cosineAngleOfIncidence{iSource, jLevel}{kImageSource, lReflection};
                % Surface impedance at corresponfing surface
                
                
                ZS = impedanceValues(reflectionSequences{iSource, jLevel}(kImageSource, lReflection), :);

                % Equivalent plane wave reflection coefficient
                if specificImpedance == 0 && angledep == 1
                    
                    
                    R = (ZS * thisCosineAngleOfIncidence - z0) ./ (ZS * thisCosineAngleOfIncidence + z0);
                elseif specificImpedance == 1 && angledep == 1
                   
                    R = (ZS * thisCosineAngleOfIncidence - 1) ./ (ZS * thisCosineAngleOfIncidence + 1);
                                        
                    
                else
                    R = (ZS - 1) ./ (ZS + 1);
                end

                % Multiply by corresponding R
                totalReflectionPerImageSource{iSource, jLevel}(kImageSource, :) = ...
                    totalReflectionPerImageSource{iSource, jLevel}(kImageSource, :) .* R;


            end
        end
    end
end

end