% nSources,sourcesNames,Ism,walls,thisReceiverPosition,planesDirections


function [Ism2, nodeDepth, pointsOfReflection, cosinesAnglesOfIncidence, reflectionSequences] = check_rule_2(nSources, sourcesNames, Ism, walls, thisReceiverPosition, planesDirections)

% Work In Progress!


% Function to make matrix with plane vectors (Pv2) and plane surface
% normals (PSN). Based on the point representations of planes in PL cell array.


[~, planeSurfaceNormals] = generate_plane_surface_normals(walls);


for iSource = 1:1:nSources
    % Make a tree indicating the depth of each node.
    nodeDepth = Ism.(sourcesNames{iSource}).isn.depthtree;
    % display(nodeDepth.tostring)

    % Depth of the tree (take last entry this is deepest level)
    nImageSourcesTotal = length(Ism.(sourcesNames{iSource}).isn.Node);
    nReflectionLevels = nodeDepth.Node{nImageSourcesTotal};

    imageSourcePositions = Ism.(sourcesNames{iSource}).isp;
    imageSourcePlanes = Ism.(sourcesNames{iSource}).isn;
    imageSourceNaming = Ism.(sourcesNames{iSource}).isn;
    namingTree = Ism.(sourcesNames{iSource}).isnms;


    for iReflectionlevel = 1:1:nReflectionLevels
        % Find all Image sources in level nr:(iReflectionLevel).
        imageSources = find(nodeDepth == iReflectionlevel);

        % Determine the amount of image sources in level.
        nImageSourcesInLevel = length(imageSources);

        % (Re)set counter for image source number in level.
        i = 1;
        for iImageSource = 1:1:nImageSourcesInLevel
            % Take image source number
            thisImageSourceNumber = imageSources(iImageSource);
            ImageSourceLineage = thisImageSourceNumber;
            x = 1;
            % Determine the lineage of the given image source

            while ImageSourceLineage(x) ~= 1
                ImageSourceLineage(x+1) = getparent(Ism.(sourcesNames{iSource}).isn, ImageSourceLineage(x));
                x = x + 1;

            end
            % Set counter for reflection array
            j = 1;

            % Loop to check if the path from receiver to source is possible

            nImageSourceLineage = length(ImageSourceLineage);

            for iImageSourceLineage = 1:1:nImageSourceLineage - 1
                thisLineageNode = ImageSourceLineage(iImageSourceLineage);
                thisImageSourcePlane = Ism.(sourcesNames{iSource}).isn.Node{thisLineageNode, 1};


                if iImageSourceLineage == 1
                    % 1st Reverse trace is image source to image source parent.
                    position1 = thisReceiverPosition;
                    position2 = Ism.(sourcesNames{iSource}).isp.Node{thisLineageNode, 1};
                else
                    % Later Reverse traces are from intersection to image source parent.
                    position1 = thisWallIntersection;
                    position2 = Ism.(sourcesNames{iSource}).isp.Node{thisLineageNode, 1};
                end

                %             imageSourcePlaneisClosestPlane,thisPointOfIntersection
                % Check if the plane through which was reflected is the first plane
                % encountered upon reverse trace
                [imageSourcePlaneIsClosestPlane, thisWallIntersection] = calculate_line_plane_intersections(planeSurfaceNormals, walls, position1, position2, thisImageSourcePlane);

                if imageSourcePlaneIsClosestPlane == false
                    imageSourcePlanes = imageSourcePlanes.set(thisImageSourceNumber, 0);
                    imageSourcePositions = imageSourcePositions.set(thisImageSourceNumber, [0, 0, 0]);
                    clear intersectionsInSequence thisCosineAngleOfIncidence thisImageSourceReflectionSequence
                    break

                else
                    l = position1 - position2;
                    nl = l / norm(l);
                    thisPlaneSurfaceNormal = planeSurfaceNormals{thisImageSourcePlane, 1} / norm(planeSurfaceNormals{thisImageSourcePlane, 1});
                    thisCosineAngleOfIncidence{1, j} = dot(nl, thisPlaneSurfaceNormal) * planesDirections(thisImageSourcePlane) * -1;
                    thisImageSourceReflectionSequence(1, j) = thisImageSourcePlane;
                    intersectionsInSequence{1, j} = thisWallIntersection;
                    j = j + 1;
                end

            end


            % Naming of the Image Sources (Maybe take this out of the code
            % somehow since it will slow down the code).

            imageSourcePlaneNo = flipud(cell2mat(imageSourceNaming.Node(ImageSourceLineage)));
            imageSourcePlaneNo = imageSourcePlaneNo(2:end);

            string1 = sprintf('.%d', imageSourcePlaneNo);
            string2 = strcat('$$ IS_{', string1, '} $$');
            namingTree = namingTree.set(thisImageSourceNumber, string2);
            clear string2


            if exist('intersectionsInSequence', 'var') ~= 0

                pointsOfReflection{iSource, iReflectionlevel}(i, :) = intersectionsInSequence;
                cosinesAnglesOfIncidence{iSource, iReflectionlevel}(i, :) = thisCosineAngleOfIncidence;
                reflectionSequences{iSource, iReflectionlevel}(i, :) = thisImageSourceReflectionSequence;
                i = i + 1;
            end

        end

    end

    Ism2.(sourcesNames{iSource}).isn = imageSourcePlanes;
    Ism2.(sourcesNames{iSource}).isp = imageSourcePositions;
    Ism2.(sourcesNames{iSource}).isnms = namingTree;

end


end
