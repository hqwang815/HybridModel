function [reflectionSequences] = is_reflection_inside_polygon(reflectionSequences, polygonInPlane, pointsOfReflection, intersectionHasPolygon, poly, nPlanes, inPlane)

%inPLane is important
[mSources, nLevels] = size(polygonInPlane);

planePosition = (1:1:size(poly, 1));


for i = 1:1:mSources

    for j = 1:1:nLevels

        hasPolygon = intersectionHasPolygon{i, j};
        % Take the points of reflection that are in a wall with a polygon
        whichIntersections{i, j} = pointsOfReflection{i, j}(hasPolygon);
        % Take the wall numbers associated with these points of reflections
        whichWalls{i, j} = reflectionSequences{i, j}(hasPolygon);

        whichPolygon{i, j} = polygonInPlane{i, j}(hasPolygon);

    end

end

%
% 	counter=1;
% 	drawpolygons(poly,1)

for i = 1:1:mSources


    for j = 1:1:nLevels

        nPolygonWallIntersections = size(whichIntersections{i, j}, 1);

        indexen = find(intersectionHasPolygon{i, j});

        for k = 1:1:nPolygonWallIntersections


            %
            thisPolygonNo = whichPolygon{i, j}(k);
            thisWall = whichWalls{i, j}(k);

            isInPlane = inPlane == thisWall;

            %Hier moet ik even een aanpassing maken hier moet ik meerdere polygonen
            %kunnen checken misschien door ze per wand op te slaan in de cell matrix
            nPolygonsInPlane = nnz(isInPlane);
            thesePolygons = poly(isInPlane, 1);
            thesePolygonNo = planePosition(isInPlane);

            for l = 1:1:nPolygonsInPlane

                thisPolygon = thesePolygons{l, 1};
                thisPolygonNo = thesePolygonNo(l);

                thisIntersection = whichIntersections{i, j}{k};

                [aa] = is_point_in_polygon(thisPolygon, thisIntersection);

                if aa == true

                    break
                end


            end


            % 				scatter3(thisIntersection(1),thisIntersection(2),thisIntersection(3),'xk')
            % 				hold on
            % 				if aa==true
            % 				string=sprintf('%0.0f hit',counter);
            % 				counter=counter+1;
            % 				else
            % 				string=sprintf('%0.0f',counter);
            % 				counter=counter+1;
            % 				end
            %
            % 				text(thisIntersection(1),thisIntersection(2),thisIntersection(3),string)
            %


            if aa == true
                % Add the Number of the correct polygon where they are in the reflectionmatrix
                % 					thisPosition=planePosition(isInPlane);
                thisPosition = thisPolygonNo;


                reflectionSequences{i, j}(indexen(k)) = nPlanes + 1 + thisPosition;


            end


        end


    end


end


% 			[aa]=isPointInPolygon(p,S)


end