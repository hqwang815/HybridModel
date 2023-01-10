%% Syntax
% #  [room,verticesPolygons] = convert_obj_file_for_ism_model(objFileName,swapYzyIsUp)

%% Description
% Input arguments:
% # objFileName  name of the object file (name.obj)
% # swapYzyIsUp change position of y and z coordinates and flip y
% coordinates
% Output arguments:
% # room   vertices of convex polygon describing room
% # verticesPolygons vertices of the polygons on the boundary

%%  Examples
% Example

%% Related functions %Related functions

%%
function [room, verticesPolygons] = convert_obj_file_for_ism_model(objFileName, swapYzyIsUp)


% Code to extract component vertices from obj. file used with sketchup 2021
% professional

fid = fopen(objFileName);


iLine = 1;
iComponent = 1;

tline = 'dummy';

% As long as tline outputs a character array keep reading lines from the
% obj. file

while ischar(tline)

    [tline, ltout] = fgets(fid);
    textArray{iLine, 1} = tline;

    % Check if the current line contains the word Component. I so line
    % number is saved in lineNos
    k = strfind(tline, 'Component');
    if isempty(k) == false

        lineNos(iComponent) = iLine;
        iComponent = iComponent + 1;

    end


    iLine = iLine + 1;


end

%% Extract the parts of the obj files for each component

nLines = length(textArray);
nComponents = length(lineNos);


for iComponent = 1:nComponents


    from = lineNos(iComponent);
    if iComponent ~= nComponents
        till = lineNos(iComponent+1);
    else
        till = nLines;
    end
    textArrayPerComponent{iComponent, 1} = textArray(from:till, 1:1);


end

% Find the lines where the vertices are defined

iVertex = 1;
k = 1;

for iComponent = 1:nComponents


    thisTextArray = textArrayPerComponent{iComponent, 1};
    nLines = length(thisTextArray);


    for iLine = 1:nLines

        if k ~= iComponent
            iVertex = 1;
        end

        if iLine == 1
            k = iComponent;
        end


        thisText = thisTextArray{iLine, 1};
        p = strfind(thisText, 'v ');
        if p == 1

            position = sprintf('%0.0f,%0.0f', iComponent, iLine);
            verticesLocation{iComponent, 1}(iVertex, 1) = iLine;
            iVertex = iVertex + 1;


        end


    end


end

% Extract the verticescoordinates and make them into doubles.

for iComponent = 1:nComponents


    theseVerticesLocations = verticesLocation{iComponent, 1};
    thisComponent = textArrayPerComponent{iComponent, 1};
    nVertex = length(theseVerticesLocations);


    for iVertex = 1:nVertex
        thisLocation = theseVerticesLocations(iVertex);
        thisVertexCoordinatesChar = thisComponent{thisLocation, 1};

        thisVertexCoordinates = extractAfter(thisVertexCoordinatesChar, 'v ');
        thisVertexCoordinates = strsplit(thisVertexCoordinates, ' ');

        thisVertexCoordinates = str2double(thisVertexCoordinates);

        if swapYzyIsUp == true
            thisVertexCoordinatesNew(1) = thisVertexCoordinates(1);
            thisVertexCoordinatesNew(2) = thisVertexCoordinates(3) * -1;
            thisVertexCoordinatesNew(3) = thisVertexCoordinates(2);

            vertices{iComponent, 1}(iVertex, :) = thisVertexCoordinatesNew;
        else
            vertices{iComponent, 1}(iVertex, :) = thisVertexCoordinates;


        end


    end


end


for iComponent = 1:nComponents

    [~, vol] = boundary(vertices{iComponent, 1});
    if vol ~= 0
        k = iComponent;
        room = vertices{iComponent, 1};

    end

end

vertices(k) = [];

nPolygons = length(vertices);


for iPolygon = 1:nPolygons

    thisPolygonVertices = vertices{iPolygon, 1};
    [theseSortedPolygonVertices] = sort_polygon_vertices(thisPolygonVertices);
    verticesPolygons{iPolygon, 1} = theseSortedPolygonVertices;

end


% figureNumber = 15;
% plotpolyhedronboundary(room, figureNumber)
% drawpolygons(verticesPolygons, figureNumber)


end