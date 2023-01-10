%% function name
% Short description
%% Syntax
% # [output argument(s)] = functionname(input arguments)
%% Description
% Description
% Input arguments
% Output arguments
%% Example(s)
% Related functions
function [polygonAreas] = calculate_polygon_areas(polygon2d)


nPolygons = size(polygon2d, 1);

for i = 1:1:nPolygons


    nVertices = size(polygon2d{i, 1}, 1);

    thisPolygon = polygon2d{i, 1};

    for j = 1:1:nVertices

        if j == nVertices

            x1 = thisPolygon(j, 1);
            y1 = thisPolygon(j, 2);

            x2 = thisPolygon(1, 1);
            y2 = thisPolygon(1, 2);

        else

            x1 = thisPolygon(j, 1);
            y1 = thisPolygon(j, 2);

            x2 = thisPolygon(j+1, 1);
            y2 = thisPolygon(j+1, 2);
            
        end

        v(j)=(x1*y2)-(x2*y1);
        
        
        

    end
    
    polygonAreas(i,1)=sum(v)/2;
    


end




end