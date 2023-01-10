% GENERAL DESCRIPTION
% Checks in which planes which polygons are. the Inputs to the function are
% the sets of 3 points from each polygon and from each wall, this function
% will use the calculateplanereducedrowform function to determine the
% reduced row echelon form description of the planes in which the polygons
% and walls lie, if a polygon and a wall have the same description then
% the polygon lies in the wall.




function [inPlane]=find_polygon_in_plane(polygon,walls)

% [wallCoefficients]=calculateplanereducedrowform(walls);
% [polygonCoefficients]=calculateplanereducedrowform(poly);


% Extract 3 points from polygon

nPolygons = size(polygon, 1);
polyPlane = cell(nPolygons, 1);

for i = 1:1:nPolygons
    thisPolygon = polygon{i, 1};
    % Taking three points out of each polygon to describe the plane in which th
    % e polygon lies.
    poly{i, 1} = thisPolygon(1:3, 1:3);
end

















[wallCoefficients]=calculate_plane_coefficients(walls);
[polygonCoefficients]=calculate_plane_coefficients(poly);


% coefficie
% coefficientsrrefpoly

% coefficientsrrefplane(1:size(coefficientsrrefplane,1),1:3)=abs(coefficientsrrefplane(1:size(coefficientsrrefplane,1),1:3));
% coefficientsrrefplane=round(coefficientsrrefplane,4);
% 
% coefficientsrrefpoly(1:size(coefficientsrrefpoly,1),1:3)=abs(coefficientsrrefpoly(1:size(coefficientsrrefpoly,1),1:3));
% coefficientsrrefpoly=round(coefficientsrrefpoly,4);

% Check which polygons lie in which walls.






[~,inPlane]=ismember(polygonCoefficients,wallCoefficients,'rows');
% [~,inPlane]=ismember(coefficientsrrefpoly,coefficientsrrefplane,'rows');



end