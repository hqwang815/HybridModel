% HEADER COMMENTS

% INPUT:
% Vertices of convex polyhedron

% OUTPUT:
% Cell array with 3-point description of unique planes that make up convex
% polygon.
% Area per polygon making up polyhedron

% This function uses the Matlab boundary function to generate a 3d convex hull
% of the polyhedron it is described in planes that make up the convex hull (3 points per plane) it
% however outputs duplicate planes because it uses triangulation. In order to remove the duplicate planes
% this function calculates the reduced row echelon form of the vector description of all the planes
% then finds the duplicate planes and removes them. It also calculates the surface area for the unique planes.

function [polyhedronPolygons, wallSurfaceArea, volumeRoom] = create_unique_planes_description(PolyhedronVertices)


% Make a convex hull description of the vertex point cloud (which vertices belong to which plane)
[planeIndeces, volumeRoom] = boundary(PolyhedronVertices, 0);

% Amount of planes in convex hull description
nPlanes = size(planeIndeces, 1);

%% Calculating the areas of the triangles
triangulationAreas = zeros(nPlanes, 1);

for i = 1:1:nPlanes

    iPoint1 = planeIndeces(i, 1);
    iPoint2 = planeIndeces(i, 2);
    iPoint3 = planeIndeces(i, 3);

    A = PolyhedronVertices(iPoint1, :);
    B = PolyhedronVertices(iPoint2, :);
    C = PolyhedronVertices(iPoint3, :);

    v1 = B - A;
    v2 = C - A;

    triangulationAreas(i) = (1 / 2) * norm(cross(v1, v2));

end

%% Cell array containing the sets of 3-points describing each plane

triangulatedPlanes = cell(nPlanes, 1);

for i = 1:1:nPlanes
    for j = 1:1:3
        triangulatedPlanes{i, 1}(j, :) = PolyhedronVertices(planeIndeces(i, j), :);
    end
end


% Coefficientsrref=zeros(nPlanes,4);
[coefficientsPlane] = calculate_plane_coefficients(triangulatedPlanes);
% [coefficientsPlane]=calculateplanereducedrowform(triangulatedPlanes);

%%
% TURNED THIS PART INTO A SEPERATE FUNCTION ! ! WHEN FUNCTIONING PROPERLY
% THIS CAN BE REMOVED
% for i=1:1:nPlanes
%
%     % Take three points describing the plane
%     p1=triangulatedPlanes{i,1}(1,:);
%     p2=triangulatedPlanes{i,1}(2,:);
%     p3=triangulatedPlanes{i,1}(3,:);
%
%     % Take the plane normal
%     normal=cross(p1-p2,p1-p3);
%     syms x y z
%     P=[x,y,z];
%
% 	functionOfPlane=dot(normal, P-p1);
%     [planeCoefficients(1:1,1:3),planeCoefficients(1,4)]=equationsToMatrix(functionOfPlane,P);
%     planerref=rref(planeCoefficients);
%
% 	for j=1:1:4
%
%         thisPlaneCoefficientrref=sym2poly(planerref(1,j));
%         Coefficientsrref(i,j)=thisPlaneCoefficientrref;
%
% 	end
%
%     clear x y z
%
% end

%%

% 	Coefficientsrref(1:size(Coefficientsrref,1),1:3)=abs(Coefficientsrref(1:size(Coefficientsrref,1),1:3));
% 	Coefficientsrref=round(Coefficientsrref,4); Already done in the
% 	make a new function to determine if a plane is unique

[~, UniquePlanes, wallTriangulation] = unique(coefficientsPlane, 'rows');
nUniquePlanes = length(UniquePlanes);
wallSurfaceArea = zeros(nUniquePlanes, 1);

for i = 1:1:nUniquePlanes

    thisWallTriangles = wallTriangulation == i;
    wallSurfaceArea(i) = sum(triangulationAreas(thisWallTriangles));

end

polyhedronPolygons = cell(nUniquePlanes, 1);

for i = 1:1:nUniquePlanes

    thisPoint1 = planeIndeces(UniquePlanes(i), 1);
    thisPoint2 = planeIndeces(UniquePlanes(i), 2);
    thisPoint3 = planeIndeces(UniquePlanes(i), 3);

    polyhedronPolygons{i, 1}(1, :) = PolyhedronVertices(thisPoint1, :);
    polyhedronPolygons{i, 1}(2, :) = PolyhedronVertices(thisPoint2, :);
    polyhedronPolygons{i, 1}(3, :) = PolyhedronVertices(thisPoint3, :);

end


end