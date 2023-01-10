function [planeTwoVectors,planeSurfaceNormals]=generate_plane_surface_normals(walls)
% Make plane surface normals.

[m,n]=size(walls);

planeTwoVectors=cell(m,n);
planeSurfaceNormals=cell(m,n);

for i=1:1:m

% Make vector representations of the planes
planeTwoVectors{i,1}(1,:)=walls{i,1}(1,:)-walls{i,1}(2,:);
planeTwoVectors{i,1}(2,:)=walls{i,1}(1,:)-walls{i,1}(3,:);
                    
% Make surface normals of the vector representations
planeSurfaceNormals{i,1}=cross(planeTwoVectors{i,1}(1,:),planeTwoVectors{i,1}(2,:));
end

end