% Check for impossible reflections. Reflections are impossible if the wall
% through which the Image source was created is not the closest wall from
% receiver position to (Image) source.

% pic   =   Plane through which the image source was created.
% rp    =   Receiver position.
% sp    =   Source position.
% PSN   =   Cell array with plane surface normals
% PL    =   Cell array with sets of 3 points defining planes.
% p0    =   Position on the plane [x,y,z].

function [imageSourcePlaneisClosestPlane,thisWallIntersection]=calculate_line_plane_intersections(planeSurfaceNormals,walls,position1,position2,thisImageSourcePlane)


nPlanes=size(walls,1);

% Preallocate array to hold d values.
d=zeros(1,nPlanes);
% Preallocate matrix to hold points of intersection with all planes.
pointsOfIntersection=zeros(nPlanes,3);

planeNumbers=(1:1:nPlanes);


for iPlane=1:1:nPlanes
% Calculate d's for intersection determination
p0=walls{iPlane,1}(1,:);

l=position1-position2;
l0=position1;
n=planeSurfaceNormals{iPlane,1};
d(1,iPlane)=dot((p0-l0),n)/(dot(l,n));
% Points of intersection
pointsOfIntersection(iPlane,:)=d(1,iPlane)*l+l0;

% scatter3(poi(z,1),poi(z,2),poi(z,3),'.k')
% 
% str=sprintf('%d',z);
% text(poi(z,1),poi(z,2),poi(z,3),str)
% scatter3(sp(1),sp(2),sp(3),'hr','filled')


end
                             
%If a reflection is parallel to a certain wall "d" will give an infinitive 
%or minus infinitive value. In this case label it as "Not a Number".

d(isinf(d)) = NaN;
%Check if there are wall interactions in the positive direction
% yy=sum(d(:) <= 0);
% If the sum of the number of elements in the array d that are negative 
% is not zero there is a plane in which a reflection is possible. 

reflectionIsPresent=sum(d(:) <=0) > 0;
d=round(d,4);


if reflectionIsPresent==true
	% find the closest plane intersection in the positive direction (smallest negative d)
	
% [thisClosestPlane]=find(d==max(d(d<0)));
% thisClosestPlane=thisClosestPlane(1);

thisClosestPlane=d==max(d(d<0));
thisClosestPlane=planeNumbers(thisClosestPlane);
thisClosestPlane=thisClosestPlane(1);

end

% Check if the if the wall through which the Image source was created is 
% not the closest wall from receiver position to (Image) source position.
if thisClosestPlane~=thisImageSourcePlane
imageSourcePlaneisClosestPlane=false;
end

if thisClosestPlane==thisImageSourcePlane
imageSourcePlaneisClosestPlane=true;
end

thisWallIntersection=pointsOfIntersection(thisClosestPlane,:);


end