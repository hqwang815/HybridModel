% Any rotation can be given as a composition of rotations about three axes (Euler's rotation theorem), and thus can be represented by a 3×3 matrix operating on a vector, 


function [rotatedSphericalCoordinatesSets]=rotate_points_around_axes(RotatePoints)
%function rotates a set of points described in spherical coordinates around
%x y or z axis

% The row size of v is the amount of sources or sets of unitsphere
% coordinares

sphericalCoordinateSets=RotatePoints.PointCoordinates;
rotationAroundAxesDeg=RotatePoints.RotationDegrees;
rotationAroundAxesOrder=RotatePoints.RotationOrder;


nPointSets=size(sphericalCoordinateSets,1);

for i=1:1:nPointSets
    
    % Degrees to radians
    xAxisAngleRad=(rotationAroundAxesDeg(i,1)*pi)/180;
    yAxisAngleRad=(rotationAroundAxesDeg(i,2)*pi)/180;
    zAxisAngleRad=(rotationAroundAxesDeg(i,3)*pi)/180;
    
    % Make rotation matrices
    xAxisRotation=[1,0,0;0,cos(xAxisAngleRad),-sin(xAxisAngleRad);0,sin(xAxisAngleRad),cos(xAxisAngleRad)];
    yAxisRotation=[cos(yAxisAngleRad),0,sin(yAxisAngleRad);0,1,0;-sin(yAxisAngleRad),0,cos(yAxisAngleRad)];
    zAxisRotation=[cos(zAxisAngleRad),-sin(zAxisAngleRad),0;sin(zAxisAngleRad),cos(zAxisAngleRad),0;0,0,1];
    
	% Allocate rotation matrix.
    rotationMatrix=cell(3,1);
	
    % Put rotation matrices in the desired order
    rotationMatrix{rotationAroundAxesOrder(i,1),1}=xAxisRotation;
    rotationMatrix{rotationAroundAxesOrder(i,2),1}=yAxisRotation;
    rotationMatrix{rotationAroundAxesOrder(i,3),1}=zAxisRotation;
    
    thisSphericalCoordinateSet=sphericalCoordinateSets{i,1};
    nSphericalCoordinates=size(thisSphericalCoordinateSet,1);
	
    % Perform rotations on all points
      
    

    for j=1:1: nSphericalCoordinates
        % Rotate around first axis
        rotatedSphericalCoordinatesSets{i,1}(j,:)=thisSphericalCoordinateSet(j,:)*rotationMatrix{1,1};
        % Rotate around second axis
        rotatedSphericalCoordinatesSets{i,1}(j,:)=rotatedSphericalCoordinatesSets{i,1}(j,:)*rotationMatrix{2,1};
        % Rotate around third axis
        rotatedSphericalCoordinatesSets{i,1}(j,:)=rotatedSphericalCoordinatesSets{i,1}(j,:)*rotationMatrix{3,1};
    end

end


end

