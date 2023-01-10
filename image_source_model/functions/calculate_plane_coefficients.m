% GENERAL DESCRIPTION
% input: plane descriptions in 3 points
% output: rref of plane description ax+by+cz+d

function [planeCoefficients] = calculate_plane_coefficients(planes)


nPlanes = size(planes, 1);

planeCoefficients = zeros(nPlanes, 4);

for i = 1:1:nPlanes

    % Take three points describing the plane
    p1 = planes{i, 1}(1, :);
    p2 = planes{i, 1}(2, :);
    p3 = planes{i, 1}(3, :);

    %     p1=[1,-2,0];
    %     p2=[3,1,4];
    %     p3=[0,-1,2];


    v1 = p2 - p1;
    v2 = p3 - p1;


    % Take the plane normal
    normal = cross(v1, v2);

    %     signCoefficients=sign(normal);
    %     idx=find(signCoefficients~=0,1);
    %     thisSignCoefficient=signCoefficients(idx);
    %         if thisSignCoefficient==-1
    %             normal=normal*-1;
    %         end


    % Unity plane normal
    % normal=normal/norm(normal);

    d = normal(1) * p1(1) + normal(2) * p1(2) + normal(3) * p1(3);
    planeCoefficient = [normal(1), normal(2), normal(3), d];


    planeCoefficient = planeCoefficient / max(abs(planeCoefficient));
    planeCoefficient = round(planeCoefficient, 4);
    signCoefficients = sign(planeCoefficient);
    idx = find(signCoefficients ~= 0, 1);
    thisSignCoefficient = signCoefficients(idx);


    if thisSignCoefficient == -1
        planeCoefficient = planeCoefficient * -1;
    end

    planeCoefficients(i, :) = planeCoefficient;


end

% planeCoefficients=round(planeCoefficients,4);

end