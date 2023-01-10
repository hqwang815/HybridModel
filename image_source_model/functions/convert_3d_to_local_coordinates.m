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
function [polygon2d] = convert_3d_to_local_coordinates(polygon3d)


for i = 1:1:size(polygon3d, 1)
    
    pol=polygon3d{i,1};

    p0 = pol(1, :);
    p1 = pol(2, :);
    p2 = pol(3, :);

    loc0 = p0;
    locx = p1 - loc0; 
  
    normal = cross(locx, p2-loc0);

    locy = cross(normal, locx);


    locx = locx / norm(locx);
    locy = locy / norm(locy);

    nVertices = size(polygon3d{i, 1}, 1);
    
    for j = 1:1:nVertices

        
        p1=dot(pol(j, :)-loc0, locx);
        p2=dot(pol(j, :)-loc0, locy);
        
        
        polygon2d{i, 1}(j, 1) = p1; 
        polygon2d{i, 1}(j, 2) = p2;

    end

end


end




% To go back from local coordinates (Lx, Ly) to 3D the transformation is
% p = loc0 + Lx*locx + Ly*locy
