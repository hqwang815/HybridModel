function [SOP] = check_point_side_of_plane(R, S)
% Test at which side of the plane the point lies

A = R(1, :);
B = R(2, :);
C = R(3, :);

Bd = B - A;
Cd = C - A;

Xd = S - A;
Smat = [Bd; Cd; Xd];
SOP = sign(det(Smat));

end
