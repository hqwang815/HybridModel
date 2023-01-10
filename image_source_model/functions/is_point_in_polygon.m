% General Description
% Function to check whether a point on a plane lies within a polygon on
% that plane.


function [aa] = is_point_in_polygon(p, S)

c = sum(p, 1) / size(p, 1);

for i = 1:1:size(p, 1)

    v3 = p(i, :) - c;

    if i == size(p, 1)
        v4 = p(1, :) - c;
    else
        v4 = p(i+1, :) - c;
    end

    v3 = v3 / norm(v3);
    v4 = v4 / norm(v4);

    n = cross(v3, v4);
    n = n / norm(n);
    A = p(i, :);

    if i == size(p, 1)
        B = p(1, :);
    else
        B = p(i+1, :);
    end

    C = p(i, :) + n;

    PLpoly{i, 1}(1, :) = A;
    PLpoly{i, 1}(2, :) = B;
    PLpoly{i, 1}(3, :) = C;

end


% centroid signs check to determine the interior sides of the planes.
for i = 1:1:size(PLpoly, 1)

    R = PLpoly{i, 1};
    [SOP] = check_point_side_of_plane(R, c);
    centroidsgns(i) = SOP;

end

for i = 1:1:size(PLpoly, 1)

    R = PLpoly{i, 1};
    [SOP] = check_point_side_of_plane(R, S);
    test(i) = SOP;

end

test = test .* centroidsgns;

if abs(sum(test)) == size(p, 1)
    aa = true;
else
    aa = false;
end


end