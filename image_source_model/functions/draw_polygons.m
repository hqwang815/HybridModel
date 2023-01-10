function [] = draw_polygons(poly, figureNumber)


figure(figureNumber)

nPolygons = size(poly, 1);

for i = 1:1:nPolygons

    thisPolygon = poly{i, 1};

    nVertices = size(thisPolygon, 1);
    % markcolor=rand(1,3);
    markcolor = [0, 0, 0];


    polygonName = sprintf('PG:%d', i);
    centroidPolygon = sum(thisPolygon, 1) / nVertices;
    text(centroidPolygon(1), centroidPolygon(2), centroidPolygon(3), polygonName, 'FontSize', 9, 'Color', 'b')


    for j = 1:1:nVertices

        if j == nVertices
            p2 = thisPolygon(1, :);
        else
            p2 = thisPolygon(j+1, :);
        end
        p1 = thisPolygon(j, :);


%         str = sprintf('$--P_{%d}:[%0.2f,%0.2f,%0.2f]$', j-1,p1(1),p1(2),p1(3));
%         text(p1(1), p1(2), p1(3), str, 'Interpreter', 'Latex')
%         hold on


        line([p1(1), p2(1)], [p1(2), p2(2)], [p1(3), p2(3)], 'Color', markcolor, 'LineWidth', 2);
        hold on


    end


end