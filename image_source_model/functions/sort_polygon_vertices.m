% % General Discription

% This functions can be used to order the vertices of a convex polygon
% in clockwise or counterclockwise succession. In the image source model this is
% important for defining the orthogonal planes to the polygon. Example if
% there are 6 vertices in a polygon then then there are a 720 ways to
% string these vertices together, the shortest path is the correct one. there
% are in this case 12 possible correct cases starting at either vertices (6 pieces)
% moving in counter clockwise or clockwise direction (2 directions).
% making 6 * 2 possibilities.


function [sortedPolygonVertices]=sort_polygon_vertices(polygonVertices) 

nVertices=size(polygonVertices,1);

pc=(1:1:nVertices);
permutations = perms(pc);

[m,n]=size(permutations);

distances=zeros(m,1);
pathLengths=zeros(n,1);

for i=1:1:m
	
	
	for j=1:1:n
		
		p1=polygonVertices(permutations(i,j),:);
		
		if j==size(permutations,2)
		p2=polygonVertices(permutations(i,1),:);	
		else
		p2=polygonVertices(permutations(i,j+1),:);
		end
		pathLengths(j)=norm(abs(p1-p2));
		
	end
	distances(i,:)=sum(pathLengths);
	
end



[~,idx]=min(distances);
rr=permutations(idx,:);
sortedPolygonVertices=polygonVertices(rr,:);

end



