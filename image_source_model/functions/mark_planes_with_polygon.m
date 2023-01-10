


% planesWithPolygon,reflectionSequences


function [intersectionHasPolygon,polygonInPlane]=mark_planes_with_polygon(planesWithPolygon,reflectionSequences)

nLevels=size(reflectionSequences,2);

intersectionHasPolygon=cell(1,nLevels);
polygonInPlane=cell(1,nLevels);

for iLevel=1:1:nLevels
	iLevelreflectionSequences=reflectionSequences{1,iLevel};
	[planeHasPolygon,planeIndeces]=ismember(iLevelreflectionSequences,planesWithPolygon);

	intersectionHasPolygon{1,iLevel}=planeHasPolygon;
	polygonInPlane{1,iLevel}=planeIndeces;
	
end


end