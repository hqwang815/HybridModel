% Description




function []=plot_reflections(pointsOfReflection,sourcePositions,thisReceiverPosition,reflectionSequences,nPlanes,nPolygons)



nSources=size(sourcePositions,1);
nOrder=size(pointsOfReflection,2);
hasPolygon=(nPlanes+2:nPlanes+1+nPolygons);




for iSource=1:1:nSources

	thisSourcePosition=sourcePositions(iSource,:);
	

	for jOrder=1:1:nOrder
	
	nRays=jOrder+1;
	nPointsOfReflection=size(pointsOfReflection{iSource,jOrder},1);


	for kPointOfReflection=1:1:nPointsOfReflection
	
	for lRays=1:1:nRays
		
		
		
		
		
switch lRays
    
	case 1
		p1=thisReceiverPosition;
		p2=pointsOfReflection{iSource,jOrder}{kPointOfReflection,lRays};
		
		thisWallNo=reflectionSequences{iSource,jOrder}(kPointOfReflection,lRays);
		
		check=thisWallNo==hasPolygon;
		if any(check)==true
			scatter3(p2(1),p2(2),p2(3),'ok')
			hold on
			scatter3(p2(1),p2(2),p2(3),'xr')
			hold on
		else
			scatter3(p2(1),p2(2),p2(3),'.g')
		end
		
    case nRays
        p1=pointsOfReflection{iSource,jOrder}{kPointOfReflection,lRays-1};
		p2=thisSourcePosition;
				
	otherwise
        p1=pointsOfReflection{iSource,jOrder}{kPointOfReflection,lRays-1};
		p2=pointsOfReflection{iSource,jOrder}{kPointOfReflection,lRays};
				
		
		thisWallNo=reflectionSequences{iSource,jOrder}(kPointOfReflection,lRays);
		
		check=thisWallNo==hasPolygon;
		if any(check)==true
			scatter3(p2(1),p2(2),p2(3),'ok')
			hold on
			scatter3(p2(1),p2(2),p2(3),'xr')
			hold on
		else
			scatter3(p2(1),p2(2),p2(3),'.g')
		end
		
		
end

plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'Color','b','LineStyle','--','LineWidth',0.5)




	end
	
	end
end





end