
% sourcePositions,pointsOfReflection,thisReceiverPosition




function [angleOfIncidenceSource,angleOfIncidenceReceiver]=calculate_angles_of_incidence(sourcePositions,pointsOfReflection,thisReceiverPosition)
% This function uses the intersections calculated with the reverse
% raytracing function to determine the angles of incidence at the receiver
% and at the source
[mSources,nLevels]=size(pointsOfReflection);
angleOfIncidenceSource=cell(mSources,nLevels);
angleOfIncidenceReceiver=cell(mSources,nLevels);



for i=1:1:mSources
	
	for j=1:1:nLevels+1     
        
        
        % For direct sound (angle of incidence = norm of position source - position receiver)
		if j==nLevels+1
            
            point2=thisReceiverPosition;
            % source
            vSource=point2-sourcePositions(i,:);
            % normalize to unit sphere. angle of incidence is now expressed
            % as a point on the unit sphere
            
            % receiver
            vReceiver=sourcePositions(i,:)-point2;
            
            % Angle of incidence source
            angleOfIncidenceSource{i,j}(1,:)=vSource/norm(vSource);
            % Angle of incidence receiver
            angleOfIncidenceReceiver{i,j}(1,:)=vReceiver/norm(vReceiver);
        
            
        % For reflected sound    
        else    
            [ni,lp]=size(pointsOfReflection{i,j});
			for k=1:1:ni      
            
            % take the first sound ray wall intersection from the source 
            point2=pointsOfReflection{i,j}{k,lp};
            % take the last sound ray wall intersection to the receiver 
            poir=pointsOfReflection{i,j}{k,1};
            % take the difference between the position vector of the wall
            % intersection and the position vector of the source
            vSource=point2-sourcePositions(i,:);
            
            vReceiver=poir-thisReceiverPosition;
                       
            % normalize source angle of incidence to unit sphere. angle of incidence is now expressed
            % as a point on the unit sphere
            angleOfIncidenceSource{i,j}(k,:)=vSource/norm(vSource);
            % normalize receiver angle of incidence to unit sphere. angle of incidence is now expressed
            % as a point on the unit sphere
            angleOfIncidenceReceiver{i,j}(k,:)=vReceiver/norm(vReceiver);
                                                
			end
			
		end
		
	end
	
end

end