
% sourcePositions,walls,reflectionLevel,planesDirections,sourcesNames





function [Ism]=check_rule_1(sourcePositions,walls,reflectionLevel,planesDirections,sourcesNames)

 nSources=size(sourcePositions,1);
 nWalls=size(walls,1);
 
for i=1:1:nSources 
    % Setup tree to hold plane numbers
    t=tree(i);  
    p4=sourcePositions(i,:);
    % Setup tree to hold image source position vectors
    ist=tree(p4);
    namestree=tree(i);

    
   for j=1:1:reflectionLevel
       
       % Make a depth tree this indicates the depth (level) of each node in
       % the tree.
      
       dt = t.depthtree;
              
       % Find all nodes in the parent level
       parentNodes=find(dt==j-1);

       % Setup a loop with length equal to the amount of nodes in parent level.       
       for k=1:1:length(parentNodes)
            thisParentNode=parentNodes(k);
            
            %take the position of the kth node image source of the parent level
             thisImageSource=ist.Node{thisParentNode,1};
       
            % Loop through the walls of the enclosure
            for l=1:1:nWalls
                
                % Take three points representing the l-th plane 
                R=walls{l,1};   
                % Function to check if (Image) Source is on the interior side of the reflecting wall
                [sideOfPlane]=check_point_side_of_plane(R,thisImageSource);

                % Multiply by vector designating interior sides of plane
                sideOfPlaneIntExt=sideOfPlane*planesDirections(l);
                
                % If chk2 returns a -1 perform reflection in the wall. 
                % When chk2 returns a 1 (image) source is on the exterior 
                % side of the given plane and no reflection should be performed. 

                if sideOfPlaneIntExt==-1 
                % Reflect the (Image) source in wall nr: d. p5 is the Image source.
   
                [reflectedImageSource]=calculate_mirror_reflection(R,thisImageSource);    
                
                % Add the new image source position to the tree holding th
                % e image source position vectors
                
                ist=ist.addnode(thisParentNode,reflectedImageSource);
                
                % Add the plane nr (l) to the plane number tree
                t=t.addnode(thisParentNode,l);

                namestree=namestree.addnode(thisParentNode,0);
                
                
                end

            end
       
       end
   
   end

   Ism.(sourcesNames{i}).isn=t;
   Ism.(sourcesNames{i}).isp=ist;
   Ism.(sourcesNames{i}).isnms=namestree;   

end


end