% sourcesNames,Ism2,nodeDepth,nSources,reflectionLevel


function [Ism2]= remove_invalid_image_sources(sourcesNames,Ism2,nodeDepth,nSources,reflectionLevel)


for i=1:1:nSources
    
	
	% Preallocate cell arrays
    checkPositionsAll=cell(1,reflectionLevel);
    checkPlaneNoAll=cell(1,reflectionLevel);
    checkImageSourceNamesAll=cell(1,reflectionLevel);





    for j=1:1:reflectionLevel
        
        imageSourcesInLevel=find(nodeDepth==j);
        imageSourcePositions=Ism2.(sourcesNames{i}).isp.Node(imageSourcesInLevel);
        imageSourcesPlaneNo=Ism2.(sourcesNames{i}).isn.Node(imageSourcesInLevel);
        imageSourcesNames=Ism2.(sourcesNames{i}).isnms.Node(imageSourcesInLevel);
    
		nImageSources=length(imageSourcesInLevel);
		
        % Preallocate arrays
        checkPositions=zeros(nImageSources,3);
        checkPlaneNo=zeros(nImageSources,1);
        checkImageSourceNames=cell(nImageSources,1);   
    
        for k=1:1:size(imageSourcePositions,1)
        
            checkPositions(k,:)=imageSourcePositions{k,1};
            checkPlaneNo(k,1)=imageSourcesPlaneNo{k,1};
            checkImageSourceNames{k,1}=imageSourcesNames{k,1};
        
        end
        
        checkPositions( ~any(checkPositions,2), : ) = [];
        checkPositionsAll{1,j}=checkPositions;
        I=checkPlaneNo~=0;
        checkPlaneNo(checkPlaneNo==0) = [];
        checkPlaneNoAll{1,j}=checkPlaneNo;
        checkImageSourceNames=checkImageSourceNames(I);
        checkImageSourceNamesAll{1,j}=checkImageSourceNames;
     
        clear checkPositions
    

    end

   Ism2.ispclean(i,:)=checkPositionsAll;
   Ism2.isnclean(i,:)=checkPlaneNoAll;
   Ism2.isnmsclean(i,:)=checkImageSourceNamesAll; 
  
end



end