
%% Syntax
% # [Ism,Ism2,pointsOfReflections,cosinesAnglesOfIncidence,reflectionSequences]=determineimagesources(walls,sourcePositions,thisReceiverPosition,reflectionLevel)

%% Description
% Input arguments:

% # walls										- a cell-array containing 3-point description of each wall. 
% # sourcePositions					   - matrix with source positions (each row one source). 
% # thisReceiverPosition			- array with receiver position.
% # reflectionLevel						 - positive integer (order of reflection)

% Output arguments:

% # Ism											   - a structure containing the image source positions, 
%														   plane numbers and names after applying rule I
% # Ism2										  - a structure containing the image source positions, 
% 														  plane numbers and names after applying rule I & rule II
% # pointsOfReflections				   - a cell matrix containing all wall reflection positions
% # cosinesAnglesOfIncidence	 - a cell matrix containing all the cosines of wall incedences
% # reflectionSequences				  -	a cell matrix containing the reflection sequence for each image source

%%  Examples
% Give some examples.

%% Related functions
% Describe related functions


function [Ism,Ism2,pointsOfReflection,cosinesAnglesOfIncidence,reflectionSequences]=determine_image_sources(walls,sourcePositions,thisReceiverPosition,reflectionLevel)

% Determine the amount of sources.
nSources=size(sourcePositions,1);

% Determine the amount of planes.
nPlanes=size(walls,1);

% Create cell names
 sourcesNames=cell(nSources,1);
 for h=1:1:nSources
     sourcesNames{h,1}=sprintf('S%d',h);   
 end    

% Added the tree class as made by Jean-Yves Tinevez (http://tinevez.github.io/matlab-tree/). 
% Using a tree class instead of growing arrays making it easy to visualize the image source
% evolution. In version one of the code growing arrays where used using tracking arrays to determine the
% lineage of the image sources. 

% Preallocate array for signs of the directions of each plane towards an
% interior point.
planesDirections=zeros(1,nPlanes);

 for i=1:1:nPlanes
 %% Using the receiver position to check how the planes are oriented regarding the receiver position
 % This makes sure that the planes orientation is defined independent of how the planes are defined. 
 R=walls{i,1};    
 planesDirections(i)=check_point_side_of_plane(R,thisReceiverPosition);
 
 end
 
 planesDirections=planesDirections*-1;

 %% Check based on rule I (right side of plane check)

 [Ism]=check_rule_1(sourcePositions,walls,reflectionLevel,planesDirections,sourcesNames);
 


%% Check based on rule II (Image source lineage check)
% Ism2,nodeDepth,pointsOfReflection,cosinesAnglesOfIncidence,reflectionSequences
[Ism2,nodeDepth,pointsOfReflection,cosinesAnglesOfIncidence,reflectionSequences]=check_rule_2(nSources,sourcesNames,Ism,walls,thisReceiverPosition,planesDirections);


%% Remove image sources marked as ineffective.

[Ism2]= remove_invalid_image_sources(sourcesNames,Ism2,nodeDepth,nSources,reflectionLevel);


end