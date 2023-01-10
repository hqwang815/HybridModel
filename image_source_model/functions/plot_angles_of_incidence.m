%% Syntax
% # [output]=functionName(input)
%% Description
% Input arguments:
% # variable name  comment
% Output arguments:
% # variable name   comment
%%  Examples
% Example
%% Related functions %Related functions
%%
function [] = plot_angles_of_incidence(anglesOfIncidence,grid,thisPosition)

scale=0.25;
% 
nLevels=size(anglesOfIncidence,2);
% gridAtRec=grid*scale+thisPosition;
% scatter3(gridAtRec(:,1),gridAtRec(:,2),gridAtRec(:,3),'MarkerEdgeColor',[0.5,0.5,0.5],'Marker','.')
%  hold on

    for iLevel=1:nLevels
    thisIncidenceArray=anglesOfIncidence{1,iLevel}*(scale*1.05);
    thisColor=rand(3,1);
    thisIncidenceArray = thisIncidenceArray+thisPosition;
    scatter3(thisIncidenceArray(:,1),thisIncidenceArray(:,2),thisIncidenceArray(:,3),'MarkerFaceColor',thisColor,'MarkerEdgeColor','w','Marker','d')
   hold on
    end
    
 end