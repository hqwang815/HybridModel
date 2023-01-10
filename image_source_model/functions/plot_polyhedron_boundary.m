%% function name
% Short description

%% Syntax
% # [output argument(s) ] = functionname(input arguments)

%% Description
% Description

% Input arguments
% * argument1                       -a 'filetype':	short description

% Output arguments
% * argument1                       -a 'filetype':	short description

%% Example(s)
% Example of function

%% Related functions



function []=plot_polyhedron_boundary(crnmat,fn)


figure(fn)

k = boundary(crnmat,0);
hold on
trisurf(k,crnmat(:,1),crnmat(:,2),crnmat(:,3),'Facecolor','blue','FaceAlpha',0.2,'EdgeColor',[1 1 1])
axis equal

end