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
function [] = lay_out_figure(fontSize,lineWeight,pictureWidth,aspectRatio,nSubPlots,logXScale)


if logXScale==true
    scale='log';
else
    scale='lin';
end
 
set(gca,...
    'FontWeight','normal',...
    'XGrid','off',...
    'XMinorGrid','off',...
    'Fontsize',fontSize,...
    'XScale',scale,'YGrid','off',...
	'fontname','Arial','lineWidth',lineWeight)

a = nSubPlots;
set(gcf, 'units', 'centimeters', 'position', [0, 0, pictureWidth, (pictureWidth / aspectRatio(1) * aspectRatio(2) * a) + 1])

 end