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
function [] = set_figure(pictureWidth,aspectRatio,nSubPlots,nthOctave)



fontSize = 10;



fCentre = 15.625*2.^((1/nthOctave)*(0:1:10*nthOctave));
faChar = num2str(fCentre,'% 5.0f ');
faChar2 = strsplit(faChar);
lineWeight = 1;


set(gca,...
    'FontWeight','normal',...
    'XGrid','off',...
    'XMinorGrid','off',...
    'Fontsize',fontSize,...
    'XTick',fCentre,...
    'XTickLabel',faChar2,...
    'XScale','log','YGrid','off',...
	'fontname','Arial','lineWidth',lineWeight)

a = nSubPlots;
set(gcf, 'units', 'centimeters', 'position', [0, 0, pictureWidth, (pictureWidth / aspectRatio(1) * aspectRatio(2) * a) + 1])




 end