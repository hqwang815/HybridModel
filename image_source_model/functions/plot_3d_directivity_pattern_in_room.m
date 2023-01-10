
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

function [] = plot_3d_directivity_pattern_in_room(cartCoord, tf, fv, plotFrequency, figNum,fontSize,receiverPosition)


thisBin = find(fv >= plotFrequency, 1);

thisTf = tf(:, thisBin);
thisTfMagNorm = abs(thisTf) / max(abs(thisTf));
scaledCoord = bsxfun(@times, cartCoord, thisTfMagNorm);
scaleFactor = 0.5;
scaledCoord = scaledCoord * scaleFactor;

scaledCoord = scaledCoord + receiverPosition;

figure(figNum)

scatter3(scaledCoord(:, 1), scaledCoord(:, 2), scaledCoord(:, 3), '.k')


view(0, 90)

pbaspect([1, 1, 1])


set(gca,... 
    'YDir', 'normal',...
    'Fontsize',fontSize,...
    'fontname','Arial',...
    'XGrid','off',...
    'XMinorGrid','off',...
    'YGrid','off',...
    'YMinorGrid','off','YTickLabel',[],'XTickLabel',[],'ZTickLabel',[])
 
axis off



end