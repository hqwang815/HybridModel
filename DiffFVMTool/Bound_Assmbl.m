function[Belements] = Bound_Assmbl(mesh, Belements, aCoef, thisBandNo)
% Boundary assembly function with already some data stored


% WW: runs through all boundary elements?

%% Loop to extract the absorption coefficients from the file names

nSurface = size(mesh.ENT, 1) - 2;
nBoundaryElement = mesh.NEB;


for iSurface = 1:nSurface

    thisSurfName = mesh.ENT{iSurface, 2}{1, 1};
    absorptionData = extractAfter(thisSurfName, '$');
    bandAbsorption = strsplit(absorptionData, ',');

    thisAlphaString = bandAbsorption(1, thisBandNo);
    thisAlpha = str2double(thisAlphaString);

    % Acoef is a function handle
    thisAcoef = aCoef(thisAlpha);
    surfNames{iSurface, 1} = thisSurfName;
    surfCoefs(iSurface) = thisAcoef;

end

surfaceNo=extractfield(Belements,'elset');

% Make an array with absorption coefficient for each boundary element

for iSurface = 1:nSurface
    
    thisAcoef = surfCoefs(iSurface);
    thisIndex = surfaceNo==iSurface; 
    aCoeffArray(thisIndex) = thisAcoef; 
    
end
    
x = num2cell(aCoeffArray');
% Assign the correct absorption coefficients at each boundary element
[Belements.Abn] = deal(x{:});



% for iBoundaryElement = 1:nBoundaryElement
% 
%     
%     thisSurfaceNo=Belements(iBoundaryElement).elset;
%     thisSurfaceName = mesh.ENT{thisSurfaceNo, 2}{1, 1};
% 
%     for iSurface = 1:1:nSurface
% 
%         thisSurfName = surfNames{iSurface, 1};
%         thisCoef = surfCoefs(iSurface);
% 
% 
%         surfaceIsSame = strcmp(thisSurfaceName, thisSurfName);
% 
% 
%         if surfaceIsSame == true
%             
%             Belements(iBoundaryElement).Abn = thisCoef;
% 
%         else
%     %else what
% 
%         end
% 
%     end
% 
% 
%     if mesh.ENT{Belements(iBoundaryElement).elset, 2}{1, 1} == "Facades"
%         thisAlpha = 0.01;
%         Belements(iBoundaryElement).Abn = aCoef(thisAlpha);
% 
%        
%     elseif isempty(Belements(iBoundaryElement).Abn)
%         text(0, 0.5, 'Err: No Corresponding BC material in database', 'fontsize', 16, 'color', 'red')
%         return
%     end
% end