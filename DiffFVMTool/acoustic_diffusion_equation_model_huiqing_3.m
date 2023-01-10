% Code developped by Quentin Goestchel for the analysis of the FVM method
% on diffusion equation.
%--- CAUTION !---% The information below are important:
% To run the code below, the file newid.m must be added in the directory
% MATLAB_<version>/toolbox/matlab/uitools and the gptoolbox folder must be
% in the current folder.
% For Mac (Unix) user, Matlab must be opened from the
% Terminal in order to enable the command line control of Gmsh
% For Windows user, Gmsh.exe must be in the work folder
% To download Gmsh: http://gmsh.info/

clear all
close all

%% Settings


% General settings

% Sound particle velocity [m.s^-1]
c0 = 343;
% Air density [Kg.m^-3] at 20°C
rho = 1.21;
% Reference pressure
pRef = 2 * 10^-5;

% Frequency resolution
fcLow = 125;
fcHigh = 2000;
nthOctave = 1;

% Point power source [W]
% Source power dB
x = 0;

Ws = (10^(x / 20) * pRef / (sqrt(rho) * c0))^2;

% Ws = 2.809870992500524e-15;

% Ws = 10^-2;

%Scenario
scenario = 2;
scenarioName = sprintf('scenario%0.0f_barewall_abs_01', scenario);
% scenarioName = 'valuelongroom';


%Boundary condition
boundaryCondition = 'Sabine';

%Mesh length
meshLength = 0.5;
% Transducer positioning
% Source position 1 (1500 mm)
% sourcePosition = [3.01, 2.64, 1.62];


% sourcePosition = [2, 2, 2];


% Source position 2 (3500 mm)
sourcePosition = [1.36, 3.76, 1.62];

% Receiver
% receiverPosition = [4.26, 1.76, 1.62];

receiverPositions(1, :) = [0.5, 0.5, 1.62];


distance_receiver=norm(sourcePosition-receiverPositions)
% receiverPositions(2, :) = [2.03, 0.5, 1.62];
% receiverPositions(3, :) = [3.56, 0.5, 1.62];
% receiverPositions(4, :) = [5.09, 0.5, 1.62];
% receiverPositions(5, :) = [0.5, 1.83, 1.62];
% receiverPositions(6, :) = [2.09, 1.74, 1.62];
% receiverPositions(7, :) = [3.68, 1.64, 1.62];
% receiverPositions(8, :) = [5.28, 1.54, 1.62];
% receiverPositions(9, :) = [0.5, 3.17, 1.62];
% receiverPositions(10, :) = [2.15, 2.97, 1.62];
% receiverPositions(11, :) = [3.80, 2.77, 1.62];
% receiverPositions(12, :) = [5.46, 2.57, 1.62];
% receiverPositions(13, :) = [0.5, 4.5, 1.62];
% receiverPositions(14, :) = [2.21, 4.21, 1.62];
% receiverPositions(15, :) = [3.92, 3.91, 1.62];
% receiverPositions(16, :) = [5.64, 3.61, 1.62];
% receiverPositions(17, :) = [4.26 ,1.76, 1.62];



% receiverPosition(1,:) = sourcePosition + [3, 0, 0];
% receiverPosition(2,:) = sourcePosition + [13, 0, 0];
% receiverPosition(3,:) = sourcePosition + [23, 0, 0];

% Plot settings
pictureWidth = 20;
fontSize = 9;
aspectRatio = [5, 3, 1];

delta = 0.1;
distance = norm(sourcePosition-mean(receiverPositions,1));
distanceMm = round(distance, 1) * 10^3;

%% Pre-processing

nBands = nthOctave * log(fcHigh/fcLow) / log(2) + 1;


% Create center frequency array
fMin = (2^(1 / 2) * 2^(1 / (2 * nthOctave)) * fcLow) / 2;
fMax = (2^(1 / 2) * fcHigh) / 2^(1 / (2 * nthOctave));

x = log(fMax/fMin) / log(2) * nthOctave;
centerFrequencies = fMin * 2.^((0:x) / nthOctave);
centerFrequencies = round(centerFrequencies);

%
spl = @(input) (20 * log10(abs(input)));

%%

% Add working folder to the path
% cd 'D:\OneDrive - TU Eindhoven\Documents\phd\matlab\code_from_others\acoustic_simulation\de_quentin\DiffFVMTool'

% Import the gptoolbox, download it and change the path
gp_subdirs = split(genpath(strcat(pwd, '/gptoolbox/')), ';');
addpath(strjoin(gp_subdirs(~contains(gp_subdirs, '.git')), ';'));
savepath
% Configuration of the plot interface
% set(gcf, 'Position', get(0, 'Screensize'));
ax = subplot(2, 2, 1);
set(ax, 'visible', 'off')
% figure(1)
% subplot(2, 2, 1)
% text(0, 1, 'Initialising modelisation...', 'fontsize', fontSize);

%% Import Mesh

% Manual input of the filename, characteristic length of the mesh,...
% Boundary conditions, importations, recording time and source emition duration
prompt = {'Enter the name of the .geo file [Without extension]: ', ...
    'Enter the characteristic length of the mesh [in meter]: ', ...
    'Enter the boundary conditions type desired [1:Sabine,2:Eyring,3:Modified] ', ...
    'Is the mesh already imported? [Y/N](default value Y): ', ...
    'Do you want a 3D plot of the mesh orthogonality? [Y/N](default value N): ', ...
    'Enter recording time [in seconds]: '};
defaultAnswer = {scenarioName, num2str(meshLength), '1', 'N', 'Y', '8'};
numLines = 1;
answer = newid(prompt, 'Code inputs', numLines, defaultAnswer);
options.Resize = 'on';
filename = num2str(answer{1});
lengthOfMesh = num2str(answer{2});
boundaryConditionType = str2double(answer{3});
meshIsImported = num2str(answer{4});
plotMeshOrthogonality = num2str(answer{5});
recordingTime = str2double(answer{6});
% sourceOnTime = str2double(answer{7});

% To save time, the mesh importation can be stored in mesh_data.mat
if isempty(meshIsImported)
    str = 'Y';
end
if meshIsImported == 'N' % Condition if nothing has been stored before:
    % The lines below command Gmsh to do the 3D meshing
    if lengthOfMesh == 'H'
        strimp = strcat({'gmsh -3 -format msh2'}, {pwd}, {'/GeoModels/'}, {filename}, {'.geo'});
        system(strimp{1});
        lengthOfMesh = 2;
    else
        strimp = strcat({'gmsh -3 -format msh2 -clscale '}, {lengthOfMesh}, {' '}, {'"'}, {pwd}, {'"'}, {'/GeoModels/'}, {filename}, {'.geo'});
        system(strimp{1});
        lengthOfMesh = str2double(lengthOfMesh);

    end
    % The lines below open the .msh file, call the connectivity function
    % and save all the data in a .mat file
    msh_file = strcat(pwd, '/GeoModels/', filename, '.msh');
    %         mesh_T = load_gmsh_WHQ(msh_file);
    [mesh, Velements, Belements, Cell_center, Vcell] = Input_mesh(msh_file, lengthOfMesh);
    % Saving of the mesh connectivity data:
    save('mesh_data.mat', 'mesh', 'Velements', 'Belements', 'Cell_center', 'Vcell')
else
    load('mesh_data.mat');
    lengthOfMesh = str2double(lengthOfMesh);
end

figure(1)
subplot(2, 2, 1)
text(0, 0.9, 'Volume meshed by Gmsh', 'fontsize', fontSize);
text(0, 0.8, 'Connectivity process completed', 'fontsize', fontSize);

nElementAll = mesh.NEV;
nElementBoundary = mesh.NEB;

%% Mesh orthogonality evaluation

% Call of the orthogonnality index computation function
[orthogonalityAllElements, ~, orthogonalityInteriorElements, ~, orthogonalityBoundaryElements] = Orth_idx(mesh, Velements, Cell_center);

% Sorting of the axis and data for volume and boundary related elements
lowerLimit = (0:delta:1 - delta);
upperLimit = (delta:delta:1);
nLimit = length(lowerLimit);

% Generate xlabels for sorted orthogonality indeces
for iLimit = 1:nLimit
    thisLowerLimit = lowerLimit(iLimit);
    thisUpperLimit = upperLimit(iLimit);
    xLabels{iLimit, 1} = sprintf('[%0.1f,%0.1f]', thisLowerLimit, thisUpperLimit);
end


% orthogonalityInteriorElements = sort(Chi_elV);
% orthogonalityBoundaryElements = sort(Chi_elB);

nInteriorElements = length(orthogonalityInteriorElements);
nBoundaryElement = length(orthogonalityBoundaryElements);

for iLimit = 1:nLimit

    thisIndexVol = orthogonalityInteriorElements >= lowerLimit(iLimit) & orthogonalityInteriorElements < upperLimit(iLimit);
    thisIndexBnd = orthogonalityBoundaryElements >= lowerLimit(iLimit) & orthogonalityBoundaryElements < upperLimit(iLimit);
    nElementVol = sum(thisIndexVol);
    nElementBnd = sum(thisIndexBnd);
    orthoVolume(iLimit) = nElementVol / nInteriorElements * 100;
    orthoBoundary(iLimit) = nElementBnd / nBoundaryElement * 100;

end

% Y = [Chi_elV; Chi_elB];
figure(1)
subplot(2, 2, 2)


xV = (0:1:nInteriorElements - 1);
xB = (0:1:nBoundaryElement - 1);

bar(orthoVolume, 'k')

pbaspect([aspectRatio(1), aspectRatio(2), aspectRatio(3)])
xtickangle(45)
set(gca, 'XTickLabel', xLabels, 'FontSize', fontSize)
ylim([0, 100])


title('Elements (interior)', 'fontsize', fontSize)
xlabel('Mesh orthogonal quality', 'fontsize', fontSize)
ylabel('Percentage of elements', 'fontsize', fontSize)

% legend({'Volume elements', 'Boundary related elements'}, 'FontSize', fontSize)


subplot(2, 2, 4)
bar(orthoBoundary, 'r')

pbaspect([aspectRatio(1), aspectRatio(2), aspectRatio(3)])
xtickangle(45)
set(gca, 'XTickLabel', xLabels, 'FontSize', fontSize)
ylim([0, 100])


title('Elements (boundary)', 'fontsize', fontSize)
xlabel('Mesh orthogonal quality', 'fontsize', fontSize)
ylabel('Percentage of elements', 'fontsize', fontSize)


if isempty(plotMeshOrthogonality)
    str = 'N';
end
if plotMeshOrthogonality == 'Y'
    T = zeros(nElementAll, 4);
    for i = 1:nElementAll
        T(i, :) = Velements(i).nodes;
    end
    figure(1)
    subplot(2, 2, 3)
    tetramesh(T, mesh.POS, orthogonalityAllElements, 'FaceAlpha', 0.4)
    title('Orthogonal quality of the mesh', 'fontsize', fontSize)
    colormap winter
    light;
    lighting phong;
    camlight('left');
    axis equal
    caxis([0, 1]) %Set the scale in order to have an interpretable graphic result
    scale = colorbar; % Create a colorbar
    scale.Label.String = 'Orthogonality index';
    scale.Label.FontSize = fontSize;
    view(0, 90)
    drawnow

end


figure(1)
subplot(2, 2, 1)


% mean and median of orthogonality indeces
meanVolume = mean(orthogonalityInteriorElements);
medianVolume = median(orthogonalityInteriorElements);

meanBoundary = mean(orthogonalityBoundaryElements);
medianBoundary = median(orthogonalityBoundaryElements);

meanAll = mean(orthogonalityAllElements);
medianAll = median(orthogonalityAllElements);

text1 = sprintf('(mean ; median) = (%0.02f ; %0.02f) Elements (interior)', meanVolume, medianVolume);
text2 = sprintf('(mean ; median) = (%0.02f ; %0.02f) Elements (boundary)', meanBoundary, medianBoundary);
text3 = sprintf('(mean ; median) = (%0.02f ; %0.02f) Elements (All)', meanAll, medianAll);

text(0, 0.7, text1, 'fontsize', fontSize);
text(0, 0.6, text2, 'fontsize', fontSize);
text(0, 0.5, text3, 'fontsize', fontSize);

a = 1;
set(gcf, 'units', 'centimeters', 'position', [0, 0, pictureWidth, (pictureWidth / aspectRatio(1) * aspectRatio(2) * a) + 1])

%% Pre-processing


% Below, temporal discretization (time step) [s] which is set to fulfill
% Adapted empirical Navarro criterion
dt = sqrt(lengthOfMesh/2*(1 * 10^-8));
% fm = 1 / dt;

% Extreme limit of the consistance condition:
if lengthOfMesh / (2 * dt) < c0
    figure(1)
    subplot(2, 2, 1)
    text(0, 0.3, 'Error: lc/(2dt) < to c.', 'fontsize', fontSize, 'color', 'red')
    text(0, 0.2, 'The propagation of the numerical scheme is slower than sound speed', 'fontsize', fontSize);
    return
end

figure(1)
subplot(2, 2, 1)
text(0, 0.4, ['"Adapted Navarro" empirical criterion = ', num2str(2*dt^2/lengthOfMesh)], 'fontsize', fontSize');

%% Set up anonymous functions

% Absorption coefficients using the  1. sabine A, 2. eyring A, and 3. modified Xiang A

if boundaryConditionType == 1
    aCoef = @(alpha) (c0 * alpha) / 4;
elseif boundaryConditionType == 2
    aCoef = @(alpha) (c0 * (-log(1-alpha))) / 4;
elseif boundaryConditionType == 3
    aCoef = @(alpha) (c0 * alpha) / (2 * (2 - alpha));
else
    close all hidden
    disp('Error : Wrong boundary condition input')
    return
end

% Create anonymous function to calculate area of tetrahedron face
triangleArea = @(v1, v2)((1 / 2) * norm(cross(v1, v2)));
% Create anonymous function to calculate volume of a tetrahedral element
tetrahedralVolume = @(v1, v2, v3)(norm(dot(v1, cross(v2, v3))) / 6);

% Geometrical parameters for diffusive fluxes
% Sum of the 4 connectivity face area over the distance to the neighbour per
% element:
fel = zeros(nElementAll, 1); % For the volume elements

% Global matrix of the faces areas over the distance to the neighbour:
%     Fmat = sparse(nElementVolume, nElementVolume);
Fmat = zeros(nElementAll, nElementAll);

%% Run difussion equation model for each band

for iBand = 1:nBands

    thisBandNo = iBand;
    thisFc = centerFrequencies(iBand);

    % freq is the (center) frequency that DE model is solving
    % So the source power increases with frequency.

    % Formula derived by Huiqing
%         Qs = 1 / (2 * pi * thisFc * rho);
%         Ws = Qs^2 * (2 * pi * thisFc)^2 * rho / (8 * pi * c0);
    % Simplified formula below
        Ws = 1/(8*pi*c0*rho);

    %% Boundary condition assembly %%

    % th type of constant A on boundary conditions. options:
    % Condition for meshes without any BC data (WW do I really need this?)

    if mesh.ENT{1, 1} == 0
        figure(1)
        subplot(3, 2, 1)
        alpha = 0.1;
        str = sprintf('No mesh BC absorption data: BC are set to alpha = %0.1f', alpha);
        text(0, 0.6, {str}, 'fontsize', 14);
        Abn = aCoef(alpha);
    else

        % All the boundary assembly process is made in the following function
        % To add a material, open the source code of the Bound_Assmbl function.


        Belements = Bound_Assmbl(mesh, Belements, aCoef, thisBandNo);

    end

    %% Geometrical parameters computation %%
    % Total boundary surface
    boundaryAreaTotal = 0;
    %     S = 0;
    %     Serr = 0;

    %     fb = zeros(nElementAll, 1); % Boundary vector (For the boundary faces)

    %     dtMat = sparse(nElementVolume, nElementVolume);
    % Fmat = zeros(mesh.NEV,mesh.NEV);
    %     wtBC = waitbar(0, 'Attribution of the boundary conditions');
    %     BCcount = 0;


    % Extract the tetrahedron neighbour indeces from the Velements structure
    neighbourVolume = extractfield(Velements, 'neighbours');
    neighbourVolume = reshape(neighbourVolume, 4, nElementAll)';


    % These are the faces not at the boundary
    % % Index the boundary elements
    nonBoundaryIndex = neighbourVolume ~= 0;

    % Extract the tetrahedron nodes indeces from the Velements structure
    nodesVolume = extractfield(Velements, 'nodes');
    nodesVolume = reshape(nodesVolume, 4, nElementAll)';


    nodePositions = mesh.POS;
    nodesBoundary = extractfield(Belements, 'nodes');
    nodesBoundary = reshape(nodesBoundary, 3, nElementBoundary)';


    %     combinations = nchoosek((1:4), 3);
    % The faces of a tetrahedral element are the four possible combinations
    % of 3 of the 4 nodes
    combinations = [1, 2, 3; 2, 3, 4; 1, 2, 4; 1, 3, 4];

    nCombination = size(combinations, 1);
    nNode = size(combinations, 2);

    % Extract all possible combinations of nodes. There are 4 possible
    % combinations of three nodes. a tetrahedral volume can have a face on the boundary
    % that can consist of any combination of three nodes taken from the set of 4 of that tetrahedron.

    for iCombination = 1:nCombination

        for iNode = 1:nNode

            thisIndex = combinations(iCombination, iNode);
            thisVolumeArray = nodesVolume(1:end, thisIndex:thisIndex);
            nodesVolumeComb{iCombination, 1}(:, iNode) = thisVolumeArray;

        end

    end


    % Need to also do this part for the elements not on the boundary
    nodesVolume1 = nodesVolumeComb{1, 1};
    nodesVolume2 = nodesVolumeComb{2, 1};
    nodesVolume3 = nodesVolumeComb{3, 1};
    nodesVolume4 = nodesVolumeComb{4, 1};


    nodesVolumeSorted1 = sort(nodesVolume1')';
    nodesVolumeSorted2 = sort(nodesVolume2')';
    nodesVolumeSorted3 = sort(nodesVolume3')';
    nodesVolumeSorted4 = sort(nodesVolume4')';

    nodesBoundarySorted = sort(nodesBoundary')';

    % Boundary elements in a corner can share the same tetrahedron
    [elementIsOnFace(:, 1), elementIndexBoundary(:, 1)] = ismember(nodesBoundarySorted, nodesVolumeSorted1, 'rows');
    [elementIsOnFace(:, 2), elementIndexBoundary(:, 2)] = ismember(nodesBoundarySorted, nodesVolumeSorted2, 'rows');
    [elementIsOnFace(:, 3), elementIndexBoundary(:, 3)] = ismember(nodesBoundarySorted, nodesVolumeSorted3, 'rows');
    [elementIsOnFace(:, 4), elementIndexBoundary(:, 4)] = ismember(nodesBoundarySorted, nodesVolumeSorted4, 'rows');

    elementArrayBoundary = sum(elementIndexBoundary, 2);

    nodesElementInBoundary = nodesVolume(elementArrayBoundary, :);

    %%

    % Extract the cartesian coordinates of each node for all tetrahedral
    % elements that have a face on the boundary

    for iElement = 1:nElementBoundary

        thisIndex1 = nodesElementInBoundary(iElement, 1);
        thisIndex2 = nodesElementInBoundary(iElement, 2);
        thisIndex3 = nodesElementInBoundary(iElement, 3);
        thisIndex4 = nodesElementInBoundary(iElement, 4);

        volumeNodeInBoundaryPos1(iElement, :) = nodePositions(thisIndex1, :);
        volumeNodeInBoundaryPos2(iElement, :) = nodePositions(thisIndex2, :);
        volumeNodeInBoundaryPos3(iElement, :) = nodePositions(thisIndex3, :);
        volumeNodeInBoundaryPos4(iElement, :) = nodePositions(thisIndex4, :);

    end

    %% The elements on the boundary


    faceNo = (1:4);

    % Allocate boundary vector
    boundaryV = zeros(nElementAll, 1);

    nodeOrder = [3, 1, 2, 1; 4, 2, 3, 2; 4, 1, 2, 1; 4, 1, 3, 1];


    for iElement = 1:nElementBoundary

        thisAbsorption = Belements(iElement).Abn;

        % Positions of the 4 nodes per tetrahedral element in cartesian coordinates
        thisNodePos(1, :) = volumeNodeInBoundaryPos1(iElement, :);
        thisNodePos(2, :) = volumeNodeInBoundaryPos2(iElement, :);
        thisNodePos(3, :) = volumeNodeInBoundaryPos3(iElement, :);
        thisNodePos(4, :) = volumeNodeInBoundaryPos4(iElement, :);

        % Find the face on which the boundary element lies
        [thisElementIsOnFace] = elementIsOnFace(iElement, :) == 1;

        % Select the right face number
        thisFaceNo = faceNo(thisElementIsOnFace);

        % Index of boundary element in the set of all elements
        thisVolumeElementNo = elementArrayBoundary(iElement);

        % Assign the absorption of the boundary element to a variable
        thisAbsorption = Belements(iElement).Abn;


        thisNodeOrder = nodeOrder(thisFaceNo, :);

        thisNodePos1 = thisNodePos(thisNodeOrder(1), :);
        thisNodePos2 = thisNodePos(thisNodeOrder(2), :);
        thisNodePos3 = thisNodePos(thisNodeOrder(3), :);
        thisNodePos4 = thisNodePos(thisNodeOrder(4), :);

        v1 = thisNodePos1 - thisNodePos2;
        v2 = thisNodePos3 - thisNodePos4;

        thisBoundaryArea = triangleArea(v1, v2);
        thisTotalAbsorption = thisBoundaryArea * thisAbsorption;
        boundaryAreaTotal = boundaryAreaTotal + thisBoundaryArea;

        % note: Boundary nodes can have absorption on different sides (e.g. in a corner), the
        % total absorption is summed per tetrahedral element !

        boundaryV(thisVolumeElementNo) = boundaryV(thisVolumeElementNo) + thisTotalAbsorption;


    end

    S = boundaryAreaTotal;

    %% Elements in the interior


    for iElement = 1:nElementAll

        thisIndex1 = nodesVolume(iElement, 1);
        thisIndex2 = nodesVolume(iElement, 2);
        thisIndex3 = nodesVolume(iElement, 3);
        thisIndex4 = nodesVolume(iElement, 4);

        volumeNodePos1(iElement, :) = nodePositions(thisIndex1, :);
        volumeNodePos2(iElement, :) = nodePositions(thisIndex2, :);
        volumeNodePos3(iElement, :) = nodePositions(thisIndex3, :);
        volumeNodePos4(iElement, :) = nodePositions(thisIndex4, :);

    end

    interior = zeros(nElementAll, nElementAll);


    nodeNo = (1:4);

    for iElement = 1:nElementAll

        % Positions of the 4 nodes per tetrahedral element in cartesian coordinates
        thisElementPos(1, :) = volumeNodePos1(iElement, :);
        thisElementPos(2, :) = volumeNodePos2(iElement, :);
        thisElementPos(3, :) = volumeNodePos3(iElement, :);
        thisElementPos(4, :) = volumeNodePos4(iElement, :);


        thisCellCenter = Cell_center(iElement, :);
        thisBoundaryIndex = nonBoundaryIndex(iElement, :);
        nNeighbour = sum(thisBoundaryIndex);
        nodeArray = nodeNo(thisBoundaryIndex);

        % Run through all the neighbouring elements
        for iNeighbour = 1:nNeighbour

            thisNeighbourNo = nodeArray(iNeighbour);

            thisNodeOrder = nodeOrder(thisNeighbourNo, :);
            thisNeighbour = Velements(iElement).neighbours(thisNeighbourNo);
            thisCellCenterNeighbour = Cell_center(thisNeighbour, :);

            v1 = thisElementPos(thisNodeOrder(1), :) - thisElementPos(thisNodeOrder(2), :);
            v2 = thisElementPos(thisNodeOrder(3), :) - thisElementPos(thisNodeOrder(4), :);

            thisFaceArea = triangleArea(v1, v2);
            thisDistance = norm(thisCellCenter-thisCellCenterNeighbour);

            % Area of the face between two adjacent nodes normalized over the distance between these nodes.
            interior(iElement, thisNeighbour) = thisFaceArea / thisDistance;

        end


    end

    Fmat = sparse(interior);
    fb = boundaryV;
    fel = sum(interior, 2);

    %% Quinten

    spoint = sourcePosition;


    % W.W.

    % Find closest cell center to receiver position

    nReceivers = size(receiverPositions, 1);


    for iReceiver = 1:nReceivers

        thisRecPos = receiverPositions(iReceiver, :);
        [~, recIndex] = min(vecnorm(bsxfun(@minus, Cell_center, thisRecPos), 2, 2));
        % Set this closest position as receiver position
        Rpel(iReceiver) = recIndex;

    end


    % Source element list initialisation for future interpolation
    sel_list = [];
    lsum = 0;
    li = [];
    Vs = 0;

    nodePositions = mesh.POS;
    distanceVectorSource = bsxfun(@minus, nodePositions, spoint);


    nodeToSourceDistance = vecnorm(distanceVectorSource');
    [~, sourcePosIndex] = min(nodeToSourceDistance);


    %     thisSourcePos = mesh.POS(sourcePosIndex,:);
    %     thisRecPos = Cell_center(recIndex,:);
    %
    %     figure(11)
    %     scatter3(thisSourcePos(1),thisSourcePos(2),thisSourcePos(3),'ko')
    %     hold on
    %     scatter3(thisRecPos(1),thisRecPos(2),thisRecPos(3),'ro')
    %     norm(thisSourcePos-thisRecPos)


    for i = 1:nElementAll

        if ismember(sourcePosIndex, Velements(i).nodes)

            sel_list = [sel_list; i]; %#ok<*AGROW> WW this is a bit weird
            lsum = lsum + norm(Cell_center(i)-mesh.POS(sourcePosIndex, :));
            % distance from the source node
            li = [li; norm(Cell_center(i)-mesh.POS(sourcePosIndex, :))];
            Vs = Vs + Vcell(i);
        end
    end

    % [W/m^3]
    w1 = Ws / Vs;

    %% WW Check volumes

    %     % Dit lijkt te kloppen
    %
    %     thisElement = 25;
    %     % Take the nodes that make a tetrahedral element
    %     thisNodeArray = Velements(thisElement).nodes;
    %
    %     nPoint = length(thisNodeArray);
    %
    %     for iPoint = 1:nPoint
    %
    %         thisIndex = thisNodeArray(iPoint);
    %         thisPoint = mesh.POS(thisIndex, :);
    %         pointMatrix(iPoint, :) = thisPoint;
    %
    %     end
    %
    %     % Calculate the volume of an irregular tertrahedral element
    %
    %     p1 = pointMatrix(1, :);
    %     p2 = pointMatrix(2, :);
    %     p3 = pointMatrix(3, :);
    %     p4 = pointMatrix(4, :);
    %
    %     v1 = p1 - p4;
    %     v2 = p2 - p4;
    %     v3 = p3 - p4;
    %
    %     % with function
    %     testVolume = tetrahedralVolume(v1, v2, v3);
    %
    %     % So Vcell are the volumes of the different irregular tetrahedron in the
    %     % mesh

    %%

    % Duration of the simulation [s]
    nSamples = ceil(recordingTime/dt);
    % Interrupted source, time before switching of [s]
    %     sourceon_steps = ceil(sourceOnTime/dt);

    %     s1 = w1 .* ones(1, sourceon_steps);
    %     source = [s1, zeros(1, nSteps-sourceon_steps)]; %for source on time, interrupted source
    %
    %     s = zeros(nElementAll, 1);
    %     % Interpolation of the source point to the surronding elements.
    %     s(sel_list) = li * source(1) / lsum;

    %% Diffusion parameters %%

    % Total  mesh volume [m^3]
    V = sum(Vcell);

    if strcmp(filename, 'sphere') || strcmp(filename, 'NewTransfSphere') || strcmp(filename, 'sphereTransfinite')
        Vr = 4 / 3 * pi * (dLx / 2)^3;
        Sr = 4 * pi * (dLx / 2)^2;
        mfp = 4 * Vr / Sr; % Mean-free-path for a sphere
    elseif strcmp(filename, 'cube') || strcmp(filename, 'Cube_attractors')
        mfp = 4 * V / S; % Mean-free-path for a sphere
    else
        mfp = 4 * V / S;
    end

    % Mean-free-time
    mft = mfp / c0;
    % Diffusion coefficient
    D = mfp * c0 / 3;
    % Atmospheric attenuation [m^-1]
    m_atm = 1e-5;

    %% Time loop initialisation %%
    w_new = zeros(nElementAll, 1);
    w = w_new;
    w_old = w;
    w_old(sel_list) = w1;
    
    gamma_zero = dt * ((D .* fel) + fb) ./ Vcell;

    w_rec = zeros(nReceivers, nSamples);


    for iSample = 1:nSamples
      
        w_new = (w_old .* (1 - gamma_zero) - c0 * m_atm * 2 * dt * w + ...
            2 * dt * D * (Fmat * w) ./ Vcell) ./ (1 + gamma_zero);
        w_old = w;
        w = w_new;

        w_rec(1:nReceivers, iSample) = w_new(Rpel);

    end


    for iReceiver = 1:nReceivers

        figure(3)
        tv = (0:1:nSamples-1)*dt;

%         plot(tv, 20*log10(sqrt(abs(w_rec(iReceiver, :))*rho)*c0/(2e-5)))
        plot(tv, 20*log10(sqrt(abs(w_rec(iReceiver, :))*rho)*c0))
%         plot(tv,10*log10(abs(w_rec(iReceiver, :))))

        hold on

    end

    xlim([0.02, 1])
    ylim([-100, 0])
    %     axis([0, recording_time, 0, max(10*log10(abs(w_rec(j, :)*rho*c0^2/(2e-5)^2))) + 10]);
    title('SPL evolution at the receivers', 'fontsize', 9, 'fontsize', 9)
    xlabel('time in s', 'fontsize', 9)
    ylabel('Sound pressure level [dB] at receiver elements', 'fontsize', 9)
    %     hold off


    envelopes{iBand, 1} = w_rec;
    envelopes{iBand, 2} = dt;


end


%% Post processing
%
nSec = str2num(answer{6, 1});
% nSec=nSec+0.01;
thisDt = envelopes{1, 2};
thisFs = round(1/thisDt);


r = norm(abs(sourcePosition-receiverPositions));

t0 = r / c0;

fs = round(1/dt);
t0Samples = round(fs*t0);
% pDirect = 1 / (4 * pi * r);


lowerLimit = centerFrequencies(1) * 2^((-1 / 2) * (1 / nthOctave));
upperLimit = centerFrequencies(end) * 2^((1 / 2) * (1 / nthOctave));
% totalBandWidth = upperLimit - lowerLimit;


maxBandWidth = centerFrequencies(end) * (2^(1 / nthOctave) - 1) / (2^(1 / (2 * nthOctave)));

% Normalize envelopes to bandwidth per max bandwidth, and cut-off samples
% before arrival time t0

for iBand = 1:nBands

    for iReceiver = 1:nReceivers

        w_rec = envelopes{iBand, 1}(iReceiver, :);
        thisEnvelopeP = sqrt(abs(w_rec*rho)) * c0 / (2e-5);

        thisFc = centerFrequencies(iBand);
        thisBandWidth = thisFc * (2^(1 / nthOctave) - 1) / (2^(1 / (2 * nthOctave)));

        thisEnvelopeP = thisEnvelopeP * thisBandWidth / maxBandWidth;
        envelopesNew{iBand, 1}(iReceiver, :) = thisEnvelopeP(t0Samples:end);


    end

end





% thisBand = 1;
% 
% thisEnvelopeOld = envelopes{thisBand,1};
% thisEnvelopeNew = envelopesNew{thisBand, 1};
% 
% nSamplesOld = length(thisEnvelopeOld);
% nSamplesNew = length(thisEnvelopeNew);
% 
% tvNew = (0:1:nSamplesNew - 1) * thisDt;
% tvOld = (0:1:nSamplesOld - 1) * thisDt;
% 
% 
% figure(5)
% plot(tvNew, 20*log10(thisEnvelopeNew), 'lineWidth', 1, 'Color', 'b')
% hold on
% plot(tvOld, 20*log10(thisEnvelopeOld), 'lineWidth', 1, 'Color', 'r')
% 
% 
% 
% xlabel('time (s)')
% ylabel('Amplitude (dB)')
% 
% 
% aspectRatio = [10, 3, 1];
% pictureWidth = 7.5;
% a = 1;
% set(gcf, 'units', 'centimeters', 'position', [0, 0, pictureWidth, (pictureWidth / aspectRatio(1) * aspectRatio(2) * a) + 1])
% set(gca, 'fontSize', 9)
% 
% 
% xlim([0, 1])
% ylim([-60, -20])
%%

envelopes(:, 1) = envelopesNew;


prompt = 'Save envelopes ?';
saveYesOrNo = input(prompt, 's');

if strcmpi(saveYesOrNo, 'y')

    nameSave = sprintf('de_scenario_%0.0f_distance_%0.0f_mm_mesh_%0.0fm_%s', scenario, distanceMm, meshLength, boundaryCondition);
    save(nameSave, 'envelopes')


end


%
% load('polygons_scenario_1_de.mat');
% drawpolygons(polygons, 3)
% axis off

% view(45, 45)
