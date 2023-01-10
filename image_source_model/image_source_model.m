clearvars
close all

folders=addpath('data','functions');

%% 01.0 Loading Data

% Surface impedance of bouding surfaces

% If you want to use discrete Impedance data
% load('Carpet.mat')
% load('Panel.mat')

% nSamplesCarpet = length(Carpet.Measurement_frequency_vector);
% fsCarpet  = (Carpet.Measurement_frequency_vector(2)-Carpet.Measurement_frequency_vector(1))*nSamplesCarpet;
% carpetAbsorptionMeasured  = Carpet.Measurement_spectra_Z(:,3);


% Geometric data

% Function to extract geometrical data from .obj file
swapYzyIsUp = true;
% objFileName = 'Scenario_2_simplified_IMS.obj';
objFileName = 'Scenario_1_simplified_IMS.obj';

[roomGeometryVertices, polygon] = convert_obj_file_for_ism_model(objFileName, swapYzyIsUp);

% Directivity data

% Receiver (Read the hrtf)
load('hrtf_aachen_2522_points_left_and_right_44k_phase_corrected_fd.mat')

% Take the complex conjugate due to matlab convention
hrtfLeftClfSelected = conj(hrtfLeftClfSelected(1:end, 1:48000));

recDirectivityPositions = cartCoord;
recDirectivityKnownValues = hrtfLeftClfSelected;

% recDirectivityKnownValues=abs(recDirectivityKnownValues);

% load('HATS_DIRECTIVITY_BNK.mat')
% recDirectivityPositions = cartCoord;
% recDirectivityKnownValues = HATSP4K;

% Sound source
load('b_k_dodecahedron_loudspeaker_directivity_one_active_driver.mat')

srcDirectivityPositions = sourceDirectivityMeasurementSphericalCoordinates;
srcDirectivityKnownValues = sourceDirectivityMeasurementPressure;

%% 02.0 General settings

% General
% Speed of sound in air (m/s)
c0 = 343;
% Density of air (Pa)
rho0 = 1.213;

% Frequency range and resolution:
% Lowest frequency
fStart = 0;
% Frequency bin size
df = 0.5;
% Highest frequency
fStop = 24000 - df;
% Generate frequency array
frequency = (fStart:df:fStop);

% Image Source Method
% Image source order
reflectionOrder = 2;

% Transducers
% Directivity
% Source directivity true or false
directivityIsOnSrc = false;
% Receiver directivity true or false
directivityIsOnRec = false;
% Directivity interpolation settings (nDataPoints: amount of known directivities closest to the angle of incidence used for interpolation, powerValue: )
powerValue = 1;
nDataPoints = 1;

% Boundary conditions

% Scattering on (true) or off (false)
scatteringIsOn = false;

%  Phase reversal upon wall interaction on(true) or off(false) (if using random
%  incidence absorption coefficients not if using complex surface impedance
%  with angle formula (then angle dependent complexwall impedance is used and phase reversal should not be used))
phaseReverseIsOn = false;

% Impedance
%  Specific impedance (On/Off)
specificImpedance = true;
% Turn on or off angle dependence
angleDependance = 1;

% Bandpass filter settings
% Blackman windowed sinc on or off
windowedSincIsOn = true;
% Blackman windowed sinc filter width (for ifft).
windowLength = 500;
% Blackman windowed sinc sample frequency
fsBlackmanWindow = 48000;

% Set lower and upper frequency limits for model circumscribed by the
% available resolution of absorption characteristic data
firstBandFmin = 125 / (2^(1 / 2));
lastBandFmax = 2000 * (2^(1 / 2));

fs = (length(frequency)) * 2 * df - 1;
% Characteristic impedance of air
z0 = c0 * rho0;

%% 03.0 Transducers

% General
% Receiver is placed at ear position with respect to the center of the head
receiverIsAtEarPosition = false;
leftOrRight = 'right';

%% 03.1 Source(s)

% Positions
% Rows sources columns x, y, z coordinates

% Source receiver distance 1.5 m
% sourcePositions = [3.04, 2.59, 1.62];
% Source receiver distance 3.5 m
sourcePositions = [1.36, 3.76, 1.62];
% sourcePositions(2,:) = [1,2,2];

% Source Directivity

% Using the same directivity file for all sources. comment out if you want to
% use different directivity per source, and assign each source directivity
% accordingly.
nSources = size(sourcePositions, 1);
sourceIsPresent = true(nSources, 1);
srcDirectivityKnownPositions(sourceIsPresent, 1) = {srcDirectivityPositions};


% Source rotation

% sourceRotationsDegrees gives the rotation around x,y and z axis per source (rotation needs to be defined for each source)
% Inititalize variable
sourceRotationsDegrees = zeros(nSources, 3);
sourceRotationsDegrees(1, :) = [170, 0, -30];
% Set the sequence of rotations [1,2,3] = first x then y then z-axis rotation.
srcRotationOrder = [1, 2, 3];
sourcesRotationOrder = repmat(srcRotationOrder, nSources, 1);

% Place variables in structure
RotatePointsSource.RotationDegrees = sourceRotationsDegrees;
RotatePointsSource.RotationOrder = sourcesRotationOrder;
RotatePointsSource.PointCoordinates = srcDirectivityKnownPositions;

% This is the rotated set of known directivity values of the source.
[srcDirectivityPositionsRotated] = rotate_points_around_axes(RotatePointsSource);


%% 03.2 Receiver

% Position(s)

% Receiver positions (rows= receivers columns = x, y, z position)
receiverPositions(1, :) = [4.26, 1.76, 1.62];

iReceiverPosition = 1;
thisReceiverPosition = receiverPositions(iReceiverPosition, :);

% Directivity

% Rotation
% Rotate receiver to face source (Rotation around vertical axis). 
thisPos = sourcePositions - thisReceiverPosition;
x = thisPos(1);
y = thisPos(2);
angle = -atan(y/x) * (180 / pi) - 180;

% Rotate receiver
angleRot = 270;
angle = angle + angleRot;

% The Head And Torso Simulator has two transducers, there is only one set of measured
% directivity responses for the left ear, the set for the right ear can be
% made from the left ear due to bi-lateral symmetry of the HATS.

flipAlong = 'x';
[recDirectivityPositionsHats, earPos] = create_hats_set(recDirectivityPositions, flipAlong);

nEars = size(recDirectivityPositionsHats, 1);

% Initialize variables (no rotation situation)
receiverRotationsDegrees = zeros(nEars, 3);

% Set x, y and z axis rotation
% E.g. receiverRotationsDegrees(i,:) = [0,60,120];
receiverRotationsDegrees(:, 3) = angle;

% Set order of rotation for x, y and z axis.
% E.g. receiverRotationOrder(i,:)=[1,3,2];
receiverRotationOrder = repmat([1, 2, 3], nEars, 1);

% Place variables in structure
RotatePointsReceivers.RotationDegrees = receiverRotationsDegrees;
RotatePointsReceivers.RotationOrder = receiverRotationOrder;
RotatePointsReceivers.PointCoordinates = recDirectivityPositionsHats;

% Rotate the directivity measurement coordinates around the origin
[recDirectivityPositionsRotated] = rotate_points_around_axes(RotatePointsReceivers);

% Take the rotated left or right ear directivity measurement coordinates
if strcmpi(leftOrRight, 'left')
    [thisRecDirectivityPositionsRotated] = recDirectivityPositionsRotated{1, 1};
elseif strcmpi(leftOrRight, 'right')
    [thisRecDirectivityPositionsRotated] = recDirectivityPositionsRotated{2, 1};
end

% Place variables in structure
RotatePointsEars.RotationDegrees = receiverRotationsDegrees;
RotatePointsEars.RotationOrder = receiverRotationOrder;
RotatePointsEars.PointCoordinates = earPos;

% Rotate the directivity measurement coordinates around the origin
[earPosRot] = rotate_points_around_axes(RotatePointsEars);

if receiverIsAtEarPosition == true

    if strcmpi(leftOrRight, 'left')
        thisReceiverPosition = earPosRot{1, 1} + thisReceiverPosition;
    elseif strcmpi(leftOrRight, 'right')
        thisReceiverPosition = earPosRot{2, 1} + thisReceiverPosition;
    end

end

%% 04.0 Geometry (pre-processing)
% % 4.1 Bounding Planes

% Function to make a unique description of all planes from a convex hull
[walls, wallSurfaceAreas, volumeRoom] = create_unique_planes_description(roomGeometryVertices);

S = sum(wallSurfaceAreas);
% The mean free path length
r = 4 * volumeRoom / S;
nPlanes = size(walls, 1);

%% 05.0 Material properties

% Polygons in plane 
% Polygons are defined by the vertices that make up the polygon
% Which bounding surfaces have a polygon on it (e.g. a partially covered wall (surface) by a carpet (sub-surface))


% Convert polygon coordinates (3D) to local coordinates (2D)
[polygon2d] = convert_3d_to_local_coordinates(polygon);
% Calculate the areas of the polygons.
[polygonAreas] = calculate_polygon_areas(polygon2d);


% The vertices of the polygons can be defined in any order. It is important
% to have them in clockwise or counterclockwise order for calculation
% purposes the sortpolygonvertices calculates the possible clockwise or
% counterclockwise ordered pairs of vertices.

% [polygon1Sorted]=sortpolygonvertices(polygon{1,1});

% polygon{1,1}=polygon1Sorted;


% Find which polygons lie in which plane.
[inPlane] = find_polygon_in_plane(polygon, walls);

uniquePlanes = unique(inPlane);
isInPlane = unique(inPlane) ~= 0;
planesWithPolygon = uniquePlanes(isInPlane);


%% 05.1 Surface impedance

% Generate Impedance using the miki model (frequency dependent)
% parametersMiki=[air flow resistivity (Pa.s/m),thickness of the material (m)]

parametersMiki(1, :) = [63000, 37 * 10^-3];
[ZsPanel] = miki_model_impedance(parametersMiki, frequency, c0, specificImpedance, z0);
parametersMiki(1, :) = [6.5 * 10^3, 8.5 * 10^-3];
[ZsCarpetPor] = miki_model_impedance(parametersMiki, frequency, c0, specificImpedance, z0);

% Generate impedance using the Cox method as described by Baltazar
% Thickness of backing air layer (m), Mass per unit area of the membrane (kg/m^2), and
% dampening factor (Pa.s/m).

parametersCox = [2.6 * 10^-3, 2.3, 5 * 10^3];
[zsCarpetRes] = ideal_membrane_resonator_impedance(parametersCox, frequency, z0, c0, specificImpedance);
ZsCarpet = (zsCarpetRes .* ZsCarpetPor) ./ (zsCarpetRes + ZsCarpetPor);

% Formula for normal incidence absorption coefficient to impedance
% aNorm = 0.007984;
aNorm = 0.2;

% Impedance for hardwall boundary
if specificImpedance == true
    Zval = -(sqrt(1-aNorm) + 1) / (sqrt(1-aNorm) - 1);
elseif specificImpedance ~= true
    Zval = z0 * -(sqrt(1-aNorm) + 1) / (sqrt(1-aNorm) - 1);
end

% ZvalPolygon=19;
ZsHard = ones(1, length(frequency)) * Zval;


% ZsPolygon=ones(1,length(frequencies))*ZvalPolygon;

% Assign impedance to wall surface
ZsWalls = repmat(ZsHard, size(walls, 1), 1);
nFrequencies = length(frequency);
nPolygons = size(polygon, 1);
ZsPolygons = repmat(ZsPanel, size(polygon, 1), 1);


% Assign the right impedance to the polygons

% % Scenario I
ZsPolygons(1, :) = ZsHard;
ZsPolygons(2, :) = ZsPanel;
ZsPolygons(3, :) = ZsPanel;
ZsPolygons(4, :) = ZsHard;
ZsPolygons(5, :) = ZsPanel;
ZsPolygons(6, :) = ZsPanel;
ZsPolygons(7, :) = ZsPanel;
ZsPolygons(8, :) = ZsPanel;
ZsPolygons(9, :) = ZsPanel;
ZsPolygons(10, :) = ZsCarpet;

% % Scenario II
% ZsPolygons(1, :) = ZsHard;
% ZsPolygons(2, :) = ZsHard;
% ZsPolygons(3, :) = ZsPanel;
% ZsPolygons(4, :) = ZsPanel;
% ZsPolygons(5, :) = ZsPanel;
% ZsPolygons(6, :) = ZsCarpet;

ZsWallsAndPolygons = ZsWalls;
% Impedance values for both the walls and the polygons
ZsWallsAndPolygons(end+2:end+1+nPolygons, 1:nFrequencies) = ZsPolygons;


% Calculate bare wall area by deducting the polygon areas from the full
% wall surfaces

nPolygon = size(polygon, 1);

for iPolygon = 1:1:nPolygon

    thisPlane = inPlane(iPolygon);
    thisWallArea = wallSurfaceAreas(thisPlane);
    thisPolygonArea = polygonAreas(iPolygon);
    wallSurfaceAreas(thisPlane) = thisWallArea - thisPolygonArea;

end

%% 06.0 Surface scattering
% (used 62.5 Hz to 4000 Hz band because then the desired range (125 Hz to 2000 Hz band) can be interpolated)

% % scenario I scattering coefficients
scatteringCoefficients = [0.12, 0.15, 0.25, 0.50, 0.53, 0.55, 0.60];

% % scenario II scattering coefficients
% scatteringCoefficients = [0.1, 0.12, 0.15, 0.25, 0.33, 0.35, 0.4];

% Giving the same scattering coefficients to each wall, can also be set per
% wall rows are walls columns are used to keep coefficients per frequency
% band.
scatteringCoefficientsWalls = repmat(scatteringCoefficients, size(walls, 1), 1);
scatteringCoefficientsPolygons = repmat(scatteringCoefficients, size(polygon, 1), 1);

% Place the scattering coefficients of the walls floor and ceiling together
% with the scattering coefficients of the polygons on the walls floor and
% ceiling.
scatteringWallsAndPolygons = scatteringCoefficientsWalls;

nBands = size(scatteringCoefficientsWalls, 2);

scatteringWallsAndPolygons(end+2:end+1+nPolygons, 1:nBands) = scatteringCoefficientsPolygons;

% Increase the resolution of the scattering coefficients
nSurfaces = size(scatteringWallsAndPolygons, 1);
scatteringWallsAndPolygonsHighRes = nan(nSurfaces, nFrequencies);


% Generate Octave bands (Base-2)
fcLow = 62.5;
fcHigh = 4000;
nthOctave = 1;

[Bands] = generate_octave_bands(fcLow, fcHigh, nthOctave);
centerFrequenciesBands = Bands.frequenciesCenter;
lowerLimitsBands = Bands.frequenciesLower;
upperLimitsBands = Bands.frequenciesUpper;


% Interpolate scattering coefficients to fine resolution
for iSurface = 1:nSurfaces

    theseScatteringCoefficients = scatteringWallsAndPolygons(iSurface, :);
    scatteringWallsAndPolygonsHighRes(iSurface, :) = interp1(centerFrequenciesBands, theseScatteringCoefficients, frequency, 'pchip', NaN);

end


[aUniWallsAndPolygons] = calculate_random_incidence_absorption(ZsWallsAndPolygons);

nBands = length(centerFrequenciesBands);
nSurfaces = size(aUniWallsAndPolygons, 1);
absorptionAverageCoefficients = nan(nSurfaces, nBands);


for iBand = 1:nBands

    for jSurface = 1:nSurfaces

        theseCoefficients = aUniWallsAndPolygons(jSurface, :);
        thisLowerLimit = lowerLimitsBands(iBand);
        thisUpperlimit = upperLimitsBands(iBand);

        index = frequency >= thisLowerLimit & frequency < thisUpperlimit;

        thisBandAbsorptionCoefficients = theseCoefficients(index);
        thisAverageCoefficient = mean(thisBandAbsorptionCoefficients);
        absorptionAverageCoefficients(jSurface, iBand) = thisAverageCoefficient;


    end

end

areaWallsandPolygons = [wallSurfaceAreas; 0; polygonAreas];
[absorptionCoefficientBandAverage] = calculate_average_coefficients(areaWallsandPolygons, absorptionAverageCoefficients);
[scatteringCoefficientBandAverage] = calculate_average_coefficients(areaWallsandPolygons, scatteringWallsAndPolygons);

%% 07.0 Hybrid model (transition order)

diffuseFieldPercentage = 90;
orderPerBand = (log(1-(diffuseFieldPercentage / 100)) ./ log(1-scatteringCoefficientBandAverage));

%% 08.0 Calculation

spl = @(input)(20 * log10(abs(input)));


% Function that generates the image sources (Open function for more comments on its functioning)
[~, Ism2, pointsOfReflection, anglesOfIncidenceWall, reflectionSequences] = determine_image_sources(walls, sourcePositions, thisReceiverPosition, reflectionOrder);

% Mark which reflections are on surfaces with a polygon in them.
[intersectionHasPolygon, polygonInPlane] = mark_planes_with_polygon(planesWithPolygon, reflectionSequences);

% Check if Polygon is hit
[reflectionSequences] = is_reflection_inside_polygon(reflectionSequences, polygonInPlane, pointsOfReflection, intersectionHasPolygon, polygon, nPlanes, inPlane);

% Cell array with all effective image sources per reflection level
effectiveImageSourcesPositions = Ism2.ispclean;
effectiveImageSourcesPlaneNo = Ism2.isnclean;
% Cell array with the sequence of reflection for each image sources
effectiveImageSourcesNames = Ism2.isnmsclean;

% 5.5 DETERMINE ANGLES OF INCIDENCE AT TRANSDUCER POSITIONS
% Function to determine the angles of incidence with the sources and the
% receivers
[angleOfIncidenceSource, angleOfIncidenceReceiver] = calculate_angles_of_incidence(sourcePositions, pointsOfReflection, thisReceiverPosition);

% 5.6 INTERPOLATE ABSORPTION DATA TO MATCH MODELS RESOLUTION
% Function that increases the absorptive characteristics of bounding surfaces resolution to the models resolution

% [interpolatedReflectionCoefficients] = interpolatecoefficients(centerFrequenciesBands, reflectionCoefficients, frequencies);
% [vq2]=FUN_MIN_IAR(fvZval,refcfts,freqarr);


% 5.7 CREATE MATRIX WITH TOTAL REFLECTION COEFFICIENT PER IMAGE SOURCE
% Function to make a matrix of the total absorption per image source

[totalReflectionPerImageSource] = calculate_total_reflection_per_image_source(reflectionSequences, ...
    ZsWallsAndPolygons, anglesOfIncidenceWall, specificImpedance, angleDependance, z0);


[totalScatteringPerImageSource] = calculate_total_scattering_per_image_source(reflectionSequences, ...
    scatteringWallsAndPolygonsHighRes);

% 5.8 DETERMINE (IMAGE) SOURCE RECEIVER DISTANCES
% Function to determine all (image)source receiver distances
[imageSourceReceiverDistances] = calculate_source_receiver_distance(effectiveImageSourcesPositions, sourcePositions, thisReceiverPosition);


% Place receiver directivity in cell array
positionsMeasuredReceiverDirectivity{1, 1} = thisRecDirectivityPositionsRotated;

% 5.9 DETERMINE ORTHODROMIC DISTANCES (DIRECTIVITY)
%Determine orthodromic distances for the source between the angles of incidence at the
%source and the knownpoints of the measured directivity.
[orthodromicDistancesSources] = calculate_orthodromic_distances(angleOfIncidenceSource, srcDirectivityPositionsRotated);

%Determine orthodromic distances for receiver
[orthodromicDistancesReceiver] = calculate_orthodromic_distances(angleOfIncidenceReceiver, positionsMeasuredReceiverDirectivity);

%extract the right range of frequency bands from the B&K Hats data. 33
%means the 33rd 1/12th octave band in this case since this corresponds to
%the lower limit of the 250 Hz octave band 250/(2^(1/2))
receiverHatsMeasurement4kExtracted = recDirectivityKnownValues(1:102, 33:end);

% 5.10 INTERPOLATE TRANSDUCER DIRECTIVITY DATA USING ORTHODROMIC DISTANCES
%Interpolate source directivity data
[interpolatedDirectivitySource] = interpolate_unknown_directivities(orthodromicDistancesSources, srcDirectivityKnownValues, powerValue, nDataPoints);

%Interpolate receiver directivity data
[interpolatedDirectivityReceiver] = interpolate_unknown_directivities(orthodromicDistancesReceiver, hrtfLeftClfSelected, powerValue, nDataPoints);

% 5.11 CREATE Q VALUES PER (IMAGE) SOURCE
%Make interpolated pressures into Q values for source
[interpolatedQSource] = create_q_values(interpolatedDirectivitySource, sourceDirectivityMeasurementPressure);

%Make interpolated pressures into Q values for receiver

% [interpolatedQReceiver] = createqvalues(interpolatedDirectivityReceiver, receiverHatsMeasurement4kExtracted);

interpolatedQReceiver = interpolatedDirectivityReceiver;

% 5.12 INCREASE TRANSDUCER DIRECTIVITY RESOLUTION TO MODELS RESOLUTION
%Increase frequency resolution to models frequency resolution for source
[interpolatedQHighResolutionSource] = increase_resolution(interpolatedQSource, fvds, frequency);

%Increase frequency resolution to models frequency resolution for receiver
% [interpolatedQHighResolutionReceiver] = increaseResolution(interpolatedQReceiver, fvdr, frequency);
interpolatedQHighResolutionReceiver = interpolatedQReceiver;

% 5.13 CREATE BLACKMAN WINDOW


%Function to make a blackman window
[blackmanWindowedSinc] = generate_blackman_windowed_sinc(firstBandFmin, lastBandFmax, fsBlackmanWindow, windowLength);

% 5.14 CALCULATION OF THE ROOM TRANSFER FUNCTION
%The frequency domain construction of the impulse response. the last three
%digits of the function inputs are: Source directivity on(1)/off(0); Receiver directivity on(1)/off(0); Pi degree phase shift upon wall interaction on(1)/off(0)
% [p2,sumtot,sumdirectfield]=FUN_KEY_FDC(distmat,freqarr,sos,absarr,intQhr,intQhrr,dirsrc,dirrec,prev);

[impulseResponsePerImageSource, totalTf, directFieldTf] = ...
    generate_transfer_function(imageSourceReceiverDistances, ...
    frequency, c0, totalReflectionPerImageSource, interpolatedQHighResolutionSource, ...
    interpolatedQHighResolutionReceiver, directivityIsOnSrc, directivityIsOnRec, ...
    phaseReverseIsOn, totalScatteringPerImageSource, scatteringIsOn);


totalTf(isnan(totalTf)) = 0;
directFieldTf(isnan(directFieldTf)) = 0;


totalTf = totalTf';
directFieldTf = directFieldTf';
% 5.15 TIME DOMAIN IMPULSE RESPONSE

%Function to create the impulse response (Inverse Fourier transform)
nZeros = 0;
zeroArray = complex(zeros(1, nZeros));
totalTfSymmetric = [totalTf(1:end), zeroArray, conj(fliplr(totalTf(2:end)))];
directFieldTfSymmetric = [directFieldTf(1:end), zeroArray, conj(fliplr(directFieldTf(2:end)))];


% Conjugate symmetric check
% isequal(tf, conj(ft([1, end:-1:2])))
% Unwindowed impulse response
totalIr = ifft(totalTfSymmetric);
directFieldIr = ifft(directFieldTfSymmetric);

% Filter with Blackman windowed sinc filter
scaleFactor = (frequency(end) - frequency(1)) / (lastBandFmax - firstBandFmin);
% scaleFactor = 1;


if windowedSincIsOn == true

    totalIrWindowed = cconv(totalIr, blackmanWindowedSinc) * scaleFactor;
    directFieldIrWindowed = cconv(directFieldIr, blackmanWindowedSinc) * scaleFactor;

    shift = -round(windowLength*(1 / 2));
    totalIrWindowed = circshift(totalIrWindowed, shift);
    totalTfWindowed = fft(totalIrWindowed) / fs;

    directFieldIrWindowed = circshift(directFieldIrWindowed, shift);
    directFieldTfWindowed = fft(directFieldIrWindowed) / fs;
    % remove artifacts from end of impulse response
    tdrbv = totalIrWindowed;
    tdrbv(end-1000:end) = [];

elseif windowedSincIsOn ~= true

    tdrbv = totalIr;


end


% use a double cosine window
Nleft = 50;
Nright = 1000;
plt = 0;

[tdrbv, ~] = time_windowing(tdrbv, Nleft, Nright, plt);

%% 09.0 Visualisation

pictureWidth = 15;
aspectRatio = [5, 3, 1];
fontSize = 9;
lineWeight = 0.5;

nSubPlots = 1;
nthOctave = 1;


% Plot of scattering coefficients
figure(1)
plot(frequency, scatteringWallsAndPolygonsHighRes(1, :))

title('Averaged scattering coefficient')
xlabel('Frequency (Hz)')
ylabel('Scattering coefficient (-)')
set_figure(pictureWidth, aspectRatio, nSubPlots, nthOctave)

% Plot of the windowed sinc filter
figure(2)

logXScale = false;
plot(blackmanWindowedSinc)

title('Blackman windowed sinc filter')
xlabel('Sample no. (-)')
ylabel('Gain (normalized to unity) (-)')

nSamplesWindow = length(blackmanWindowedSinc);
yMax = max(abs(blackmanWindowedSinc));

xlim([1, nSamplesWindow])
ylim([-yMax, yMax])
lay_out_figure(fontSize, lineWeight, pictureWidth, aspectRatio, nSubPlots, logXScale);


% Plot of the transfer function
nSamples = length(totalTfSymmetric);
df = fs / nSamples;
fv = (0:1:nSamples - 1) * df;

if directivityIsOnRec == true

    if strcmpi(leftOrRight, 'right')
        color = 'r';
    elseif strcmpi(leftOrRight, 'left')
        color = 'k';
    end

else

    color = 'k';

end


figure(3)


logXScale = false;

plot(fv, spl(totalTfSymmetric), color)
hold on
plot(fv, spl(directFieldTfSymmetric), 'b')

nameLegend = sprintf('IR (ISM order:%0.0f)', reflectionOrder);
lgnd = legend(nameLegend, 'IR (direct field)', 'location', 'south east');
lgnd.ItemTokenSize = [10, 5];
set(lgnd, 'Box', 'off')
set(lgnd, 'color', 'none');

xlim([0, 2000*2^(1/2)])
ylim([-50, 0])

xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)')

lay_out_figure(fontSize, lineWeight, pictureWidth, aspectRatio, nSubPlots, logXScale);


% Plot of the impulse response
nSamples = length(tdrbv);
tv = (0:nSamples - 1) * (1 / fs);

figure(4)
logXScale = false;


plot(tv, spl(tdrbv))
% xlim([0, 0.04])
% ylim([-0.06 0.06])
lay_out_figure(fontSize, lineWeight, pictureWidth, aspectRatio, nSubPlots, logXScale);


figure(41)
plot(tv, tdrbv)
xlim([0, 0.04])
xlabel('time (s)')
ylabel('Pressure (Pa)')
lay_out_figure(fontSize, lineWeight, pictureWidth, aspectRatio, nSubPlots, logXScale);


% Plot Geometry
figNum = 5;
plotFrequency = 2000;
fontSize = 9;

plot_polyhedron_boundary(roomGeometryVertices, figNum)
hold on
scatter3(sourcePositions(1), sourcePositions(2), sourcePositions(3), 'ok', 'filled')
scatter3(thisReceiverPosition(1), thisReceiverPosition(2), thisReceiverPosition(3), 'or', 'filled')


leftEarPos = earPosRot{1, 1} + thisReceiverPosition;
rightEarPos = earPosRot{2, 1} + thisReceiverPosition;
scatter3(leftEarPos(1), leftEarPos(2), leftEarPos(3), '+k')
scatter3(rightEarPos(1), rightEarPos(2), rightEarPos(3), '+r')

draw_polygons(polygon, figNum)
plot_reflections(pointsOfReflection, sourcePositions, thisReceiverPosition, reflectionSequences, nPlanes, nPolygons)
plot_angles_of_incidence(angleOfIncidenceReceiver, thisRecDirectivityPositionsRotated, thisReceiverPosition)
view(45, 35)
axis off

% Loop to check the numbers of the walls
%     nWall = size(walls, 1);
%
%     for iWall = 1:nWall
%         thisWall = walls{iWall, 1};
%         vn = sum(thisWall, 1) / 3;
%
%         str = sprintf('w:%d', iWall);
%         text(vn(1), vn(2), vn(3), str)
%         hold on
%     end

a = 1;
set(gcf, 'units', 'centimeters', 'position', [0, 0, pictureWidth, (pictureWidth / aspectRatio(1) * aspectRatio(2) * a) + 1])

surfNo = 12;
aUni = aUniWallsAndPolygons(surfNo,:);
aUniAvg = absorptionAverageCoefficients(surfNo,:);

figure(6)
plot(frequency, aUni, 'r')
hold on
scatter(centerFrequenciesBands, aUniAvg, 'sb', 'filled')


xlim([125 * 2^(-.5), 2000 * 2^(.5)])
ylim([0, 1])

xlabel('Frequency [Hz]');
ylabel('\alpha [-]');
lgnd = legend('$\alpha_{uni}$', '$\bar{\alpha}_{uni}$', 'location', 'northwest', 'interpreter', 'latex');


lgnd.ItemTokenSize = [10, 5];
set(lgnd, 'Box', 'off')
set(lgnd, 'color', 'none');

set_figure(pictureWidth, aspectRatio, nSubPlots, nthOctave)
set(gcf, 'Color', 'w')



%% 10.0 Saving data for the hybrid model

sourceReceiverDistance = imageSourceReceiverDistances{1, end};
arrivalTime = sourceReceiverDistance / c0;

prompt = 'Save (BR)IR ?';
saveYesOrNo = input(prompt, 's');

if strcmpi(saveYesOrNo, 'y')

    if directivityIsOnRec == false
        nameSave = sprintf('ism_ir_order_%d_no_polygons_%d_src_rec_dis%d.mat', reflectionOrder, nPolygons, round(sourceReceiverDistance, 1)*1000);
        save(nameSave, 'tdrbv', 'scatteringCoefficientBandAverage', 'absorptionCoefficientBandAverage', 'arrivalTime','directFieldTfSymmetric','totalTfSymmetric')
    elseif directivityIsOnRec == true

        xRot = receiverRotationsDegrees(1, 1);
        yRot = receiverRotationsDegrees(1, 2);
        zRot = receiverRotationsDegrees(1, 3);

        xOrd = receiverRotationOrder(1, 1);
        yOrd = receiverRotationOrder(1, 2);
        zOrd = receiverRotationOrder(1, 3);

        nameSave = sprintf('ism_brir_%s_order_%d_no_polygons_%d_src_rec_dis%d_rotation_%0.0f.mat', leftOrRight, reflectionOrder, nPolygons, round(sourceReceiverDistance, 1)*1000, angleRot);
        save(nameSave, 'tdrbv', 'scatteringCoefficientBandAverage', 'absorptionCoefficientBandAverage', 'arrivalTime')
    end

end
