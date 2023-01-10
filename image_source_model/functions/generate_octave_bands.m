
%% Syntax
% #  [Bands] = generate_octave_bands(fcLow, fcHigh, nthOctave)

%% Description
% Input arguments:
% # fcLow           lowest 1/1 octave band center frequency
% # fcHigh          highest 1/1 octave band center frequency
% # nthOctave   generate 1/nth octaves
% Output arguments:
% # Bands           structure with 1/nth base-2 band center frequencies, lower
% band limits and upper band limits 

%%  Examples
% Example

%%
function [Bands] = generate_octave_bands(fcLow, fcHigh, nthOctave)

% Create an nth octave center frequency array
fMin = (2^(1 / 2) * 2^(1 / (2 * nthOctave)) * fcLow) / 2;
fMax = (2^(1 / 2) * fcHigh) / 2^(1 / (2 * nthOctave));

a = nthOctave * log2(fMax/fMin);
frequenciesCenter = fMin * 2.^((1 / nthOctave) .* ([0:1:a]));
frequenciesLower = frequenciesCenter .* 2^-((1 / 2) * (1 / nthOctave));
frequenciesUpper = frequenciesCenter .* 2^((1 / 2) * (1 / nthOctave));

Bands.frequenciesCenter = frequenciesCenter;
Bands.frequenciesLower = frequenciesLower;
Bands.frequenciesUpper = frequenciesUpper;

end