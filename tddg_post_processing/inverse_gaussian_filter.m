clear all
close all


scenario = 2;
distance = 3500;
frequencyRange = 'hf';


if strcmpi(frequencyRange, 'lf')
    nameTddg = sprintf('simulation_take2_scenario_%0.0f_lf_%0.0fmm', scenario, distance);


    % Limits where the simulation is valid
    lowerFrequency = 31.25 * 2^(-1 / 2);
    upperFrequency = 800;

    % Windowed sinc filter before IFFT
    lowerLimitSinc = 62.5 * 2^(-1 / 2);
    upperLimitSinc = 500 * 2^(1 / 2);

    saveName = sprintf('irs_spectrum_corrected_scenario_%0.0f_distance_%0.0f_mm_lf', scenario, distance);


elseif strcmpi(frequencyRange, 'hf')
    nameTddg = sprintf('simulation_take2_scenario_%0.0f_hf_%0.0fmm', scenario, distance);

    % Limits where the simulation is valid
    lowerFrequency = 100;
    upperFrequency = 2100;

    % Windowed sinc filter before IFFT
    lowerLimitSinc = 200;
    upperLimitSinc = 2000;

    saveName = sprintf('irs_spectrum_corrected_scenario_%0.0f_distance_%0.0f_mm_hf', scenario, distance);

end

IrsTddg = load(nameTddg);

%% Settings

% halfWidth = 0.17;
c0 = 343;
S0 = 1;
% r = 1.4756;


% Plotting
checkPos = 1;

% Tranfer functions
pictureWidth = 15;
aspectRatio = [5, 3, 1];
nthOctave = 1;

fontSize = 9;
lineWeight = 1;
nSubPlots = 1;
logXScale = 'log';

xMin = 62.5 * 2^(-1 / 2);
xMax = 2000 * 2^(1 / 2);
yMin = -100;
yMax = 0;

%% Anonymous functions
spl = @(input)(20 * log10(abs(input)));

%%

irs = IrsTddg.p_audio;
nPositions = size(irs, 1);


dtTddg = IrsTddg.dt_audio;
fsTddg = round(1/dtTddg);

halfWidth = IrsTddg.halfwidth;

posSrc(1) = IrsTddg.xs;
posSrc(2) = IrsTddg.ys;
posSrc(3) = IrsTddg.zs;


for iPosition = 1:nPositions


    irTddg = double(irs(iPosition, :));


    % Nleft = round(0.001 * fsTddg);
    % Nright = round(0.1 * fsTddg);
    % plt = 0;
    %
    % [irTddg, ~] = time_windowing(irTddg, Nleft, Nright, plt);


    posRec(1) = IrsTddg.rec_x(iPosition);
    posRec(2) = IrsTddg.rec_y(iPosition);
    posRec(3) = IrsTddg.rec_z(iPosition);

    r = norm(posSrc-posRec);

    pDirect = 1 / (4 * pi * r);
    splDirect = spl(pDirect);

    nSamplesTddg = length(irTddg);

    % Length of FFT use a power of two to increase fft efficiency
    nFft = 2^nextpow2(nSamplesTddg);

    % Take FFT
    tfTddg = fft(irTddg, nFft);

    % Length of one sided-spectrum
    nSamplesOneSided = ceil((nFft + 1)/2);


    % Calculate freefield Gaussian pulse response

    FWHM = 2 * halfWidth;
    alpha = 4 * log(2) / FWHM^2;

    tv = (0:1:nFft - 1) * dtTddg;
    pr = S0 * (r - c0 * tv) / (2 * r) .* exp(-alpha*(r - c0 * tv).^2);


    tfTddgFreefield = fft(pr);


    % Take one-sided spectra of transfer functions
    tfTddg = tfTddg(1:nSamplesOneSided);
    tfTddgFreefield = tfTddgFreefield(1:nSamplesOneSided);


    % Calculate the Greens function of the delta source
    df = fsTddg / nFft;
    fv = (0:1:nSamplesOneSided - 1) * df;
    omega = 2 * pi * fv;
    k0V = omega / c0;
    Pd = exp(-1i*k0V*r) / (4 * pi * r);

    % Windowed sinc filter
    M = nFft;
    sincFilter = generateBlackmanWindowedSinc(lowerLimitSinc, upperLimitSinc, fsTddg, M);
    tfSinc = fft(sincFilter);
    tfSinc = tfSinc(1:nSamplesOneSided);

    tfTddgNormalized = Pd ./ tfTddgFreefield .* tfTddg;


    % Use a set value below lower limit
    indexLow = find(fv >= lowerFrequency, 1);
    tfTddgNormalized(1:indexLow) = pDirect;

    % Use a set value above higher limit
    indexHigh = find(fv >= upperFrequency, 1);
    tfTddgNormalized(indexHigh:end) = pDirect;

    % Apply sinc filter
    tfNormFilt = tfTddgNormalized .* tfSinc;


    if iPosition == checkPos

        figure()
        plot(fv, spl(tfTddg))
        hold on
        plot(fv, spl(tfTddgFreefield))
        plot(fv, spl(tfTddgNormalized))
        plot(fv, spl(tfNormFilt))
        
        ylabel('SPL (dB)')
        xlabel('Frequency (Hz)')

        dSpl = 30;

        xlim([xMin, xMax])
        ylim([splDirect - dSpl, splDirect + dSpl])

        lay_out_figure(fontSize, lineWeight, pictureWidth, aspectRatio, nSubPlots, logXScale)

        lgnd = legend('TF (uncorrected)', 'TF (free-field Gaussian)', 'TF (corrected)', 'TF (corrected + bandpassed)');
        lgnd.ItemTokenSize = [10, 5];
        set(lgnd, 'Box', 'off')
        set(lgnd, 'color', 'none');

    end

    tfTddgNormalized = tfNormFilt(:);


    % Construct double sided spectrum

    if rem(nFft, 2)
        % Odd excludes Nyquist point
        tfTddgNormalizedTwoSided = [tfTddgNormalized; conj(tfTddgNormalized(end:-1:2))];
    else
        % Even includes Nyquist point
        tfTddgNormalizedTwoSided = [tfTddgNormalized; conj(tfTddgNormalized(end-1:-1:2))];
    end

    % Conjugate symmetric check
    isequal(tfTddgNormalizedTwoSided, conj(tfTddgNormalizedTwoSided([1, end:-1:2])))


    % Perform IFFT
    irFiltered = ifft(tfTddgNormalizedTwoSided, 'symmetric');
    irFiltered = circshift(irFiltered, round(M/2));


    Nleft = round(0.001*fsTddg);

    if strcmpi(frequencyRange, 'lf')
        Nright = round(1*fsTddg);
        plt = 0;
    elseif strcmpi(frequencyRange, 'hf')
        Nright = round(0.1*fsTddg);
        plt = 0;
    end


    [irFiltered, ~] = time_windowing(irFiltered', Nleft, Nright, plt);
    irsFiltered(iPosition, :) = irFiltered;


end


figure()

nSamples = length(irsFiltered(1, :));
tv = (0:1:nSamples - 1) * dtTddg;


plot(tv, irsFiltered(1, :))

saveOrNot = input('save file (y/n) ?', 's');


if strcmpi(saveOrNot, 'y')

    save(saveName, "irsFiltered", "fsTddg", "lowerLimitSinc", "upperLimitSinc")

elseif strcmpi(saveOrNot, 'n')

    disp('File not saved')
end



