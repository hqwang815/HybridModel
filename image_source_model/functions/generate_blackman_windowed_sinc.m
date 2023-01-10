function [B, sit, Bmw, A, x, x2] = generate_blackman_windowed_sinc(lowerLimit, upperLimit, sampleRate, M)


% lowerLimit, upperLimit, fsBlackmanWindow, filterWidth

%M should be even
x = (-.5 * M:1:.5 * M);
x2 = (0:1:M);

%Place the lower and upper limits in fc vector
fc(1) = lowerLimit / sampleRate;
fc(2) = upperLimit / sampleRate;

for i = 1:1:2
    %Calculate the tappered sync function
    sit = sin(2*pi*fc(i).*x) ./ (x);
    %because of discretisation there will be a NAN value at the symmetry point
    %in the filter so find this value.
    del = isnan(sit);
    %Change this value to the value it should approach
    sit(del) = 2 * pi * fc(i);
    %Create a Blackman Window
    Bmw = (0.42 - 0.5 * cos(2*pi*x2/M)) + 0.08 * cos(4*pi*x2/M);
    %Multiply Blackman Window with Sync function.
    A(i, :) = sit .* Bmw;
    %Normalize to unity gain
    A(i, :) = A(i, :) / sum(A(i, :));

    if i == 2
        %Perform spectral inverse on the upper limit window
        A(i, :) = A(i, :) * -1;
        %Add 1 to the central point
        A(i, del) = A(i, del) + 1;
        %Add low pass to high pass
        B = sum(A, 1);
        %Use spectral inverse to change band reject filter to bandpass filter.
        B = B * -1;
        %Add 1 to the central point
        B(1, del) = B(1, del) + 1;
    end

end

end