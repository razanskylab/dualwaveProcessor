%% compute one-sided fft of time-domain signal x
function [P1, f] = get_1D_fft(x, fs, flagDisplay)

    % make sure x is column vector
    x = x(:);
    
    % remove DC offset
    x = x - mean(x);

    % make sure x has even number of elements
    if size(x,1)/2 ~= 0
        x = x(1:end-1);
    end
    L = size(x, 1);

    y = fft(x);
    % two-sided amplitude spectrum
    P2 = abs(y ./ L);
    % one-sided amplitude spectrum
    P1 = P2(1 : (L/2)+1);
    P1(2:end-1) = 2 .* P1(2:end-1);

    % frequencies
    f = fs .* ((0:(L/2))./L);
    f = f(:);

    if flagDisplay
        figure('Name', 'FFT plot');
        plot(f/1e6, P1);
        xlabel('f [MHz]');
        ylabel('Amplitude');
        title('Single-sided Amplitude Spectrum');
    end

end