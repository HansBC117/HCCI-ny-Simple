%function filtsignal = fftBPfilter(time,signal,bandpasses,bw,plotflag)

% time:       Column vector were each element contain a time [s] that
%             correspond to elemnts in the signal vector. Also if the
%			  signal is multidimensional array, the time has to a single
%			  column vector that applies to all the signal columns.
%
% signal:     Column vector with even number of elements containing the
%             signal to be processed. The signal may also be two or three
%			  dimensional arrays but then the filtering will be performed
%			  columnwise.
%
% bandpasses: Matrix with 4 columns and a number of rows corresponding to
%             the number bandpasses to be used. Each row defines a bandpass
%             by first and second element being the lower and upper band
%             frequencies [Hz]. Third and fourth element defines the
%             passfactors (between 0 and 1) at lower and upper band
%             frequencies. Passfactors between the band limits will be
%             interpolated linearly.
%
% bw:         Bandwidth [Hz] of the filter kernel as well as the bandwidth
%             of the transitions in the frequency respons. The transition 
%             bandwidth will be symetrical around the specified bandpass
%             limit.
%             
% plotflag:   Options: 'plot_on' or 'plot_off'. Determines whether or not
%             plots that show the filter performance will be generated.
%			  The plot option will only display the filterperformance of
%			  the signal in the first column if a multidimensional signal
%			  array is used.

function filtsignal = fftBPfilter(time, signal, bandpasses, bw, plotflag)

% Automatically truncate the signal to even length
L = size(signal, 1);
if mod(L, 2) ~= 0
    signal = signal(1:end-1);
    time   = time(1:end-1);
    L = L - 1;
end

fmax = (L / 2) / diff(time([1,end]));
M = 2 * round(2 / (bw / (2 * fmax)));
if M > L
    M = L;
end
f = [0:L/2-1]' / diff(time([1,end]));

KERNEL = zeros(L/2,1);
for i = 1:size(bandpasses,1)
    index = f >= bandpasses(i,1) & f < bandpasses(i,2);
    KERNEL(index) = interp1([1 find(index,1,'last')]', bandpasses(i,3:4)', find(index)', 'linear');
end
KERNEL = [KERNEL; KERNEL(end:-1:1)];
KERNEL = KERNEL + 1j * KERNEL;

kernel = real(ifft(KERNEL));
black = blackman(M);
black = [black(M/2+1:end); zeros(L-M,1); black(1:M/2)];
kernel = kernel .* black;

KERNEL = fft(kernel);

% Extended signal and kernel
KERNELlong = fft([kernel(1:L/2); zeros(2*M,1); kernel(L/2+1:end)]);
KERNELlong = repmat(KERNELlong, [1, size(signal,2), size(signal,3)]);
signallong = [repmat(signal(1,:,:), [M,1,1]); signal; repmat(signal(L,:,:), [M,1,1])];
i_signal = logical([zeros(M,size(signal,2),size(signal,3)); ones(size(signal)); zeros(M,size(signal,2),size(signal,3))]);

fftsignallong = fft(signallong);
filtsignallong = ifft(fftsignallong .* KERNELlong);
filtsignal = zeros(size(signal));
filtsignal(:) = filtsignallong(i_signal);

% Plotting
if strcmp(plotflag,'plot_on')
    figure
    plot(f, abs(KERNEL(1:L/2)))
    title('Frequency response of filter kernel')
    xlabel('Frequency [Hz]')
    ylabel('Normalized amplitude')

    before = 2 * abs(fft(signal(1:L)) / L);
    after  = 2 * abs(fft(filtsignal(1:L)) / L);
    figure
    plot(f, before(1:L/2), '-k', f, after(1:L/2), '-r')
    legend({'Unfiltered','Filtered'})
    title('Single sided amplitude spectrum of signal before and after filtering')
    xlabel('Frequency [Hz]')
    ylabel('Amplitude')
    xlim([0 f(end)])
    maxy = 10 * max(mean(before), mean(after));
    ylim([-0.01 * maxy 1.2 * maxy])

    figure
    plot(time, signal(1:L), '-b', time, filtsignal(1:L), '-r')
    title('Signal in time domain before and after filtering')
    legend({'Unfiltered','Filtered'})
    xlabel('Time [s]')
    ylabel('Signal units')
end
end