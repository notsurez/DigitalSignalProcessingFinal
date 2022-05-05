clc; clear; close all;

% waveform (Sect 2:10 of text), of moderate complexity.

%Using parametrs defined in case study 2 from 2.10
 %to show chirp (M = 512, f_s = 1E6, T = 1/f_s)
M = 512; f_s = 1E6; T = 1/f_s;
k = 1:M;
%f(k) and x(k) from the chirp function in 2.10
f_k = (k.*f_s)./(2*(M-1));
x_k = sin(2*pi.*f_k.*k*T);

num_samples = 1048576;

num_chirps = 8;
num_sines = 3;

lowNoise = 30;
moderateNoise = 10;
highNoise = 1;

SineFreqs = zeros(1,num_sines);
Signal = zeros(1, num_samples);
SineFreqs = createRandomSineFreqs(SineFreqs, num_sines);
chirp_signal = addRandomChirp(Signal, num_chirps);
chirp_sine_signal = addRandomSines(chirp_signal, num_sines, SineFreqs);

subplot(4,2,1)
plot(chirp_signal)
title('5 Chirps - No Noise')
ylabel('x(k)')
xlabel('k')
xlim([0 num_samples])

subplot(4,2,2)
plot(chirp_sine_signal)
title('5 Chirps & 3 sines - No Noise - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])

ln_chirp_signal = awgn(chirp_signal, lowNoise);
ln_chirp_sine_signal = awgn(chirp_sine_signal, lowNoise);

subplot(4,2,3)
plot(ln_chirp_signal)
title('5 Chirps - Low Noise')
ylabel('x(k)')
xlabel('k')
xlim([0 num_samples])

subplot(4,2,4)
plot(ln_chirp_sine_signal)
title('5 Chirps & 3 Sines - Low Noise - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])

mn_chirp_signal = awgn(chirp_signal, moderateNoise);
mn_chirp_sine_signal = awgn(chirp_sine_signal, moderateNoise);

subplot(4,2,5)
plot(mn_chirp_signal)
title('5 Chirps - Moderate Noise')
ylabel('x(k)')
xlabel('k')
xlim([0 num_samples])

subplot(4,2,6)
plot(mn_chirp_sine_signal)
title('5 Chirps & 3 Sines - Moderate Noise - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])

hn_chirp_signal = awgn(chirp_signal, highNoise);
hn_chirp_sine_signal = awgn(chirp_sine_signal, highNoise);

subplot(4,2,7)
plot(hn_chirp_signal)
title('5 Chirps - High Noise')
ylabel('x(k)')
xlabel('k')
xlim([0 num_samples])

subplot(4,2,8)
plot(hn_chirp_sine_signal)
title('5 Chirps & 3 Sines - High Noise - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])
pause;

j = 1;
sigs_to_correlate = [chirp_signal; ln_chirp_signal; mn_chirp_signal; hn_chirp_signal];
for i = 1:(size(sigs_to_correlate, 1))
    current_plot = sigs_to_correlate(i,:);
    %Correleate the signal using dsp_companion
    r = f_corr(current_plot,x_k, 0, 0);
    r_i = 1:length(r);

    %get the maximum value of the cross correlation as well as the index
    %where it occurs
    peak_i = find(r > 1.5E-4);
    
    if(i == 2)
        j = 3;
    end
    if(i == 3)
        j = 5;
    end
    if(i == 4)
        j = 7;
    end

    subplot(4,2,j)
    plot(current_plot)
    title('Signal with Chirps')
    ylabel('x(k)')
    xlabel('k')
    xlim([0 num_samples])
    subplot(4,2,j+1)
    plot(r)
    hold on
    plot(r_i, r,'-p', 'MarkerIndices', peak_i,'MarkerFaceColor','red', 'MarkerSize',10);
    title('Linear correlation')
    ylabel('r_{xy}')
    xlabel('k')
    xlim([0 num_samples])
end
pause;

j = 1;
sigs_to_correlate = [chirp_sine_signal; ln_chirp_sine_signal; mn_chirp_sine_signal; hn_chirp_sine_signal];
for i = 1:(size(sigs_to_correlate, 1))
    current_plot = sigs_to_correlate(i,:);
    %Correleate the signal using dsp_companion
    r = f_corr(current_plot,x_k, 0, 0);
    r_i = 1:length(r);

    %get the maximum value of the cross correlation as well as the index
    %where it occurs
    peak_i = find(r > 1.5E-4);
    
    if(i == 2)
        j = 3;
    end
    if(i == 3)
        j = 5;
    end
    if(i == 4)
        j = 7;
    end

    subplot(4,2,j)
    plot(current_plot)
    title('Signal with Chirps & Sines')
    ylabel('x(k)')
    xlabel('k')
    xlim([0 num_samples])
    subplot(4,2,j+1)
    plot(r)
    hold on
    plot(r_i, r,'-p', 'MarkerIndices', peak_i,'MarkerFaceColor','red', 'MarkerSize',10);
    title('Linear correlation')
    ylabel('r_{xy}')
    xlabel('k')
    xlim([0 num_samples])
end
pause;

TFofFIRFilterBank = createWindowFilters(SineFreqs, num_sines);
[FIRN,FIRD] = tfdata(TFofFIRFilterBank);
FIRN = cell2mat(FIRN);
FIRD = cell2mat(FIRD);
fvtool(FIRN,FIRD)

TfofIIRFilterBank = createNotchFilters(SineFreqs, num_sines);
[IIRN,IIRD] = tfdata(TfofIIRFilterBank);
IIRN = cell2mat(IIRN);
IIRD = cell2mat(IIRD);
fvtool(IIRN,IIRD)
pause;

subplot(3,1,1)
plot(chirp_sine_signal)
title('Unfiltered - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])

FIR_chirp_sine_signal = filter(FIRN,FIRD,chirp_sine_signal);

subplot(3,1,2)
plot(FIR_chirp_sine_signal)
title('FIR Filtered - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])

IIR_chirp_sine_signal = filter(IIRN,IIRD,chirp_sine_signal);

subplot(3,1,3)
plot(IIR_chirp_sine_signal)
title('IIR Filtered - First 1000 Samples Only')
ylabel('x(k)')
xlabel('k')
xlim([0 1000])
pause;

%(4.c)Varying FFTs of the Signal
FFTSignal=chirp_signal;
    %2097152pt
        DoubleFFT=fft(FFTSignal,2097152);
    %1048576pt
        StraightFFT=fft(FFTSignal);
    %524288pt
        HalfFFT_pt1= fft(FFTSignal,524288);
        HalfFFT_pt2= fft(FFTSignal(524289:end),524288);
        HalfFFT=HalfFFT_pt1+HalfFFT_pt2;
    %1024pt
        TenTwentyFourFFT=fft(FFTSignal,1024);
        for i=1:(length(FFTSignal)/1024)-1
            TenTwentyFourFFT_RT=fft(FFTSignal(i*1024:(i+1)*1024),1024);
            TenTwentyFourFFT=TenTwentyFourFFT+TenTwentyFourFFT_RT;
        end            
    %32pt
        ThirtyTwoFFT=fft(FFTSignal,32);
        for i=1:(length(FFTSignal)/32)-1
            ThirtyTwoFFT_RT=fft(FFTSignal(i*32:(i+1)*32),32);
            ThirtyTwoFFT=ThirtyTwoFFT+ThirtyTwoFFT_RT;
        end    
 %Plotting The FFTs
        subplot(5,2,1)%2097152pt
        stem(DoubleFFT)
        title('2097152pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 2097152])
        subplot(5,2,3)%1048576pt
        stem(StraightFFT)
        title('1048576pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 1048576])
        subplot(5,2,5)%524288pt    
        stem(HalfFFT)
        title('524288pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 524288])
        subplot(5,2,7)%1024pt
        stem(TenTwentyFourFFT)
        title('1024pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 1024])
        subplot(5,2,9)%32pt
        stem(ThirtyTwoFFT)
        title('32pt FFT')  
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 32])

%(4.c)Varying FFTs of the Signal
FFTSignal=chirp_sine_signal;
    %2097152pt
        DoubleFFT=fft(FFTSignal,2097152);
    %1048576pt
        StraightFFT=fft(FFTSignal);
    %524288pt
        HalfFFT_pt1= fft(FFTSignal,524288);
        HalfFFT_pt2= fft(FFTSignal(524289:end),524288);
        HalfFFT=HalfFFT_pt1+HalfFFT_pt2;
    %1024pt
        TenTwentyFourFFT=fft(FFTSignal,1024);
        for i=1:(length(FFTSignal)/1024)-1
            TenTwentyFourFFT_RT=fft(FFTSignal(i*1024:(i+1)*1024),1024);
            TenTwentyFourFFT=TenTwentyFourFFT+TenTwentyFourFFT_RT;
        end            
    %32pt
        ThirtyTwoFFT=fft(FFTSignal,32);
        for i=1:(length(FFTSignal)/32)-1
            ThirtyTwoFFT_RT=fft(FFTSignal(i*32:(i+1)*32),32);
            ThirtyTwoFFT=ThirtyTwoFFT+ThirtyTwoFFT_RT;
        end    
 %Plotting The FFTs
        subplot(5,2,2)%2097152pt
        stem(DoubleFFT)
        title('2097152pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 2097152])
        subplot(5,2,4)%1048576pt
        stem(StraightFFT)
        title('1048576pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 1048576])
        subplot(5,2,6)%524288pt    
        stem(HalfFFT)
        title('524288pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 524288])
        subplot(5,2,8)%1024pt
        stem(TenTwentyFourFFT)
        title('1024pt FFT')
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 1024])
        subplot(5,2,10)%32pt
        stem(ThirtyTwoFFT)
        title('32pt FFT')  
        ylabel('X(f/f_s))')
        xlabel('f/f_s')
        xlim([0 32])

function IIRNotch = createNotchFilters(Freqs, num_sines)
    f_s = 1E6;
    wo = zeros(1,num_sines);
    bw = zeros(1,num_sines);
    FilterBank = filt([1],[1]);
    for i = 1:num_sines
        wo(i) = Freqs(i)/(f_s/2);
        bw(i) = wo(i)/35;
        [b,a] = iirnotch(wo(i),bw(i));
        filters(i) = filt(b,a);
        FilterBank = FilterBank*filters(i);
    end
    IIRNotch = FilterBank;
end

function FIRFilter = createWindowFilters(freqs, num_sines)
    f_s = 1E6;
    n = 100;
    Wn_low = zeros(1,num_sines);
    Wn_high = zeros(1, num_sines);
    FIRBank = filt([1],[1]);
    for i = 1:num_sines
        Wn_low(i) = (freqs(i)-10000)/(f_s/2);
        Wn_high(i) = (freqs(i)+10000)/(f_s/2);
        b = fir1(n,[.000001 Wn_low(i) Wn_high(i) .999999]);
        filters(i) = filt(b,[1]);
        FIRBank = FIRBank*filters(i);
    end
    FIRFilter = FIRBank;
end

function FreqsOfSines = createRandomSineFreqs(x, num_sines)
    f_s = 1E6;
    for i = 1:num_sines
        x(i) = rand*f_s/4;
    end
    FreqsOfSines = x;
end

%Function to add sin waves to a given signal,
function sigWithSines = addRandomSines(x, num_sines, Freqs)
    M = length(x); f_s = 1E6; T = 1/f_s;

    x_sine = zeros(1, M);
    Amp = zeros(1,num_sines);
    for i = 1:num_sines
        Amp(i) = rand*1;
    end
    for j = 1:num_sines
        Freqs(j) = rand*f_s/10;
        for i = 1:M
            x_sine(i) = x_sine(i) + Amp(j)*sin(2*pi*Freqs(j)*i*T);
        end
    end
    x = x + x_sine;
    sigWithSines = x;
end

function sigWithChirp = addRandomChirp(x, num_chirps)
    L = length(x); %length of input signal

    %Generate Chirp Signal
    M = 512; f_s = 1E6; T = 1/f_s;
    k = 1:M;
    %f(k) and x(k) from the chirp function in 2.10
    f_k = (k.*f_s)./(2*(M-1));
    x_k = sin(2*pi.*f_k.*k*T);

    %Generate desired number of chirps in signal
    for j = 1:num_chirps
        %Generate add chirp signal to a random point in x, make sure it fits
        chirpStart = floor(random('Uniform', 0, L-M));
        for i = 1:M
            x(i + chirpStart) = x_k(i);
        end
    end
    sigWithChirp = x; %return argument
end