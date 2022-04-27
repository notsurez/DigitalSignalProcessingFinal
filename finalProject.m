clc; clear; close all;

% (1) Select a known time domain function, for example a chirp 
% waveform (Sect 2:10 of text), of moderate complexity.

%Using parametrs defined in case study 2 from 2.10
 %to show chirp (M = 512, f_s = 1E6, T = 1/f_s)
M = 512; f_s = 1E6; T = 1/f_s;
k = 1:M;
%f(k) and x(k) from the chirp function in 2.10
f_k = (k.*f_s)./(2*(M-1));
x_k = sin(2*pi.*f_k.*k*T);
%plot(k, x_k)


%title('Chirp Signal')
%xlabel("k")
%ylabel("x(k)")


% (2) Create MATLAB signals 1M (2^20 = 1,048,576) 
% samples long consisting of several chirps and noise.

num_samples = 1048576
k = 1:num_samples;
lowNoise = 40;
moderateNoise = 20;
highNoise = 10;
num_chirps = 5;

% awgn(signal, SNR) = add white gaussian noise
chirp_signal = zeros(1, num_samples);
chirp_signal = addRandomChirp(chirp_signal, num_chirps);

%low noise signal
ln_chirp = awgn(chirp_signal, lowNoise);

%moderate noise signal
mn_chirp = awgn(chirp_signal, moderateNoise);

%heavy noise signal
hn_chirp = awgn(chirp_signal, highNoise);

%plot 
%{
subplot(4, 1, 1)
plot(chirp_signal)
title('Chirp Signal (chirps = 5)')
subplot(4, 1, 2)
plot(ln_chirp)
title('Low Noise Chirp (SNR = 40)')
subplot(4, 1, 3)
plot(mn_chirp)
title('Moderate noise chirp (SNR = 20)')
subplot(4, 1, 4)
plot(hn_chirp)
title('High noise chirp = (SNR = 5)')
%}

r = f_corr(chirp_signal,x_k, 0, 0);
r_i = 1:length(r);
%get the maximum value of the cross correlation as well as the index
%where it occurs
peak_i = find(r > 1.5E-4)

subplot(2,1,1)
plot(chirp_signal)
title('Chirp signal no noise')
subplot(2,1,2)
plot(r)
hold on
plot(r_i, r,'-p', 'MarkerIndices', peak_i,'MarkerFaceColor','red', 'MarkerSize',10);
title('Linear correlation with noiseless chirp signal')


%Chirp adding function
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
