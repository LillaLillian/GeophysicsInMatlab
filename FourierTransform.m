%% Exercise set 6 - The Fourier Transform

%% 1. Seismogram 

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load seis.dat                   % Loading the seismogram

f=seis;                         % Defining the seismogram as f
t_int=length(f);                % Finding the length of the time vector
t=0:1:t_int-1;                  % Defining a sampling vector
sf=20;                          % Sampling frequency (Hz)
time=t/sf;                      % Defining a time vector (s)

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
plot(time,f)
title('Seismogram')
xlabel('Time (s)')
ylabel('Amplitude')
hold on 

% Plotting P arrrival at 596.9 s, and S arrival at 1063.6 s in red, 
% with width 1.5:
line([596.9 596.9],[0 25000],'Color','r','LineWidth',1.5)        % P wave 
line([1063.6 1063.6],[0 25000],'Color','r','LineWidth',1.5)      % S wave

% Labeling the two arrivals:
text(600, 25000, 'P')           % Labeling the P-arrival
text(1100, 25000, 'S')          % Labeling the S-arrival
hold off


%% 2. Discrete Fourier Transform

% INTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt=1/sf;                        % Defining delta t (s)
wc=2*pi/dt;                     % Defining the period of the signal 
N=length(f);                    % Defining N (amount of samples)
dw=wc/N;                        % Defining delta omega 
w=[0:dw:wc-dw];                 % Defining the angular frequency vector (Hz)

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftseis=fft(f);                  % Fourier transform of f
amp_ftseis=abs(ftseis);         % Amplitude spectrum (modulus)
ph_ftseis=angle(ftseis);        % Phase spectrum (argument of ftseis)


% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(20)
plot(w, real(ftseis))
title('Real part of F(w) (ftseis)')
xlabel('w (Hz)')

figure(21)
plot(w,imag(ftseis))
title('Imaginary part of F(w) (ftseis)')
xlabel('w (Hz)')

figure(22)
plot(w,amp_ftseis)
title('Amplitude spectrum of F(w) (ftseis)')
xlabel('w (Hz)')
ylabel('Amplitude')

figure(23)
plot(w,ph_ftseis)
title('Phase spectrum of F(w) (ftseis)')
xlabel('w (Hz)')
ylabel('Phase (rad)')

% The amplitude spectrum shows that the most of the energy is
% constrained to the lower frequencies.

%% 3. Time Shifting

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t01=150;                        % Time offset nr 1, (s)
t02=-245.3;                     % Time offset nr 2, (s)

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftseis=fft(f);                  % FT of f; time shifting is done in Fourier domain

sh1_FT= ftseis.*exp(-i*w*t01);  % Calculate time shift nr 1
sh2_FT= ftseis.*exp(-i*w*t02);  % Calculate time shift nr 2

sh1_f=ifft(sh1_FT);             % Inverse FT to find the shifted signals in 
sh2_f=ifft(sh2_FT);             % time

im1=imag(sh1_f);                % Imaginary part of both is close to 0
im2=imag(sh2_f);

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(30)
plot(time,f,time,real(sh1_f),'r')
title('Plot of original signal and the signal shifted by 150 s')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original signal','Signal shifted 150 s')

figure(31)
plot(time,f,time,real(sh2_f),'g')
title('Plot of original signal and the signal shifted by -245.3 s')
xlabel('Time (s)')
ylabel('Amplitude')
legend('Original signal','Signal shifted -245.3 s')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I used the command ginput to verify the result. The distance between
% the first arrival of the S wave was shifted 150 s further in time in the
% first plot and approximately 245.3 s earlier in time in the second plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4. Exploring the spectrum - dominant frequencies of P and S waves 

% 4a)

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tp1= 590;                           % Startime P (s)
tp2=650;                            % Endtime P (s)
ts1= 1050;                          % Startime S (s)
ts2=1110;                           % Endtime S (s)

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boxcarP=boxcar(f,time,sf,tp1,tp2); % Creating the boxcar function for P wave, with a function
f_P=f.*boxcarP;                    % Multiplying the boxcarP with the seismogram      

boxcarS=boxcar(f,time,sf,ts1,ts2); % Creating the boxcar function for S wave, with a function
f_S=f.*boxcarS;                    % Multiplying the boxcarS with the seismogram

ftseisP=fft(f_P);                  % FT of boxcar-isolated P arrival
ftseisS=fft(f_S);                  % FT of boxcar-isolated S arrival

amp_fP=abs(ftseisP);               % Amplitude spectrum of boxcar-isolated P arrival
amp_fS=abs(ftseisS);               % Amplitude spectrum of boxcar-isolates S arrival

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(40)
subplot(2,1,1)
plot(time,boxcarP)
axis([0 4000 0 2])                 % Changing axis so it looks better
title('Boxcar function, 590-650 s (P-wave)')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
plot(time,boxcarS)
axis([0 4000 0 2])                 % Changing axis so it looks better
title('Boxcar function, 1050-1110 s (S-wave)')
xlabel('Time (s)')
ylabel('Amplitude')

figure(41)
subplot(2,1,1)
plot(time,f_P)
title('P waves isolated by boxcar function')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fP)
title('Amplitude spectrum of boxcar-isolated P wave')
xlabel('w (Hz)')
ylabel('Amplitude')

figure(42)
subplot(2,1,1)
plot(time,f_S)
title('S waves isolated by boxcar function')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fS)
title('Amplitude spectrum of boxcar-isolated S wave')
xlabel('w (Hz)')
ylabel('Amplitude')

%% 4b) Cosine Taper

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tp1= 590;                            % Startime P (s)
tp2=650;                             % Endtime P (s)
ts1= 1050;                           % Startime S (s)
ts2=1110;                            % Endtime S (s)

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taperP=taper2(f,time,80,sf,tp1,tp2); % Creating costap for P, with function
f_P2=f.*taperP;                      % Appying the taperP to the signal

taperS=taper2(f,time,80,sf,ts1,ts2); % Creating costap for S, with function
f_S2=f.*taperS;                      % Appying the taperS to the signal

ftseisP2=fft(f_P2);                  % FT of taper-isolated P wave
ftseisS2=fft(f_S2);                  % FT of taper-isolated S wave

amp_fP2=abs(ftseisP2);               % Amplitude spectrum of taper-isolated P wave
amp_fS2=abs(ftseisS2);               % Amplitude spectrum of taper-isolated S wave

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(43)
subplot(2,1,1)
plot(time,taperP)
title('Taperfunction for P wave')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
plot(time,taperS)
title('Taperfunction for S wave')
xlabel('Time (s)')
ylabel('Amplitude')

figure(44)
subplot(2,1,1)
plot(time,f_P2)
title('P wave isolated by taper function')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fP2)
title('Amplitude spectrum of taper-isolated P wave')
xlabel('w (Hz)')
ylabel('Amplitude')

figure(45)
subplot(2,1,1)
plot(time,f_S2)
title('S wave isolated by taper function')
xlabel('Time (s)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fS2)
title('Amplitude spectrum of taper-isolated S wave')
xlabel('w (Hz)')
ylabel('Amplitude')

%% 4c) Spectrum 

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Isolated P spectrum
figure(46)
subplot(2,1,1)
plot(w,amp_fP)
title('Amplitude spectrum of boxcar-isolated P')
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fP2)
title('Amplitude spectrum of taper-isolated P')
xlabel('w (Hz)')
ylabel('Amplitude')

% Zoom in on isolated P spectrum
figure(47)
subplot(2,1,1)
plot(w,amp_fP)
title('Amplitude spectrum of boxcar-isolated P')
axis([0 10 0 6*10^5])
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fP2)
title('Amplitude spectrum of taper-isolated P')
axis([0 10 0 6*10^5])
xlabel('w (Hz)')
ylabel('Amplitude')

% Isolated S spectrum
figure(48)
subplot(2,1,1)
plot(w,amp_fS)
title('Amplitude spectrum of boxcar-isolated S')
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w,amp_fS2)
title('Amplitude spectrum of taper-isolated S')
xlabel('w (Hz)')
ylabel('Amplitude')

% Zoom in on isolated S spectrum
figure(49)
subplot(2,1,1)
plot(w,amp_fS)
title('Amplitude spectrum of boxcar-isolated S')
xlabel('w (Hz)')
ylabel('Amplitude')
axis([0 6 0 4*10^6])

subplot(2,1,2)
plot(w,amp_fS2)
title('Amplitude spectrum of taper-isolated S')
xlabel('w (Hz)')
ylabel('Amplitude')
axis([0 6 0 4*10^6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The spectrums of the taper-isolated waves seems to be the superior
% ones. They are more spikey and detailed and looks more similar to the 
% original spectrum, than the boxcar-isolated waves. The spectrum of the
% boxcar-isolated waves is more bubbly, and this is a common result from
% using the boxcar function. The reason for this difference is seen when 
% looking at the amplitude spectrum of the boxcar and taper function.
% The amplitude spectrum of the boxcar is a peak, and several smaller, 
% round peaks/lobes on both sides, while the amplitude spectrum of the 
% taper function have less of these side-lobes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4d. Dominant frequencies 

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

point_freP = find(amp_fP2==max(amp_fP2));   % Finding the point where amp is max
point_freS = find(amp_fS2==max(amp_fS2));   % for P spectrum and S spectrum

domf_P=w(point_freP(1))     % Picking out the point in the frequency vector, for P
domf_S=w(point_freS(1))     % Picking out the point in the frequency vector, for S

% P waves have generally higher frequencies than the S waves and this is
% reflected in the amplitude spectrum.

%% 5. Exploring the spectrum 

% a) Downsampling

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sf_down=0.5;                            % New, lower frequency (Hz)
dt_down=1/sf_down;                      % Defining new delta t (s) for the lower frequency
wc_down=2*pi/dt_down;                   % Defining the period of the downsampled signal

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdown=downsample(f,40);                 % Downsampling the signal f, by
                                        % choosing every 40 point = 0.5 Hz
ftseis_down=fft(fdown);                 % FT of the downsampled signal
time_down=downsample(time,40);          % Downsampling the time vector for plotting

N_down=length(fdown);                   % Calculating new N, for downsampled signal
dw_down=wc_down/N_down;                 % Calculating new delta omega 
w_down=[0:dw_down:wc_down-dw_down];     % Calculating the new angular frequency vector (Hz)

amp_ftseis=abs(ftseis);                 % Amp spectrum of original signal
amp_fdown=abs(ftseis_down);             % Amp spectrum of downsampled signal
norm_amp=amp_fdown*40;                  % Normalizing the downsampled spectrum

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(50)
plot(time_down,fdown)
title('Downsampled signal')
xlabel('Time (s)')
ylabel('Amplitude')

figure(51)
subplot(2,1,1)
plot(w,amp_ftseis)
title('Amplitude spectrum for orginal signal')
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w_down,norm_amp)
title('Amplitude spectrum for downsampled signal')
xlabel('w (Hz)')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The spectral content has changed so that it does not longer go to 0
% at the Nyquist frequency. This is beacuse the new, lower sampling 
% frequency gives a smaller Nyquist frequency. Since the spectrum doesn't
% go to 0 at Nyquist, we get aliasing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5b) Upsampling

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t2=0:1/6:t_int-1;                  % Defining an upsampled sampling vector
sf=20;                             % Sampling frequency (Hz)
sf_up=120;                         % Upsampled frequency (Hz)
time2=t2/sf;                       % Defining an upsampled time vector

dt_up=1/sf_up;                     % Defining delta t (s)
wc_up=2*pi/dt_up;                  % Defining the period of the signal

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_interp=interp1(time,f,time2,'spline');% Upsampling the signal f, using spline method

N_up=length(f_interp);             % Calculating new N, for upsampled signal
dw_up=wc_up/N_up;                  % Calculating new delta omega 
w_up=[0:dw_up:wc_up-dw_up];        % Calculating the new angular frequency vector (Hz)

ftseis_interp=fft(f_interp);       % FT of the upsampled signal
amp_finterp=abs(ftseis_interp);    % Amp spectrum of upsampled signal
n_amp_finterp=amp_finterp*(1/6);   % Normalizing the upsampled spectrum

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(52)
plot(time2,f_interp)
title('Upsampled signal')
xlabel('Time (s)')
ylabel('Amplitude')

figure(53)
subplot(2,1,1)
plot(w,amp_ftseis)
title('Amplitude spectrum of original signal')
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w_up,n_amp_finterp)
title('Amplitude spectrum of upsampled signal')
xlabel('w (Hz)')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The upsampling makes the Nyquist frequency higher, because the sampling
% frequency is higher. This means that the spectral content is even further
% away from the Nyquist frequency.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 5c) Padding the signal with zeroes

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sf=20;                              % Sampling frequency
tap=taper(200,sf,f,time);           % Creating a 200s range taper
f_tap=tap.*f;                       % Applying the taper to the signal
f_0=zeros(1,length(f)+5000);        % Creating vector of 0's with length of signal + 5000
f_0(1:length(f_tap))=f_tap;         % Adding the signal to the start of emtpy vector, 
                                    % leaving the last 5000 numbers as zero

t_int0=length(f_0);                 % Finding the length new time vector
t_0=0:1:t_int0-1;                   % Defining a new sampling vector
sf=20;                              % Sampling frequency (Hz)
time_0=t_0/sf;                      % Defining a new, longer time vector for plotting

dt0=1/sf;                           % Defining delta t (s)
wc_0=2*pi/dt0;                      % Defining the period of the signal
N_0=length(f_0);                    % Defining N (amount of samples)
dw_0=wc_0/N_0;                      % Defining delta omega 
w_0=[0:dw_0:wc_0-dw_0];             % Defining the angular frequency vector (Hz)

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ftseis=fft(f);                      % Fourier transform of original signal
F_0=fft(f_0);                       % FT of the longer f signal

amp_ftseis=abs(ftseis);             % Amplitude spectrum of original signal
amp_F_0= abs(F_0);                  % Amplitude spectrum of longer signal

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(54)
plot(time_0,f_0)
title('Signal with 5000 0`s added')
xlabel('Time (s)')
ylabel('Amplitude')

figure(55)
subplot(2,1,1)
plot(w,amp_ftseis)
title('Amplitude spectrum of original signal')
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w_0,amp_F_0)
title('Amplitude spectrum of signal with added 0`s')
xlabel('w (Hz)')
ylabel('Amplitude')

figure(56)
subplot(2,1,1)
plot(w,amp_ftseis)
axis([0.5 0.7 0 10^7])
line([0.508 0.508],[3*10^6 9*10^6],'Color','r','LineWidth',1.5)
title('Zoom in on amplitude spectrum of original signal')
xlabel('w (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
plot(w_0,amp_F_0)
axis([0.5 0.7 0 10^7])
line([0.508 0.508],[3.5*10^6 9*10^6],'Color','r','LineWidth',1.5) 
title('Zoom in on amplitude spectrum of signal with added 0`s')
xlabel('w (Hz)')
ylabel('Amplitude')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding zeros to the signal means you'll get more points (N gets bigger).
% This does not add any information, but the spectral resolution is higher.
% This can be seen by the higher detail that's in the new amplitude
% spectrum. I have pointed out one example with a red line, where an extra
% spike is visible on the new ampl spectrum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6. The sampling theorem 

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load mystery.mat;            % Loading the mystery file
sf2=10;                      % Sampling frequency
N2=1000;                     % Number of points
t2=0:N2-1;                   % Sampling points
time2=t2/sf2;                % Defining time vector

% COMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_flip=conj(fliplr(myst(2:length(myst)-1)));   % Flipping and taking 
                                               % conjugate of signal, 
                                               % minus first and last point
                                               
myst_F= [myst m_flip];       % Putting original myst and flipped myst 
                             % together into one signal
myst_f=ifft(myst_F);         % Taking the inverse FT of signal

% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(60)
plot(time2,myst_f)
title('Mystery time series')
xlabel('Time (s)')
ylabel('Amplitude')
