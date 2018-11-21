%% 5. exercise set

%% 1. Complex numbers

% INPUT
a = 3.5;
b = -7.8;
z1 = a+b*i;

% COMP
rz1=real(z1);                   % Real part of z1
imz1=imag(z1);                  % Imaginary part of z1

% OUTPUT
figure(10)
plot(rz1,imz1,'*')              % Plotting z1 in the complex plane
axis([-30 30 -30 30])
hold on
plot([-30 30], [0 0])           % Create x-axis on plot
plot([0 0], [-30 30])           % Create y-axis on plot

xlabel('Re')        
ylabel('Im')
title('Plot of z1 in the complex plane')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modulus and argument

% By definition:
'Modulus of z1 calculated from definition:'
mod_z1=sqrt(a^2 + b^2)   
'Argument of z1 calculated from definition:'
arg_z1=atan(b/a)                

% By matlab commands
'Modulus of z1 calculated with matlab command:'
mod2_z1 = abs(z1)              
'Argument of z1 calculated with matlab command:'
arg2_z1= angle(z1)             

% On the plot the modulus of z1 is the length from (0,0) to z1, and the
% argument is the angle between the real axis and the line from (0,0) to z1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Complex conjugate of z1, z1_con

% INPUT
'Complex conjugate, z1*'
z1_con=a-b*i                    % Complex conjugate of z1


% COMP
rz1_con = real(z1_con);         % Real part of the conjugate
imz1_con = imag(z1_con);        % Imaginary part of the conjugate

% Modulus and argument
'Modulus of z1* calculated with matlab command:'
mod2_z1_con = abs(z1_con) 
'Argument of z1* calculated with matlab command:'
arg2_z1_con= angle(z1_con)      


% OUPUT
figure(11)
plot(rz1,imz1,'*',rz1_con,imz1_con,'r*')      % Plotting z1* and z1 in the complex plane
axis([-30 30 -30 30])
hold on
plot([-30 30], [0 0])           % Create x-axis on plot
plot([0 0], [-30 30])           % Create y-axis on plot

xlabel('Re')        
ylabel('Im')
title('Plot of z1 and z* in the complex plane')
hold off
legend('z1','z1*')

% The modulus of the complex conjugate, z1*, is the same as for z1, and the
% argument z1* has the same value, but the opposite sign of that of z1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1b.

%INPUT
z2=25.3*exp(i*5.1);

%COMP
rz2=real(z2);
imz2= imag(z2);

% Modulus and argument of z2
'Modulus of z2 calculated with matlab command:'
mod2_z2 = abs(z2)               
'Argument of z2 calculated with matlab command:'
arg2_z2= angle(z2)             

%OUTPUT
figure(12)
plot(rz2,imz2,'*')
axis([-30 30 -30 30])
hold on
plot([-30 30], [0 0])           % Create x-axis on plot
plot([0 0], [-30 30])           % Create y-axis on plot
xlabel('Re')        
ylabel('Im')
title('Plot of z2 in the complex plane')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On the form z2=a + ib, z2 is written as z2 = 9.6 - 23.4i. Where a = 9.6
% and b = -23.4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Seismogram

% 2a)
% INPUT
load col_z.txt                  % Loading the seismogram
t_int=0:(length(col_z)-1);      % Time intervals
freq=20;                        % Sampling frequency
t=t_int/freq;                   % Sampling time
f=col_z;                        % Denoting the seismogram as f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT

figure(20)
plot(t,f)
ylabel('Amplitude')
xlabel('Time (s)')
title('Seismogram of Mariana Island earthquake')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2b) Removing the mean

%COMP
fmean=mean(f);                  % Finding the mean of the seismogram
fm=f-fmean;                     % Removing the mean from the seismogram to get the mean of f = 0

%OUTPUT
figure(21)
subplot(3,1,1)
plot(t,f)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, f(t)')

subplot(3,1,2)
plot(t,fm,'r')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Timeshifted seismogram, mean=0, fm(t) ')

subplot(3,1,3)
plot(t,f,t,fm,'r')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Comparison between f(t) and fm(t)')
legend('f(t)','fm(t)')

% The timeseries f(t) is shifted downwards when the mean is removed from
% the signal. The modified signal, fm(t), is centered around 0 on the
% y-axis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2c) Detrending the signal

%COMP
fd=detrend(f);                  % Removing trends from the signal

%OUTPUT
figure(22)
subplot(3,1,1)
plot(t,f)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, f(t)')

subplot(3,1,2)
plot(t,f,t,fd,'g')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Comparison between f(t) and fd(t)')
legend('f(t)','fd(t)')

subplot(3,1,3)
plot(t,fm,'r',t,fd,'g')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Comparison between fm(t) and fd(t)')
legend('fm(t)','fd(t)')

% When using the detrend command the trends in the seismogram is removed.
% The modified signal, fd(t), is centered around 0 on the y-axis.
% In addition, an upward-going trend is removed, and the signal is
% therefore oriented more along the horizontal axis, unlike fm(t) which is
% sloping upwards with time. When using the detrend command on f(t), it is
% not necessary to remove the mean first as we did in 2b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2d) Downsampling script

%INPUT
a=length(fd);                   % Length of the signal

%COMP
d=20/0.5;                       % Ratio between orginal and new sampling
fd1=fd(1:d:a);                  % Picking out every 40 point (ratio =40)
t1=t(1:d:a);                    % Creating a time vector with every 40 point

%OUTPUT
figure(23)
subplot(2,1,1)
plot(t,fd)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, fd(t)')

subplot(2,1,2)
plot(t1,fd1)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Downsampled seismogram of Mariana Island earthquake, fd1(t)')

figure(24)
subplot(2,1,1)
plot(t,fd,'*-')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
axis([1100 1120 -1.5*10^4 1.5*10^4])
title('Zoom in on original seismogram, fd(t)')

subplot(2,1,2)
plot(t1,fd1,'*-')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
axis([1100 1120 -1.5*10^4 1.5*10^4])
title('Zoom in on downsampled seismogram, fd1(t)')

% When looking at the whole signal it is difficult to see any difference
% between the orignal signal and the downsampled one. However, when looking
% at a zoomed in section the difference is great. The signals resolution is
% greatly reduced. The smooth waves from the original signal is reduced to
% a jagged signal. The points marking the samples are also greatly reduced.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Downsampling using downsample command

%COMP
fd_ds=downsample(fd,40);        % The downsample command picks every 40 point from fd
t_ds=downsample(t,40);          % Need to downsample the time vector for plotting

%OUTPUT
figure(25)
subplot(3,1,1)
plot(t,fd)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, fd(t)')

subplot(3,1,2)
plot(t1,fd1)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Downsampled seismogram of Mariana Island earthquake, fd1(t)')

subplot(3,1,3)
plot(t_ds,fd_ds)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Downsampled seismogram, fd2(t), using the downsample command')

figure(26)
subplot(3,1,1)
plot(t,fd,'*-')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
axis([1100 1120 -1.5*10^4 1.5*10^4])
title('Zoom in on original seismogram, fd(t)')

subplot(3,1,2)
plot(t1,fd1,'*-')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
axis([1100 1120 -1.5*10^4 1.5*10^4])
title('Zoom in on downsampled seismogram, fd1(t)')

subplot(3,1,3)
plot(t1,fd_ds,'*-')
grid on
ylabel('Amplitude')
xlabel('Time (s)')
axis([1100 1120 -1.5*10^4 1.5*10^4])
title('Zoom in on downsampled seismogram, fd2(t)')

% The downsampled signal resulting from the downsample command seems to be
% similar to the one resulting from my script.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2d) Upsampling

%INPUT
a2=length(fd);                  % Length of fd
freq2=80;                       % New, upsampled frequency
c=(length(fd))*4;               % Increasing the amount of samples by a factor of 4 (20*4=80Hz)
to=0:c;                         % Time increments increased by a factor of 4
t2=to/freq2;                    % New time vector 

%COMP
fd3=interp1(t,fd,t2);           % Using the interp1 command to interpolate fd
fd4=interp1(t,fd,t2,'spline');  % Using the interp1 command to interpolate fd, with the spline method

%OUTPUT
figure(27)
subplot(3,1,1)
plot(t,fd)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, 20Hz, fd(t)')

subplot(3,1,2)
plot(t2,fd3)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Upsampled seismogram, using linear interpolation, 80Hz, fd3(t)')

subplot(3,1,3)
plot(t2,fd4)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Upsampled seismogram, using spline interpolation, 80Hz, fd4(t)')

% Zoomed in sections
figure(28)
subplot(3,1,1)
plot(t,fd,'*-')
axis([1100 1110 -10^4 10^4])
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram, fd(t), zoomed in')

subplot(3,1,2)
plot(t2,fd3,'*-')
axis([1100 1110 -10^4 10^4])
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Upsampled seismogram, fd3(t), zoomed in')

subplot(3,1,3)
plot(t2,fd4,'*-')
grid on
axis([1100 1110 -10^4 10^4])
ylabel('Amplitude')
xlabel('Time (s)')
title('Upsampled seismogram, fd4(t), zoomed in ')

% There is no visible change in the signal after upsampling. Even after
% zooming in the signal looks the same. The difference is only seen when marking
% the samples with *, then the increased amount of samples can be seen.
% Even so, the original signal is already very densly sampled and the shape
% of the signal does not really change with the increased amount of
% samples.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2e) Taper function

clear all
close all

%INPUT
% Preparing the signal
load col_z.txt                  % Loading the seismogram
ft=col_z;                       % Defining f(t)
fd=detrend(ft);                 % Detrending the signal
freq=20;                        % Frequency
np=length(fd);                  % Length fd
t_int=0:(length(fd)-1);         % Time intervals
t=t_int/freq;                   % Sampling time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Input for Taper1, 3s
range1=3;                       % 3 seconds range
s_range1=range1*freq;           % Sampling range (3 seconds * amount of samples per sec)
theta1=linspace(0,pi,s_range1); % Defining theta, 0-pi, over the length of the 3 s range

costap1=zeros(1,np);            % Creating empty array for the taper1 function

counter1=0;                     % Counter to be used in the loop

%COMP1
for i=1:np;

    if i < s_range1;            % For the first 500 s of the taper, the formula used for costap is:
        costap1(i) = 0.5*(-cos(theta1(i))+1);
        
    elseif i > np-(s_range1);   % For the last 500 s of the taper, the formula used for costap is:
        counter1=counter1+1;      % (Using counter in place of i, in the second part of the loop)
        costap1(i) = 0.5*(cos(theta1(counter1))+1);
        
    else costap1(i)=1;           % For all other times the taper is set to 1
    end
end

fd_tap1=fd.*costap1;              % Taper1 applied to the signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input for Taper2, 500s

range2=500;                     % 500 seconds range
s_range2=range2*freq;           % Sampling range
theta2=linspace(0,pi,s_range2); % Defining theta, 0-pi, over the length of the 500s range

costap2=zeros(1,np);            % Creating empty array for the taper2 function

counter2=0;                     % Counter to be used in the loop

%COMP2
for i=1:np;

    if i < s_range2;             % For the first 500 s of the taper, the formula used for costap is:
        costap2(i) = 0.5*(-cos(theta2(i))+1);
        
    elseif i > np-(s_range2);    % For the last 500 s of the taper, the formula used for costap is:
        counter2=counter2+1;      % (Using counter in place of i, in the second part of the loop)
        costap2(i) = 0.5*(cos(theta2(counter2))+1);
        
    else costap2(i)=1;           % For all other times the taper is set to 1
    end
end

fd_tap2=fd.*costap2;              % Taper applied to the signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
% Plotting the 3 s taper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(29)
subplot(4,1,1)
plot(t,costap1)
grid on
title('Taper function, range = 3s')
ylabel('Amplitude')
xlabel('Time (s)')

subplot(4,1,2)
plot(t,fd)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, 20Hz, fd(t)')

subplot(4,1,3)
plot(t,fd_tap1)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('The 3s taper function applied to the original seismogram')

subplot(4,1,4)
plot(t,fd_tap1)
axis([0 20 -10^4 10^4])
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Zoom in on the seismogram with applied taper')

% Plotting the 500s taper %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(30)
subplot(3,1,1)
plot(t,costap2)
grid on
title('Taper function, range = 500s')
ylabel('Amplitude')
xlabel('Time (s)')

subplot(3,1,2)
plot(t,fd)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('Original seismogram of Mariana Island earthquake, 20Hz, fd(t)')

subplot(3,1,3)
plot(t,fd_tap2)
grid on
ylabel('Amplitude')
xlabel('Time (s)')
title('The 500s taper function applied to the original seismogram')

% The 3s taper has a very short range and it is not possible to see the
% result of applying it to the signal, without zooming in. When zoomed in
% you can see that the signal is now starting at zero and increasing during 
% the first 3 seconds. At the 3 last seconds of the signal it tapers of to 
% zero again.  
% The 500s taper har a longer range, and the effect of applying it to the
% signal is easy to see, even without zooming. The signal starts at zero
% and increases for the first 500 seconds, and tapers off to zero during 
% the last 500s. 
% The signal between the ranges are left unchanged.































