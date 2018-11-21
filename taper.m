function [ costap ] = taper( range, freq, signalf, time  )
%Creating a taper function, which starts at beginning of signal and ends at
%the end of the signal
%
% signalf = seismogram
% time = time vector
% freq = sampling frequency
% range = length of taper interval 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np=length(signalf);            % Length of signal 
s_range=range*freq;            % Sampling range (range (s) * amount of samples per sec)
theta=linspace(0,pi,s_range);  % Defining theta, 0-pi, over the length of the s_range

costap=zeros(1,np);            % Creating empty array for the taper function

counter=0;                     % Counter to be used in the loop

%COMP1
for i=1:np;

    if i < s_range;            % For the first 500 s of the taper, the formula used for costap is:
        costap(i) = 0.5*(-cos(theta(i))+1);
        
    elseif i > np-(s_range);   % For the last 500 s of the taper, the formula used for costap is:
        counter=counter+1;     % (Using counter in place of i, in the second part of the loop)
        costap(i) = 0.5*(cos(theta(counter))+1);
        
    else costap(i)=1;          % For all other times the taper is set to 1
    end
end

end

