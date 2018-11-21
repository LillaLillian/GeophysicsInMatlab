function [ costap ] = taper2( signalf, time, range, samplingfrequency, starttime, endtime )
%Creating taper function for shorter intervals
%  The taper can start at other point than the start of the signal
%
% signalf = seismogram
% time = time vector
% samplingfrequency = sampling frequency
% starttime = time of start of boxcar = 1
% endtime = time of end of boxcar = 1
% range = length of taper interval 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np=length(signalf);
costap=zeros(1,np);                 % Creating empty array for the taper function
t1=starttime;                       % Start of interval where taper = 1
t2=endtime;                         % End of interval where taper = 1
srange=range*samplingfrequency;     % Range in frequency
theta=linspace(0,pi,srange);        % Defining theta

int_1=find(time>=t1 & time<=t2);    % Filling in the interval with 1's
costap(int_1)=1;

int_start=find(time>(t1-range) & time<t1);  % Creating the interval going from 0 to 1 (start)
int_end=find(time>(t2) & time<t2+range);    % Creting the interval going from 1 to zero (end)

counter=0;                          % counter for use in loop 1

for i=int_start;                    % The interval from 0-1 is defined 
    counter=counter+1;              % by this formula:
    costap(i) = 0.5*(-cos(theta(counter))+1);    
end

counter2=0;                         % counter for use in loop 2
for j=int_end;                      % The interval from 1-0 is defined 
     counter2=counter2+1;           % by this formula:
     costap(j) = 0.5*(cos(theta(counter2))+1);
end

end

