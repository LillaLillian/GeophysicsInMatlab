function [ boxcar ] = boxcar( signalf, time, samplingfrequency, starttime, endtime )
% Creating a boxcar function
%
% signalf = seismogram
% time = time vector
% samplingfrequency = sampling frequency
% starttime = time of start of boxcar = 1
% endtime = time of end of boxcar = 1
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np=length(signalf);

for n=1:np;
    if n < starttime*samplingfrequency;   % Multiply with sampling frequency
                                          % to get it in samples instead of sec 
                            
        boxcar(n) = 0;                    % For all n < 590 s the boxcar is 0
        
    elseif n > endtime*samplingfrequency; % For all n > 650 s the boxcar is 0
        boxcar(n)=0;
        
    else
        boxcar(n)=1;                      % For all other values of n the boxcar is 1
        
    end
end

end

