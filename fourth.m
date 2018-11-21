%% Fourth Exercise set
% Lillian Jensen

%% Acoustic wave equation script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-2 DF Soln to ACOUSTIC WAVE EQUATION (cnst density)
% (nx,nz,nt)      - input- (Horizontal,Vertical) gridpt dimens. of vel 
%                            model & # Time Steps
% FRE             - input- Peak frequency of Ricker wavelet
% BVEL            - input- NXxNZ matrix of background velocity model
% (dx,dt)         - input- (space, time) sample intervals
% (xs,zs)         - input- (x,z) coordinates of line source
% RICKER(nt)      - input- nt vector of source time histories
% (p2,p1,p0)      -calcul- (future,present,past) NXxNZ matrices of 
%                           modeled pressure field
% (p0,p1)         -output- Old and present pressure panels at time nt.
% REALDATA(nx,nt) -output- CSG seismograms at z=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

c=4.0;              % Velocity
FRE=20;             % Source freq. Hz 
nx=300;             % Define model size, size in x direction
nz=nx;              % Size in z direction. Model is 300x300
dx=c/FRE/20;        % Set spatial sampling dx

dt=.5*dx/c;         % Set dt for stability condition, time steps
xs=round(nx/2.3);   % Source coordinates, x
zs=round(nx/2);     % Source coordinates, z

nt=600;             % Length t
nxm=nx-1;           % Boundary limits
nzm=nz-1;           % Boundary limits

%( COMP 1 )
t = [0:1:nt-1]*dt-0.95/FRE;     % Time array
RICKER = zeros(length(t));      % Initiate Ricker wavelet
RICKER = (1-t .*t * FRE^2 *pi^2  ).*exp(- t.^2 * pi^2 * FRE^2 ) ; %Calculate Ricker = Source wavelet

%( OUTPUT 1 )
figure(21)
plot([0:nt-1]*dt,RICKER);
title('Ricker Wavelet');
xlabel('Time (s)')

pause(3)    

% INPUT cont.
% Setting up velocity model
BVEL=ones(nx,nz)*c; % Matrix of 1's, with size 300x300, multiplied with velocity = matrix of 4's
% BVEL(nx-round(nx/2):nx,:)= BVEL(nx-round(nx/2):nx,:)*1.4; % 2-Layer velocity model, 
                                                          % Changing the velocity model to create a layer

% 4 layers velocity model
BVEL(1:25,:)=3.8;
BVEL(26:100,:)=4.9;
BVEL(101:200,:)=3.5;
BVEL(201:300,:)=4.9;

REALDATA=zeros(nx,nt); % Define shot gather matrix (ERROR - ntt should be nt)
p0=zeros(nx,nz);       % Size of velocity model
p1=p0;                 % ---------||-----------
p2=p0;                 % ---------||-----------
cns=4*(dt/dx*BVEL).^2; % 
recdepth=7;            % Receiver line is 7 grid pts below source, Receriver Depth

figure(22)

% MAIN COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOOP OVER TIME STEPS
for it=1:1:nt     % When it = 1, all terms in p2 is zero, and all terms 
                  % except RICKER (it) is zero in p2(xs,zs)
                  
   p2 = 2*p1 - p0 +  cns.*del2(p1);    % FD of Acoustic Wave Equation
   p2(xs,zs) = p2(xs,zs) + RICKER(it); % Add bodypoint src term
   a = p2(xs+recdepth,:); 
   REALDATA(:,it) = a';                 % Save seismograms
   alpha=dt/dx;
   
%  ABC BOUNDARY CONDITIONS
%  p2(1,:)=p1(1,:)   + alpha*BVEL(1,:).*(p1(2,:)-p1(1,:));      % Top ABC
   p2(nx,:)= p1(nx,:) - alpha*BVEL(nx,:).*(p1(nx,:)-p1(nxm,:)); % Bottom ABC; stops the wave from reflecting 
                                                                % of bottom
   p2(:,nz)= p1(:,nz) + alpha*BVEL(:,nz).*(p1(:,nzm)-p1(:,nz)); % RHS ABC; stops the wave from 
   p2(:,1) = p1(:,1) + alpha*BVEL(:,1).*(p1(:,2)-p1(:,1));      % reflecting at the sides

   p0=p1;
   p1=p2; % UPDATE PRESSURE FIELDS
   
   % OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if round(it/20)*20==it;
       plt;     %Script 
       pause(1);
   end % PLOT OUT RESULTS
end;
