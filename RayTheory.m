%% Ray tracing in simple velocity models, defined by equations

%% 1. 1D Velocity model

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=3.5:0.01:6;                      % Lateral distance (km)
z1=-1.8:0.01:0;                     % Depth (km)
[X1 Z1] = meshgrid(x1,z1);          % Grid of coordinates

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v=1.6 - 0.45.*Z1;                   % Depth dependent velocity model

% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the velocity model
figure(1)       
imagesc(x1,z1,v)
set(gca,'YDir','normal')
hold on
text(4.6,-1.5,'*','FontSize',15)     % Plotting the location of the source
text(4.5,-1.55,'Source')             % Labelling the source
colorbar
ylabel(colorbar,'Velocity (km/s)')
ylabel('Depth (km)')
xlabel('Lateral distance (km)')
title('1D Velocity model (dependent only on depth)')


%% 1.2. Ray Tracing in 1D velocity model

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global pdv_x pdv_z              % Global because the function needs them

X_source=[4.6,-1.5];            % Source, x and z coordinates
theta_0=linspace(0,2*pi,100);   % Take off angle (radians)
pdv_x = 0;                      % Partial derivative of v with respect to x
pdv_z = -0.45;                  % Partial derivative of v with respect to z
z=[-1.5:0.1:0];                 % z used for interpolation of traveltime
x=[3.5:0.1:6];                  % x used for interpolation of traveltime

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v_source=1.6-0.45*(-1.5);       % Velocity at source (km/s)
p_0 = [sin(theta_0)./v_source; cos(theta_0)./v_source];   %[px_0,pz_0]

% Plotting the velcity model
figure(2)
imagesc(x1,z1,v)
set(gca,'YDir','normal')
hold on

% Solving the ray equations and plotting the ray paths
for i = 1:length(theta_0);
[t raypath] = ode45('rayeq',[0 10],[X_source, p_0(1,i), p_0(2,i)]);


plot(raypath(:,1),raypath(:,2),'k') % Plotting x and z coordinates
% axis([3.5 6 -1.8 0])
hold on

if raypath(end,2)>-1.5;
    time_vec=interp2(raypath(:,1),raypath(:,2),t,x,z);
    T_time(i)=time_vec(end);

else raypath(end,2)<-1.5;
    T_time(i)=Inf;
end
end

% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
plot(X_source(1),X_source(2),'*g')
colorbar
ylabel(colorbar,'Velocity (km/s)')
ylabel('Depth (km)')
xlabel('Lateral distance (km)')
title('Ray tracing in a 1D velocity model')

figure(3)
plot((180/pi)*theta_0,T_time)
axis([0 360 0 15])
grid on
ylabel('Traveltime (s)')
xlabel('Take-off angle (degrees)')
title('Travel time curve')

%% 2. 2D velocity model

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_source=[4.6,-1.5]; % source, x and z coordinates
x=3.5:0.01:6;
z=-1.8:0.01:0;
[X Z] = meshgrid(x,z); 
global x0 z0 r0
x0 = 4.6;
z0 = -0.6;
r0 = 0.3;
r = sqrt((X-x0).^2 + (Z-z0).^2);

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v2 = 1.6 - 0.45*Z - 0.8.*exp((-r.^2)./r0^2);

% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4)
imagesc(x,z,v2)
set(gca, 'YDir', 'normal')           % Reversing the Y-axis
hold on
text(4.6,-1.5,'*','FontSize',15)     % Plotting the location of the source
text(4.5,-1.55,'Source')             % Labelling the source
colorbar
ylabel(colorbar,'Velocity (km/s)')
ylabel('Depth (km)')
xlabel('Lateral distance (km)')
title('2D Velocity model')
hold off

%% 2.2 Ray tracing in 2D velocity model

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_source=[4.6,-1.5];                % Source, x and z coordinates
theta_0=linspace(0,2*pi,100);       % Take off angle (radians)            

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_source=1.6-0.45*X_source(2)-0.8*exp(-((X_source(1)-x0)^2 + (X_source(2)-z0)^2)/r0^2);           % Velocity at source (km/s)
p_0 = [sin(theta_0)./v_source; cos(theta_0)./v_source];   %[px_0,pz_0]

figure(5)
imagesc(x,z,v2)
set(gca, 'YDir', 'normal')          % Reversing the Y-axis
hold on
for i = 1:length(theta_0);
[t2 raypath2] = ode45('rayeq2',[0 10],[X_source, p_0(1,i), p_0(2,i)]);

%figure(2)
plot(raypath2(:,1),raypath2(:,2),'k') % Plotting x and z coordinates
axis([3.5 6 -1.8 0])
hold on

if raypath2(end,2)>0;
    z2=-1.5:0.1:0;
    time_vec2=spline(raypath2(:,2),t2,z2);
    T_time2(i)=time_vec2(end);

else raypath2(end,2)<0;
    T_time2(i)=Inf;
end
end

% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
plot(X_source(1),X_source(2),'*g')
colorbar
ylabel(colorbar,'Velocity (km/s)')
ylabel('Depth (km)')
xlabel('Lateral distance (km)')
title('Ray tracing in a 2D velocity model')

figure(6)
plot((180/pi)*theta_0,T_time2)
axis([0 360 0 15])
grid on
ylabel('Traveltime (min)')
xlabel('Take-off angle (degrees)')
title('Travel time curve, 2D velocity model')



%% 3. 3D ray tracing in 1D velocity model

global pdv_x3 pdv_y3 pdv_z3

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x3 = 3.5:0.1:6;                    % Lateral distance (km)
y3 = x3;
z3 = -1.8:0.1:0;                   % Depth (km)
[X3 Y3 Z3] = meshgrid(x3,y3,z3);   % Grid of coordinates
X_source = [4.6,4.6,-1.5];         % Source, x, y and z coordinates
theta_0 = [0:pi/8:2*pi];           % Take off angle (radians), x/z plan
phi_0 = [0:pi/8:2*pi];             % Take off angle (radians), x/y plan

pdv_x3 = 0;                        % Partial derivative of v with respect to x
pdv_y3 = 0;                        % Partial derivative of v with respect to y
pdv_z3 = -0.45;                    % Partial derivative of v with respect to z

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V3=1.6 - 0.45.*Z3;                 % Depth dependent velocity model
v_source=1.6-0.45*(-1.5);          % Velocity at source (km/s)

figure(7)
for i = 1:length(phi_0)
    
for j = 1:length(theta_0)
    p3_0 = [(sin(theta_0(j)).*cos(phi_0(i)))./v_source; (sin(theta_0(j)).*sin(phi_0(i)))./v_source;...
    cos(theta_0(j))./v_source];
    
[t3 raypath3] = ode45('rayeq3',[0 10],[X_source, p3_0(1),p3_0(2),p3_0(3)]);

plot3(raypath3(:,1),raypath3(:,2),raypath3(:,3))
grid on
axis([3.5 6 3.5 6 -1.8 0])
hold on

% Checking if x/z plane and y/z plane is symmetrical and identical
% if phi_0(i) == 0
%     figure(8)
%     subplot(2,1,1)
%     plot(raypath3(:,1),raypath3(:,3))
%     axis([3.5 6 -1.8 0])
%     ylabel('Depth (km)')
%     xlabel('Lateral distance, x (km)')
%     hold on
% end
% 
% if phi_0(i) == pi/2
%     figure(8)
%     subplot(2,1,2)
%     plot(raypath3(:,2), raypath3(:,3))
%     axis([3.5 6 -1.8 0])
%     ylabel('Depth (km)')
%     xlabel('Lateral distance, y (km)')
%     hold on
% end


end
end

figure(7)
plot3(X_source(1),X_source(2),X_source(3),'*r')
zlabel('Depth (km)')
ylabel('Lateral distance, y (km)')
xlabel('Lateral distance, x (km)')
title('Ray tracing in 3D')

%% 3.2. 3D ray tracing in 2D velocity model

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = 4.6;
z0 = -0.6;
r0 = 0.3;
X_source = [4.6,4.6,-1.5];         % Source, x, y and z coordinates
theta_0 = [0:pi/8:2*pi];           % Take off angle (radians), x/z plan
phi_0 = [0:pi/8:2*pi];             % Take off angle (radians), x/y plan

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_source=1.6-0.45*X_source(3)-0.8*exp(-((X_source(1)-x0)^2 + (X_source(3)-z0)^2)/r0^2);           % Velocity at source (km/s)

figure(9)
for i = 1:length(phi_0)
    
for j = 1:length(theta_0)
    p3_0 = [(sin(theta_0(j)).*cos(phi_0(i)))./v_source; (sin(theta_0(j)).*sin(phi_0(i)))./v_source;...
    cos(theta_0(j))./v_source];
    
[t32 raypath32] = ode45('rayeq32',[0 10],[X_source(1),X_source(2),X_source(3), p3_0(1),p3_0(2),p3_0(3)]);


plot3(raypath32(:,1),raypath32(:,2),raypath32(:,3))
grid on
axis([3.5 6 3.5 6 -1.8 0])
hold on

% Checking x/z plane and y/z plane 
% if phi_0(i) == 0
%     figure(10)
%     subplot(2,1,1)
%     plot(raypath32(:,1),raypath32(:,3))
%     axis([3.5 6 -1.8 0])
%     ylabel('Depth (km)')
%     xlabel('Lateral distance, x (km)')
%     hold on
% end
% 
% if phi_0(i) == pi/2
%     figure(10)
%     subplot(2,1,2)
%     plot(raypath32(:,2), raypath32(:,3))
%     axis([3.5 6 -1.8 0])
%     ylabel('Depth (km)')
%     xlabel('Lateral distance, y (km)')
%     hold on
% end

end
end

figure(9)
plot3(X_source(1),X_source(2),X_source(3),'*r')
zlabel('Depth (km)')
ylabel('Lateral distance, y (km)')
xlabel('Lateral distance, x (km)')
title('Ray tracing in 3D')









