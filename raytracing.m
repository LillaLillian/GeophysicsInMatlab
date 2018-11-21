%% 1. Nummerical solution to 1D model

% Horizontal dimensions are 3.5km to 6 km and the vertical dimensions are
% 0km to -1.8km.

% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global vel grad_x grad_z Z1 X1      % Defining global variables for function
% global xstart xend zstart zend    % Til stoppskript

% xstart = 3.5;                     % Til stoppskript
% xend=6;
% zstart=0;
% zend=-1.8;

x=-10:0.01:10;                       % Lateral distance (km)
z=-1.8:0.01:0;                      % Depth (km)
dx=0.01;                            % x-increment needed for gradient
dz=dx;                              % z-increment needed for gradient
X_source=[4.6,-1.5];                % Source, x and z coordinates
[X1 Z1] = meshgrid(x,z);            % Grid of coordinates
theta_0=[0:pi/20:2*pi-pi/20];       % Take off angle (radians)  
zs=0;                               % z_surface, used for interpolation

% Computations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vel=1.6 - 0.45.*Z1;                     % Creates a grid of velocity
[grad_x grad_z] = gradient(vel,dx,dz);  % Gradient needed for raytracing
v_source=1.6-0.45*X_source(2);          % Velocity at source (km/s)
p_0 = [sin(theta_0)./v_source; cos(theta_0)./v_source];   %[px_0,pz_0]

% options = odeset('Events',@stoppskript);  % Til stoppskript

% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
imagesc(x,z,vel)            % Velocity model
axis([3.5 6 -1.8 0.2])
set(gca,'YDir','normal')
hold on

for i = 1:length(theta_0);
[tn raypathn] = ode45('rayeqn',[0 10],[X_source, p_0(1,i), p_0(2,i)]);

X = raypathn(:,1);
Z = raypathn(:,2);

plot(X,Z,'k-*')
%axis([3.5 6 -2 2])
ylabel('Z, Depth (km)')
xlabel('X-coordinates (km)')
hold on

if Z(end)>0; %Interpolating x and t to the surface, for the rays reaching the surface
    n=find(Z<0);
    m=find(Z>0);
    x1=X(n(end));
    x2=X(m(1));
    z1=Z(n(end));
    z2=Z(m(1));
    
    % Interpolating for x at surface
    a=(z2-z1)/(x2-x1);  % delta y/ delta x, stigningstall
    b= z1-a*x1;
    xs=-b/a;
     if x2==x1
        xs=x2;
    end
    plot(xs,zs,'*y')
    
    if xs<6 & xs>3.5
    % Interpolating for traveltime at surface
    t1=tn(n(end));
    t2=tn(m(1));
    l1=sqrt((xs-x1)^2+(zs-z1)^2);
    l2=sqrt((x2-x1)^2+(z2-z1)^2);
    ts(i)=t1+(l1/l2)*(t2-t1);
    end
    
end

% ODE23 Kommandoen
[tn2 raypathn2] = ode23('rayeqn',[0 10], [X_source, p_0(1,i), p_0(2,i)]);
hold on
plot(raypathn2(:,1),raypathn2(:,2),'m-*')

X2 = raypathn2(:,1);
Z2 = raypathn2(:,2);

if Z2(end)>0; %Interpolating x and t to the surface, for the rays reaching the surface
    n2=find(Z2<0);
    m2=find(Z2>0);
    x1_2=X2(n2(end));
    x2_2=X2(m2(1));
    z1_2=Z2(n2(end));
    z2_2=Z2(m2(1));
    
    % Interpolating for x at surface
    a2=(z2_2-z1_2)/(x2_2-x1_2);  % delta y/ delta x, stigningstall
    b2= z1_2-a2*x1_2;
    xs_2=-b2/a2;
    if x2_2==x1_2
        xs_2=x2_2;
    end
    plot(xs_2,zs,'xr')
    
    if xs_2<6 & xs_2>3.5
    % Interpolating for traveltime at surface
    t1_2=tn2(n2(end));
    t2_2=tn2(m2(1));
    l1_2=sqrt((xs_2-x1_2)^2+(zs-z1_2)^2);
    l2_2=sqrt((x2_2-x1_2)^2+(z2-z1_2)^2);
    ts_2(i)=t1_2+(l1_2/l2_2)*(t2_2-t1_2);
    end
end

end
plot(X_source(1),X_source(2),'*r')
legend('ode45','surface ode45','ode23', 'surface ode23')

hold off

figure(2)
plot(theta_0*(180/pi),ts,'*-')
title('Traveltime vs take-off angle')
hold on
plot(theta_0*(180/pi),ts_2,'x-r')
legend('ode45','ode23')


%% 2. Smooth marmousi model

% Setting up the marmousi model
m=load('marmousi.txt');
mar=reshape(m,122,384);
xm=0:24:9192;
zm=0:-24:-2916;

figure(8)
imagesc(xm,zm,mar)
set(gca,'YDir','normal')

% Using 1 profile for the whole slice
smar=zeros(size(mar));
for i = 1:122
    smar(i,:) = mar(i,80);
end

figure(9)
imagesc(xm,zm,smar)
set(gca,'YDir','normal')
hold on

% Doing the ray tracing in simplified marmousi (1 profile)
x_source= [4500,-2000];
vel=smar;
[X1 Z1] = meshgrid(xm,zm);
theta_0=linspace(0,2*pi,10);       % Take off angle (radians) 
vn_source=interp2(X1,Z1,vel,x_source(1),x_source(2),'linear',1.6);

p_0 = [sin(theta_0)./vn_source; cos(theta_0)./vn_source];   %[px_0,pz_0]

for i =1:1  % for one ray
%for i = 1:length(theta_0)
   [tn raypathn] = ode45('rayeqn',[0 10],[x_source,p_0(1,i),p_0(2,i)]);
   plot(raypathn(:,1),raypathn(:,2),'k')
   hold on
end

