%% Exercise set 3 Gravity Anomalies

%% 1.

%a) Sketch

%b) & d)

%INPUT
G = 6.6742*10^-11;
drho = -200;

dz1=0.1;
dz2=0.3;
dz3=1;

zmin=800;
zmax=820;

z1=zmin:dz1:zmax;
z2=zmin:dz2:zmax;
z3=zmin:dz3:zmax;

nz1=length(z1);
nz2=length(z2);
nz3=length(z3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP for b)

%Integrand
g_int1=(z1.^-2) * drho*G;
g_int2=(z2.^-2) * drho*G;
g_int3=(z3.^-2) * drho*G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT b)
subplot(3,1,1)
plot(g_int1,z1)
axis ij
ylabel('Depth (m)')
xlabel('Integrand (m/s^2)')
title('Sampling 0.1')

subplot(3,1,2)
plot(g_int2,z2)
axis ij
ylabel('Depth (m)')
xlabel('Integrand (m/s^2)')
title('Sampling 0.3')

subplot(3,1,3)
plot(g_int3,z3)
axis ij
ylabel('Depth (m)')
xlabel('Integrand (m/s^2)')
title('Sampling 1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP for d)

dgz1 = dz1*sum(g_int1)-dz1*(0.5)*(g_int1(1) + g_int1(nz1)); %delta gz in m/s^2 
dgz2 = dz2*sum(g_int2)-dz2*(0.5)*(g_int2(1) + g_int2(nz2)); 
dgz3 = dz3*sum(g_int3)-dz3*(0.5)*(g_int3(1) + g_int3(nz3));  

%Analytical 
dgz_an=-drho*G*((1/zmax)-(1/zmin));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT for d)
'Gravity response for samplings 0.1,0.3 and 1, in mGal'
mdgz1=dgz1*10^5 %In mgal (cm->m:*100,gal->mgal:*1000
mdgz2=dgz2*10^5
mdgz3=dgz3*10^5

%Analytical
mdgz_an=dgz_an*10^5

%c) Gal: unit for acceleration. 1 gal = 1 centimeter per second squared (1
%cm/s^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2.

%INPUT
drho=-200;
Z = 810;
G = 6.6742*10^-11;

%Gravitmeter position
x0=0;
y0=0;
z0=0;

%Defining intervals
dx= 4000/200;       % = 20
dx2= 4000/100;      % = 40
dx3=4000/40;        % = 100

dy=4000/200;
dy2= 4000/100;
dy3=4000/40;

dz=20;

%Sampling interval=20:
x=[-2000:dx:2000];
y=[-2000:dy:2000];
[X Y]=meshgrid(x,y);

%Sampling interval=40:
x2=[-2000:dx2:2000];
y2=[-2000:dy2:2000];
[X2 Y2]=meshgrid(x2,y2);

%Sampling interval=100:
x3=[-2000:dx3:2000];
y3=[-2000:dy3:2000];
[X3 Y3]=meshgrid(x3,y3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP

term1 = Z-z0;           %Nominator

%Sampling interval = 20
term2_20= sqrt(((X-x0).^2 + (Y-y0).^2 + (Z-z0)^2).^3);  %Denominator
int_20 = (term1./term2_20)*G*drho; % Integrand in m/s^2
int_20_mgal = int_20*10^5;         % In mgal

%Sampling interval = 40
term2_40= sqrt(((X2-x0).^2 + (Y2-y0).^2 + (Z-z0)^2).^3);
int_40 = (term1./term2_40)*G*drho;
int_40_mgal = int_40*10^5;

%Sampling interval = 100
term2_100= sqrt(((X3-x0).^2 + (Y3-y0).^2 + (Z-z0)^2).^3);
int_100 = (term1./term2_100)*G*drho;
int_100_mgal = int_100*10^5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(21)
subplot(2,2,1)
mesh(X,Y,int_20_mgal)
title('Sampling interval = 20m')
xlabel('X coordinates (m)'); ylabel('Y coordinates (m)'); zlabel('delta g (mgal)')

subplot(2,2,2)
mesh(X2,Y2,int_40_mgal)
title('Sampling interval = 40m')
xlabel('X coordinates (m)'); ylabel('Y coordinates (m)'); zlabel('delta g (mgal)')

subplot(2,2,3)
mesh(X3,Y3,int_100_mgal)
title('Sampling interval = 100m')
xlabel('X coordinates (m)'); ylabel('Y coordinates (m)'); zlabel('delta g (mgal)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2b.

%INPUT

%Defining lengths
nx1=length(x);
nx2=length(x2);
nx3=length(x3);

ny1=length(y);
ny2=length(y2);
ny3=length(y3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP

sum_int_20=sum(sum(int_20_mgal(1:nx1-1,1:ny1-1))); %sum(sum(f)) - summing matrix
sum_int_40=sum(sum(int_40_mgal(1:nx2-1,1:ny2-1)));
sum_int_100=sum(sum(int_100_mgal(1:nx3-1,1:ny3-1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT

'Gravity respone in ugal, for sampling intervall=20'
dx*dy*dz*sum_int_20*10^3   %*10^3 to get ugal (microgal)
                            %No need for G and drho because it's already in
                            %int_100mgal

'Gravity respone in ugal, for sampling intervall=40'
dx2*dy2*dz*sum_int_40*10^3

'Gravity respone in ugal, for sampling intervall=100'
dx3*dy3*dz*sum_int_100*10^3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2d)

%INPUT

dz2=15;     %thickness of gas layer after 3 years (20m -5m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT

'Gravity respone in ugal after 3 years, for sampling intervall=20'
dx*dy*dz2*sum_int_20*10^3   %*10^3 to get ugal (microgal)

'Gravity respone in ugal after 3 years, for sampling intervall=40'
dx2*dy2*dz2*sum_int_40*10^3

'Gravity respone in ugal after 3 years, for sampling intervall=100'
dx3*dy3*dz2*sum_int_100*10^3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2e)

%From 2b we know that 20 m of gas gives a gravity response of -110.36 ugal.
%The thickness that gives a signal of 5 ugal is therefore calculated by:
'Min thickness change (m) ='
(20/-110.36)*-5 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2f)

%INPUT

ndx=3000/200;          %new interval
xn=[-1000:ndx:2000];
[Xn Y] = meshgrid(xn,y);
n_xn=length(xn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP

%Integrand
term1 = Z-z0;           %Nominator

term2_n= sqrt(((Xn-x0).^2 + (Y-y0).^2 + (Z-z0)^2).^3);  %Denominator
int_n = (term1./term2_n)*G*drho; % Integrand in m/s^2
int_n_mgal = int_n*10^5;         % In mgal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
'Gravity response after 3 years of horizontal influx in ugal'
sum_int_n=(sum(sum(int_n_mgal(1:n_xn-1,1:ny1-1))));
ndx*dy*dz*sum_int_n*10^3 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3.

%INPUT

%Gravimeter array positions
x0_3=[-3000:500:3000];
y0=0;
z0=0;

%Thickness of anomaly
dz1=20;
dz2=15;     %At time 2
    
N=length(x0_3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMP

term1 = Z-z0;           %Nominator

for n=1:N
    term2_20= sqrt(((X-x0_3(n)).^2 + (Y-y0).^2 + (Z-z0)^2).^3);  
    int_20 = (term1./term2_20)*G*drho;
    int_20_mgal = int_20*10^5;
    %for time 1, dz1=20
    grav_resp1(n)=dx*dy*dz1*(sum(sum(int_20_mgal(1:nx1-1,1:ny1-1))))*10^3;
    %for time 2, dz2=15
    grav_resp2(n)=dx*dy*dz2*(sum(sum(int_20_mgal(1:nx1-1,1:ny1-1))))*10^3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT

figure(31)
plot(x0_3,grav_resp1,x0_3,grav_resp2)
title('Gravity response in array')
xlabel('Gravimeter position (x-coordinates in m)')
ylabel('Gravity response (ugal)')
legend('Gravity response at t=1','Gravity response at t=2')


























