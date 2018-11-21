% Geov219 1.st Exercise Set

%% 1. Basic Plotting
%a)

%input
x=linspace(0,pi,101); %start 0, end pi, 101 numbers

%Computation
y=sin(x);

%Output
figure(11)
plot(x,y,'-x','Linewidth',2)
xlabel('x-value (0-pi)')
ylabel('sin(x)')
title('The sin function at 101 equidistant points')
hold on


%b)
plot(x,x,'o')
%c)

%Comp
y1=(x-((x.^3)/factorial(3)));
y2=y1 +((x.^5)/factorial(5));
y3=y2-((x.^7)/factorial(7));

%Output
figure(13)
plot(x,y,x,y1,'g',x,y2,'r',x,y3,'c')
title('Taylor series')
xlabel('x')
ylabel('y')
legend('y=sin(x)','y1 = x-(x^3/3!)', 'y2=y1+(x^5/5!)', 'y3=y2-(x^7/7!)')

%d)
% a&b) While x is constantly increasing, the sine function will still form
% waves
%
% c) The taylor series is an approximation of sine of x, and the more terms
% used the closer the plot gets to the originall sin(x)

%% 2. 1D velocity models

%a)
%input
z=linspace(0,1000,21);
a=1500;
b=0.1;

%Computation
v=a*z.^b;

%Output
figure(21)
plot(v,z)
axis ij
xlabel('Velocity (m/s)')
ylabel('Depth (m)')
title('1D velocity model')

%b)
figure(22)
stairs(v,z)
axis ij
xlabel('Velocity (m/s)')
ylabel('Depth (m)')
title('1D velocity model, stairs function')

%% 3. Interpolation

%a) 
% Known values: f_i, f_i+1, x_i, x_i+1, x
% Need to find: f
% Using linear interpolation the function is assumed to be of the form
% f=ax+b, where we have to find a and b. This is done by writing out  
% the known equations:
% 
%   f_i = ax_i + b
%   f_i+x = ax_i+1 + b
%
% Substracting f_i from f_i+x gives us an expression for a:
% 
% a= (f_i+1 - f_i)/(x_i+1 - x_i)
%
% When a is found, b can be found from one of the previous forumlas:
%
%   f_i = ax_i + b
%
%   b = f_i - ax_i
%
% When a and b is known, f is found from the formula:
%
%   f = ax + b

%b)
% Known values: f_i, f_i+1, f'_i, f'_i+1, x_i, x_i+1, x
% Need to find: f
% In this case we assume the function to be of the form: 
% f = ax^3 + bx^2 + cx + d, in which case the derivative is given as:
% f' = 3ax^2 + 2bx + c
% Using these equations for f'_i, f'_i+1, x_i, x_i+1 gives us 4 equations
% and 4 unknowns (a,b,c,d). Arranging this into matrixes and solving with
% linear algebra gives us a,b,c and d.
% 
% | f_i    |      | x_i^3     x_i^2    x_i     1 |       | a |
% | f_i+1  |   =  | x_i+1^3   x_i+1^2  x _i+1  1 |   *   | b |
% | f'_i   |      | 3x_i^2    2x_i     1       0 |       | c |
% | f'_i+1 |      | 3x_i+1^2  2x_i+1   1       0 |       | d |
% 
% When a,b,c and d is known, f is found from the equation:
% 
% f = ax^3 + bx^2 + cx + d
% 

%%

%c)

%Input
z_int = 25:50:975;
a = 1500;
b = 0.1;
z = linspace(0,1000,21);
zi = 275;

%computation
v=a*z.^b;

%output
'Velocity at 275m = '

vi=interp1(z,v,zi)


%%
clear all
close all

%d)
z=linspace(0,1000,21);
z_int=25:50:975;
a=1500;
b=0.1;

for n =1:1:20;
    z1=z(n);
    z2=z(n+1);
    v1=a*z1.^b;
    v2=a*z2.^b;
    
    a_n=(v2-v1)/(z2-z1);
    b_n= v1 - a_n*z1;
    v(n)=a_n*z_int(n)+b_n;
end

plot(v,z_int,'*')
title('Interpolation of velocity model')
xlabel('Interpolated velocity (m/s)')
ylabel('Depth (m)')
axis ij

%% 4. Seismic Exploration and 2D plotting

%a)

%Input
z=linspace(-2,0,100); %100 points, from -2 to 0 (Depth km)

%Comp
v=1.6 - 0.45.*z; %simplified velocity model

%Output
figure(41)
plot(v,z);
xlabel('Velocity (km/s)')
ylabel('Depth (km)')
title('Simplified velocity model Valhall')

%%
%b & c)

%Input
ro=0.3; %km
zo=-0.6; %km
xo=4.6; %km


[xx zz] = meshgrid(3.5:0.02:5.5,-1.5:0.01:0);

%Comp
r= sqrt((xx-xo).^2 + (zz-zo).^2) + ro;
v2=1.6-0.45*zz-0.8*exp(-((r-ro).^2)/ro^2); %More complicated vel model

%Output
figure(42)
contourf(xx,zz,v2)
colorbar
xlabel('x-coordinates')
ylabel('Depth (km)')
ylabel(colorbar,'Velocity (km/s)')
title('Complicated velocity model, Valhall (contourf)')


figure(43)
imagesc(xx(1,:),0:0.01:1.5,flipud(v2)) %plot only first row of x (1,:), and first 
                                    %column of z (:,1)
xlabel('x-coordinates')
ylabel('Depth (km)')
title('Complicated velocity model, Valhall (imagesc)')
colorbar
ylabel(colorbar,'Velocity (km/s)')

%% 5. Plotting a 2D Seismic Veloocity Model

clear all 
close all

%a)

%Input
load marmousi.txt

dx=24;
x=(0:dx:9192); %m
z=(0:dx:2904); %m

mar=reshape(marmousi,122,384);      %reshaping the marmousi model to a matrix

%Output
figure(51)
imagesc(x,z,mar)
colorbar
ylabel('Depth (m)')
xlabel('x-coordinates (m)')
title('Marmousi model, imagesc')
ylabel(colorbar,'Velocity (m/s)')    %label the colorbar

%%
%b

%Input
[xx zz] = meshgrid(x,z);    % contour needs matrixes, (meshgrid)

%Output
figure(52)
contour(xx,zz,mar)
axis ij
colorbar
ylabel('Depth (m)')
xlabel('x-coordinates (m)')
title('Marmousi model, contour')
ylabel(colorbar,'Velocity (m/s)')

figure(53)
contourf(xx,zz,mar)
axis ij
colorbar
ylabel('Depth (m)')
xlabel('x-coordinates (m)')
title('Marmousi model, contourf')
ylabel(colorbar,'Velocity (m/s)')

%%
%c)

Minimum_vel = min(marmousi)             %Min velocity of Marmousi model
Maximum_vel = max(marmousi)             %Max velocity of Marmousi model
Average_vel = mean(marmousi)            %Average velocity


%%
%d
%Input
mar=reshape(marmousi,122,384); 
loc1=mar(1:122,100);
loc2=mar(1:122,200);
loc3=mar(1:122,300);

%Output

%Profile 1

figure(54)

%1D velocity
subplot(3,2,1)
stairs(loc1,z)
axis ij
ylabel('Depth (m)')
xlabel('Velocity (m/s)')
title('1D depth-velocity profile, loc1')

mar=reshape(marmousi,122,384); 
mar(:,100)=[1000];   %set all numbers in row 100 to 1000 to make a line in the profile
                     
subplot(3,2,2)          % Plotting the profile with a line marking the
contourf(xx,zz,mar)     % location of the 1D depth-velocity model                    
axis ij
colorbar
ylabel('Depth (m)')
xlabel('x-coordinates (m)')
title('Marmousi model, loc1')
ylabel(colorbar,'Velocity (m/s)')

%Profile #2
subplot(3,2,3)
stairs(loc2,z)
axis ij
ylabel('Depth (m)')
xlabel('Velocity (m/s)')
title('1D depth-velocity profile, loc2')


mar=reshape(marmousi,122,384); 
mar(:,200)=[1000]; 

subplot(3,2,4)
contourf(xx,zz,mar)
axis ij
colorbar
ylabel('Depth (m)')
xlabel('x-coordinates (m)')
title('Marmousi model, loc2')
ylabel(colorbar,'Velocity (m/s)')

%Profile #3
subplot(3,2,5)
stairs(loc3,z)
axis ij
ylabel('Depth (m)')
xlabel('Velocity (m/s)')
title('1D depth-velocity profile, loc3')

mar=reshape(marmousi,122,384); 
mar(:,300)=[1000]; 

subplot(3,2,6)
contourf(xx,zz,mar)
axis ij
colorbar
ylabel('Depth (m)')
xlabel('x-coordinates (m)')
title('Marmousi model, loc3')
ylabel(colorbar,'Velocity (m/s)')









