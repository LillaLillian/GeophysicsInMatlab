%Exercise set 2 Geov219
% Nummerical derivation and interpolation

%% 1. McMahon Chapter 2

%1.
A=[-1 7 3 2];
a=sum(A.*A);        % sum of array multiplication
MagA = sqrt(a)      %Magnitude of A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2.
A2=[-1+i 7i 3 -2-2i];   % Complex #
At=conj(A2);            % Complex conjugate of A2
a2=sum(At.*A2);
MagA2 = sqrt(a2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.
v1=[1 2 3];             %Row vector
v2=[1; 2; 3]            %Column vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4.
clear all

A=[1;2;3];
B=[4;5;6];
Ans=A.*B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%5. Making an identity matrix 5x5
eye(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%6. 
clear all

A=[8 7 11; 6 5 -1; 0 2 -8];
B=[2 1 2; -1 6 4; 2 2 2];

'Array product of A and B ='
A.*B

'Matrix product of A and B ='
A*B

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%7.
clear all

A=[1 2 3; 4 5 6; 7 8 9];
B=[A(3,:);A(3,:);A(2,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%8.
clear all
a=[1 2 3; -4 1 2; 9 -8 -1]; 
b=[12;13;-1];               
c=a\b;      %solution
det(a);     %determinant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%9.
clear all
a=[1 2 3; -4 1 2; 0 9 -8];
b=[1;2;3];
c=a\b;         %Solution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%10.
clear all
a=[1 7 -9; 2 -1 4; 1 1 -7];
b=[12;16;16];
[L, U] = lu(a);
c=U\(L\b);      %Solution, using LU decomposition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. McMahon Chapter 3

%1.
%INPUT
x=[0:0.1:1];

%COMP
y=tan(x);

%OUTPUT
figure(21)
plot(x,y)
ylabel('tan(x) & sin(x)')
xlabel('x = [0:0.1:1]')
title('Plot of tan(x) and sin(x) vs x (McMahon 1.)')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%2.
%COMP
z=sin(x);

%OUTPUT
plot(x,z,'r')
legend('tan(x)','sin(x)')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%3.
x=[-pi:0.2:pi];
x2=linspace(-pi,pi,100);
x3=linspace(-pi,pi,50);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4.
clear all
x=[-3:0.1:2];
y=[-5:0.1:5];
[xx yy]= meshgrid(x,y);

x2=[-5:0.2:5];
y2=[-5:0.2:5];

[xx2 yy2]=meshgrid(x2,y2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5.
clear all
t=[0:pi/10:2*pi];
x=exp(-t).*cos(t);
y=exp(-t).*sin(t);
z=t;

figure(22)
plot3(x,y,z)
grid on
title('3D plot (McMahon 5.)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3. 1D interpolation
clear all 

%a)
%INPUT
x=rand(1,11)*10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
z=rand(1,11)+(rand(1,11)-1);
y=x+z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(31)
plot(x,y,'*')
xlabel('x-values')
ylabel('y-values')
title('Plot of random x and y values (3a)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3b)  Interp1 interpolation

%INPUT
a=rand(1,101)*10;
b=interp1(x,y,a);       %interp1(X,Y,X1); interpolates to find Y1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(32)
plot(a,b,'*')
xlabel('x-values')
ylabel('y-values')
title('Y values interpolated (3b)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3c)  Spline interpolation

%INPUT
b2=spline(x,y,a);       %interp1(X,Y,X1); interpolates to find Y1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(33)
plot(a,b2,'*')
xlabel('x-values')
ylabel('y-values')
title('Y values interpolated (3c)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3d)
%Assuming that the spline interpolated (b2) results are correct:

%COMP
Abs_err = abs(b2-b);
Rel_err= (abs(b2-b))./(abs(b2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3e) 

%OUTPUT
figure(34)
subplot(2,2,1)
plot(a,b,'*')
xlabel('x-values')
ylabel('y-values')
title('Y values interpolated with interp1 (3e)')

subplot(2,2,2)
plot(a,b2,'*')
xlabel('x-values')
ylabel('y-values')
title('Y values interpolated with spline interpolation')

subplot(2,2,3)
plot(a,Abs_err,'*')
xlabel('x-values')
ylabel('y-values')
title('Absolute error')

subplot(2,2,4)
plot(a,Rel_err,'*')
xlabel('x-values')
ylabel('y-values')
title('Relative error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4. 1D Derivatives
%a)

%INPUT
x=linspace(0,pi,11);
f=sin(x);
N=length(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP 4a)

for n=1:N-1
    num_der(n)= (f(n+1)-f(n))/(x(n+1)-x(n));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMP 4b)
%f(x) = sin(x)
%f'(x) = cos(x)

f_der=cos(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT 4a & 4b

figure(41)
plot(x,f,x(1:10),num_der,x,f_der)
legend('f(x)=sin(x)','Numerical derivative of f(x)','Analytical derivative of f(x)')
xlabel('x-values')
title('f(x) and its numerical and analytical derivative for 11 equidistant points (4a/b)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4c)

%INPUT
x2=linspace(0,pi,101);
f2=sin(x2);
N2=length(x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
%Numerical derivative
for n=1:N2-1
    num_der2(n)= (f2(n+1)-f2(n))/(x2(n+1)-x2(n));
end

% Analytical derivative
%f2(x) = sin(x2)
%f2'(x) = cos(x2)

f_der2=cos(x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT

figure(42)
plot(x2,f2,x2(1:100),num_der2,x2,f_der2)
legend('f(x)=sin(x)','Numerical derivative of f(x)','Analytical derivative of f(x)')
xlabel('x-values')
title('f(x) and its numerical and analytical derivative for 101 equidistant points (4c)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4d) Assuming the analytical derivative is correct:

%COMP
abs_err_4c=abs(f_der2(1:100)-num_der2);
rel_err_4c=abs((f_der2(1:100)-num_der2)./(f_der2(1:100)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(43)
subplot(2,2,1)
plot(x2(1:100),num_der2)
title('Numerical derivative for f(x) at 101 equidistant points (4d)')
xlabel('x-values')
ylabel('Num derivative of f(x)')

subplot(2,2,2)
plot(x2,f_der2)
title('Analytical derivative for f(x) at 101 equidistant points')
xlabel('x-values')
ylabel('Analytical derivative of f(x)')

subplot(2,2,3)
plot(x2(1:100),abs_err_4c)
xlabel('x-values')
ylabel('y-values')
title('Absolute error')

subplot(2,2,4)
plot(x2(1:100),rel_err_4c)
xlabel('x-values')
ylabel('y-values')
title('Relative error')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4e)

%INPUT
x=linspace(0,pi,101);
h=pi/101;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
f=sin(x);
f_der=diff(f)/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(44)
plot(x(1:100),f_der,x,f,x,f_der2) %x(1:100) = forward differentiation
title('Matlabs numerical derivative (diff) of sin(x)')
xlabel('x = 0:pi')
ylabel('Derivative of sin(x)')
legend('Numerical derivative of f(x), using diff','f(x)','Analytical derivative of f(x)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matlabs diff command gives a very similar result as the numerical 
% derivative with 101 points. The diff command and the nummerical 
% derivative with 101 points is both quite similar to the analytical
% derivative, and give a good approximation. The derivative
% using only 11 points is not as precise as the two others and gives a 
% poorer approximation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6. 2D interpolation

% 6a) On paper

% 6b) On paper

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6c) 

%INPUT
x1=4;
x2=20;
y1=5;
y2=2;
x=5;
y=10;
x0=2;
y0=5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
%Calculating the 4 points to interpolate between
v_11= x1 + exp(-(((x1-x0).^2)+(y1-y0).^2)); % 1. point v_11(x1,y1)
v_21= x2 + exp(-(((x2-x0).^2)+(y1-y0).^2)); % 2. point v_21(x2,y1)
v_12= x1 + exp(-(((x1-x0).^2)+(y2-y0).^2)); % 3. point v_12(x1,y2)
v_22= x2 + exp(-(((x2-x0).^2)+(y2-y0).^2)); % 4. point v_22(x2,y2)

%Calculating the interpolation, part by part
term1= 1/((y2-y1)*(x2-x1));

term_v11=v_11.*(x2-x).*(y2-y);
term_v21=v_21.*(x-x1).*(y2-y);
term_v12=v_12.*(x2-x).*(y-y1);
term_v22=v_22.*(x-x1).*(y-y1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
%Interpolated velocity
v_P= term1*(term_v11 + term_v21 + term_v12 + term_v22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6d)

%INPUT
load marmousi.txt

%Defining x and z coordinates in meter
dx=24;
x=(0:dx:9192); %m
z=(0:dx:2904); 

% x and z for interpolated model
dx2=12;
x2=(0:dx2:9192); 
z2=(0:dx2:2904); 

%Defining meshgrids
[X,Z]=meshgrid(x,z);        %Original model
[X2,Z2]=meshgrid(x2,z2);    %Interpolated model

%Reshaping the marmousi model
mar=reshape(marmousi,122,384); 

%Location of the chosen 1D velocity profile
loc=mar(1:122,100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
v_int=interp2(X,Z,mar,X2,Z2);   %Interpolating marmousi
loc_int=v_int(:,100);           %Interpolated marmousi values for chosen location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(61)
subplot(1,3,1)
stairs(loc, z)
axis ij
title('Original 1D velocity profile')
xlabel('Velocity (m)')
ylabel('Depth (m)')

subplot(1,3,2)
stairs(loc_int,z2)
axis ij
title('Interpolated 1D velocity profile, using command interp2')
xlabel('Velocity (m)')
ylabel('Depth (m)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6e)

%COMP
v_spline=interp2(X,Z,mar,X2,Z2,'spline'); %interpolation marmousi using the spline method
loc_spline=v_spline(:,100);               %spline interpolated marmousi values for chosen location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
subplot(1,3,3)
stairs(loc_spline,z2)
axis ij
title('Interpolated 1D velocity profile, using command interp2 with spline method ')
xlabel('Velocity (m)')
ylabel('Depth (m)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 6f)

%Assuming the spline interpolation is the correct one

%COMP
abs_err_6f=abs(loc_spline-loc_int);
rel_err_6f=abs(loc_spline-loc_int)./(abs(loc_spline))*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(62)
subplot(2,2,1)
plot(abs_err_6f,z2)
axis ij
title('Absolute error')
xlabel('Absolute error')
ylabel('Depth (m)')

subplot(2,2,2)
stairs(loc_spline,z2)
axis ij
title('Interpolated 1D velocity profile, using command interp2 with spline method ')
xlabel('Velocity (m)')
ylabel('Depth (m)')

subplot(2,2,3)
plot(rel_err_6f,z2)
axis ij
title('Relative error')
xlabel('Relative error')
ylabel('Depth (m)')

subplot(2,2,4)
stairs(loc_int,z2)
axis ij
title('Interpolated 1D velocity profile, using command interp2')
xlabel('Velocity (m)')
ylabel('Depth (m)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6g) It seems that both the absolute error and the relative error is
% greater where the velocity is greater

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 7. 2D Derivatives
clear all
close all

% 7a)
%INPUT
%Defining x and y values
x=[-1:0.1:1];
y=[-2:0.1:2];
x0=0;
y0=0;

%Defining the meshgrid
[X,Y]=meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
v=X+exp(-(((X-x0).^2)+(Y-y0).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(71);
subplot(2,1,1);
imagesc(x,y,v);
xlabel('x-coordinates')
ylabel('y-coordinates')
colorbar
ylabel(colorbar,'Velocity')

subplot(2,1,2);
surf(X,Y,v)
xlabel('x-coordinates')
ylabel('y-coordinates')
zlabel('Velocity')

figure(72)          %surf and imagesc will not plot nicely together in subplot
imagesc(x,y,v)
colorbar
xlabel('x-coordinates')
ylabel('y-coordinates')
ylabel(colorbar,'Velocity')

figure(73)
surf(X,Y,v)
xlabel('x-coordinates')
ylabel('y-coordinates')
zlabel('Velocity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 7b

%INPUT
%Defining x and y values
x=[-1:0.1:1];
y=[-2:0.1:2];
x0=0;
y0=0;

%Defining the meshgrid
[X,Y]=meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
%Because xo=0 and y0=0, v is given as:
%v= X + exp(-(X.^2+Y.^2));

%Analytical derivative of v:
v_der_x= 1 - 2.*X.*(exp(-(X.^2+Y.^2)));     %With respect to X
v_der_y= - 2.*Y.*(exp(-(X.^2+Y.^2)));       %With respect to Y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
figure(74)
subplot(2,2,1)
surf(X,Y,v_der_x)
title('Analytical derivative, with respect to X')
xlabel('X coordinates')
ylabel('Y coordinates')
zlabel('Velocity')

subplot(2,2,2)
surf(X,Y,v_der_y)
title('Analytical derivative, with respect to Y')
xlabel('X coordinates')
ylabel('Y coordinates')
zlabel('Velocity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%7c)

%INPUT
syms X Y      %Defining X and Y as symbols to do partial derivatives with diff
v=inline('X + exp(-((X).^2+(Y).^2))','X','Y');  %Defines v as a function of 'X' and 'Y'
der_vx=inline(diff(v(X,Y),X),'X','Y');          %Using diff to calculate derivative 
der_vy=inline(diff(v(X,Y),Y),'X','Y');          %with respect to X and then Y

[X,Y]=meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
subplot(2,2,3)
surf(X,Y,der_vx(X,Y))       %Plotting the derivative with respect to X for all values of X and Y
title('Numerical derivative, with respect to X')
xlabel('X coordinates')
ylabel('Y coordinates')
zlabel('Velocity')

subplot(2,2,4)
surf(X,Y,der_vy(X,Y))       %Plotting all the derivative with respect to Y for all values of X and Y
title('Numerical derivative, with respect to Y')
xlabel('X coordinates')
ylabel('Y coordinates')
zlabel('Velocity')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 7d)

%COMP
syms X Y
der2_vx=inline(diff(der_vx(X,Y),X),'X','Y')         %Using diff to calculate 2. derivative 
                                                    %with respect to X
                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7e)

%COMP
der2_vy=inline(diff(der_vy(X,Y),Y),'X','Y')        %Using diff to calculate 2. derivative 
                                                   %with respect to Y 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 7f)

%INPUT
x = [-1:0.1:1];
y = [-2:0.1:2];
x0 = 0;
y0 = 0;

[X,Y] = meshgrid(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COMP
v = X + (exp(-((X-x0).^2 + (Y-y0).^2)));
L = 4*del2(v,0.1,0.1);
Laplacian = sum(sum(L));

L_de=der2_vx(X,Y)+der2_vy(X,Y);
Laplacian_de=sum(sum(L_de));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OUTPUT
'del2 command yields:'
Laplacian

'Using exercise 7d and 7e yields:'
Laplacian_de




