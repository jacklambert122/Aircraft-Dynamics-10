%% Author: Jack Lambert
% ASEN 3128
% Homework 10
% Purpose: Find the dimensional derivatives given the none dimensional
% derivatives from pg.187 of Etkin for the the flight conditions of a
% Boeing 747, given on by Case 2 in table E.1. Also finds the changes in the
% y-component bofy force, roll moment, and yaw moment with changes in
% y-component of velocity, x-comp of angular velocity, and z-comp of 
% angular velocity 
% Date Modified: 4/20/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A] = Amat()
%% Airplane Parameters
% Nondimensional Lateral Derivatives
% Table 6.6 -
Cy = [-.8771, 0, 0];
Cl = [-.2797, -.3295, .304];
Cn = [.1946, -.04073, -.2737];

% Table E.1 B747 Case 2
Alt = 20000*(0.3048); % Altitude [ft] -> [m]
[~, ~, ~, rho] = atmosisa(Alt); % Standard Atmosphere Properties at Alt.
W = 6.366*10^5*4.44822; % Weight [lb]->[N]
Ix_PA = 1.82e7*1.35581795; % Moment of Interia x-PA [slug ft^2]-> [kg m^2]
Iy_PA = 3.31e7*1.35581795; % Moment of Interia y-PA [slug ft^2]-> [kg m^2]
Iz_PA = 4.97e7*1.35581795; % Moment of Interia z-PA [slug ft^2]-> [kg m^2]
Izx_PA = 9.70e5*1.35581795; % Moment of Interia zx-PA [slug ft^2]-> [kg m^2]
zeta = -6.8; % Angle between Stability Axis and PA [degrees] 
I = [Ix_PA, 0,-Izx_PA;...
    0, Iy_PA,0;...
    -Izx_PA, 0, Iz_PA]; % Inertia Matrix in PA
Q_PA_SA = [cosd(zeta), 0, -sind(zeta);...
    0, 1, 0;...
    sind(zeta), 0, cosd(zeta)]; % Transformation Matrix [PA-SA]
I_SA = Q_PA_SA * I * Q_PA_SA'; % MOI in Stability axis Frame
Ix = I_SA(1,1); % Moment of Interia x-SA [kg m^2]
Iy = I_SA(2,2); % Moment of Interia y-SA [kg m^2]
Iz = I_SA(3,3); % Moment of Interia z-SA [kg m^2]
Izx = -I_SA(1,3);
CD = .040; % Coefficient of Drag
cbar = 27.31*(0.3048); % Mean Chord Length [ft]->[m]
b = 195.68*(0.3048); % Span [ft] ->[m]
S = 5500*(0.3048)^2; % Surface Area [ft^2]->[m^2]
g = 9.81; % Gravity Constant [m/s^2]
m = W/g; % Mass of Plane [kg]

% Primed Inertias
Ix_lat = (Ix*Iz-Izx^2)/Iz; % [kg m^2]
Iz_lat = (Ix*Iz-Izx^2)/Ix; % [kg m^2]
Izx_lat = Izx/(Ix*Iz-Izx^2); % [kg m^2]

%% Trim States
Vel = 518*(0.3048);% Velocity [ft/s] -> [m/s]
u0 = Vel; % Initial Velocity in x-coord - Stability Axis Frame (Trim State)
theta0 = 0; % Initial Pitch Angle [deg]

%% State Variable Derivatives
% Y (N)
Yv = (1/2)*rho*u0*S*Cy(1);
Yp = (1/4)*rho*u0*b*S*Cy(2);
Yr = (1/4)*rho*u0*b*S*Cy(3);

Y = [Yv, Yp, Yr]';

% L (N*m)
Lv = (1/2)*rho*u0*b*S*Cl(1);
Lp = (1/4)*rho*u0*b^2*S*Cl(2);
Lr = (1/4)*rho*u0*b^2*S*Cl(3);

L = [Lv, Lp, Lr]';

% N (N*m)
Nv = (1/2)*rho*u0*b*S*Cn(1);
Np = (1/4)*rho*u0*b^2*S*Cn(2);
Nr = (1/4)*rho*u0*b^2*S*Cn(3);

N = [Nv, Np, Nr]';

%% Lateral Dynamics A matrix

row1 = [Y(1)/m, Y(2)/m, (Y(3)/m-u0), g*cosd(theta0)];
row2 = [(L(1)/Ix_lat+Izx_lat*N(1)), (L(2)/Ix_lat+Izx_lat*N(2)),(L(3)/Ix_lat+Izx_lat*N(3)),0];
row3 = [(Izx_lat*L(1)+ N(1)/Iz_lat), (Izx_lat*L(2)+ N(2)/Iz_lat), (Izx_lat*L(3)+ N(3)/Iz_lat), 0];
row4 = [0, 1, tand(theta0),0];

A = [row1;row2;row3;row4;];

