%% Author: Jack Lambert
% ASEN 3128
% Purpose: Function for ODE45 to call to calculate the State variables 
% u_dot, w_dot, q_dot, and theta_dot for the PWD Approximation. This function
% uses the simplified assumptions for the Linearized Longitudinal Dynamics Set
% Last Edited: 4/9/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dydt] = ODEcall_DR(t,y)

v_dot = y(1); % y-component of Velocity, Body Frame
p_dot = y(2); % roll-rate 
r_dot = y(3); % yaw rate
theta_dot = y(4); % Pitch Angle 

%% State Variable Matrix for Linearized Longitudinal Set
[A] = Amat(); % A matrix function based on plane and parameters
State = [v_dot, p_dot, r_dot, theta_dot]'; % Couple State Variables in Long. Set
var = A*State; % Couple State Variables in Long. Set
%% Solving for State Variables in the Linearized Longitudinal Set
dydt(1) = var(1); % v
dydt(2) = var(2); % p
dydt(3) = var(3); % v
dydt(4) = var(4); % p

dydt = dydt'; % Inverts for ODE45   
end