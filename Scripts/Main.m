%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Jack Lambert
% Dale Lawrence
% Aircraft Dynmaics Homework 10
% Purpose: Sets Initial Conditions for each Pertubation Case and Calls ODE45
% to plot the State Variables vs time
% Date Modefied: 4/20/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE45 Variable Allocation
%                     v_dot = z(1); % y-component of Velocity, Body Frame
%                     p_dot = z(2); % Angular Velocity about the z-axis  [rad/s]
%                     r_dot = z(3); % Angular Velocity about the z-axis [rad/s]
%                     phi_dot = z(4); % Bank Angle [rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Conditions
c1 = [10, 0,-3.0747, 1.1039]; % y-component of Velocity, Body Frame
c2 = [0, 0.1,-0.2556,-0.0051]; % Angular Velocity about the z-axis  [rad/s]
c3 = [0, 0, 0.0038,0.3015]; % Angular Velocity about the z-axis [rad/s]
c4 = [0, 0, 0.0287,0.4697]; % Bank Angle [rad]
for i = 1:4
    condition{i}= [c1(i) c2(i) c3(i) c4(i)]; 
end
%% State Variables vs. Time
time = [0 200]; % Set the time to be integrated [s]

string = ["Case I ","Case II","Case III","Case IV"]; % Title for Varying IC's
% Phugoid Response (Longer Time)
for i = 1:4
    % Calling ODE45 
    [t,z] = ode45('ODEcall_DR',time,condition{i});
    
   
    % V_E vs time
    figure
    subplot(4,1,1)
    plot(t ,z(:,1),'Linewidth',1)
    tit = sprintf('%s %s %s','State Variable of a B 747,',string(i));
    title(tit)
    ylabel('\Deltav_E [m/s]')
    
    
    % p vs time
    subplot(4,1,2)
    plot(t ,z(:,2),'Linewidth',1)
    ylabel('\Deltap [rad/s]')
    
    % r vs time
    subplot(4,1,3)
    plot(t ,z(:,3),'Linewidth',1)
    ylabel('\Deltar [rad/s]')
    
    % Phi vs time
    subplot(4,1,4)
    plot(t ,z(:,4),'Linewidth',1)
    ylabel('\DeltaPhi [rad]')
    xlabel('Time [s]')

end