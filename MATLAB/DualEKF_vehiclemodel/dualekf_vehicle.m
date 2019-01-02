% Dual EKF based state and parameter estimation for autonomous vehicle control 

clc
clear all
close all

% Estimation of Mass and Length using Kinematic model 
% v = longitudinal velocity
% x_dot = v*cos(heading_angle)
% y_dot = v*sin(heading_angle)
% heading_angle = 
% v = 

dt= 0.1;

% state = [x, y, theta, v]
% parameter = [L, m]
% Simulation model 
x = [0;0; 0; 0];
params = [3; 1700];

%The system model for both EKF's
xsbar = [0; 0; 0; 0]; % initial belief 
xpbar = [2; 1200]; % initial belief

Ps = diag([1, 1, 0.2, 5]); %initial belief uncertainity
Pp = diag([5, 500]);
Rs = diag([2^2 2^2 0.2^2 1^2]); %process noise
Rp = diag([1^2 10^2]);
Q = 10*diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise
Hs = eye(4);

n = 1000; %number of iteration

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;

sim_noise = 0.1;
disp('Estimating parameters from kinematic model')
for k = 1:n
    %u(1) - Throttle
    %u(2) - Steering Angle
    u = [0.75e5*abs(cos(0.1*k)),0.25*sin(0.1*k)];
    %u = [7.5e5,0.1];
    x = step_sim(x,params,u,dt,sim_noise);
    z = x + 0.1*randn(4,1); %observing all state in simulation mode
    
    %parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
    %state prediction
    xshat = step_sim(xsbar,xphat,u,dt,0);
    
    G = get_state_jacobian(xsbar,xphat, u,dt);
    Ps = G*Ps*G' + Rs;
    
    %calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
    %update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
    %calculating Kalman gain for parameters
    Hp = get_params_jacobian(xsbar,xphat, u,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
    %update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,k+1) = xsbar;
    Xpbar(:,k+1) = xpbar;
    X(:,k+1) = x;
end
plot(Xsbar(1,:),'b')
hold on;
plot(X(1,:),'r')
figure, 
plot(Xpbar(2,:),'b')

figure,
plot(Xpbar(1,:),'b')

disp('Done')
pause()
%% Estimation of Tire Co-efficients using Dynamic model 

clear all
%loading control input data
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;
pre_computed_params = [1.5037,1698];
orig_params = [1.5,1700];
%dynamic states = vx vy yaw_rate ax ay  
%observing vx yaw_rate ax ay 
x = [0.01;0;0];
%tire stiffness co-efficients 
param = [20; 5];

% The system model for both EKF's
xsbar = [0.05;0;0]; % initial belief of state
xpbar = [15; 4]; % initial belief of parameter

Ps = diag([1, 1, 0.2]); %initial belief uncertainity
Pp = diag([5 3]);
Rs = diag([0.2^2 0.2^2 0.05^2]); %process noise
Rp = diag([5.5^2 1.25^2]);

Q = diag([0.1^2 0.1^2 0.1^2]); %measurement noise
Hs = eye(3);

n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;

disp('Estimating parameters from dynamic model')

for loop = 1:n
    u = 2*u1(loop,:);
    x = predict_state(x,u,param,orig_params,dt);
    z = Hs*x + 0.01*randn(3,1); 
    
    %parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
    %state prediction
    xshat = predict_state(xsbar,u,xphat,pre_computed_params,dt);
    
    G = get_dynstate_jacobian(xsbar,u,xphat,pre_computed_params,dt);
    Ps = G*Ps*G' + Rs;
    
    %calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
    %update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
    %calculating Kalman gain for parameters
    Hp = get_dynparam_jacobian(xsbar,u,xphat,pre_computed_params,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
    %update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
end
figure,
plot(Xpbar(1,:));
figure,
plot(Xpbar(2,:));

disp('Done')