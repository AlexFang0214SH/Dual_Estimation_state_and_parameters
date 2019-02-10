% Dual EKF based state and parameter estimation for autonomous vehicle control 

clc
clear all
close all

%% Estimation of COM using Kinematic model - lf, lr
dt= 0.1;

%params = [lf,lr]
%states = [x,y,phi,v];
%u = [acc,steering_angle]

x = [0;0; 0; 0];
params = [1.5; 1.5];

%The system model for both EKF's
xsbar = [0; 0; 0; 0]; % initial belief
xpbar = [3; 2]; % initial belief

Ps = diag([1, 1, 0.2, 2]); %initial belief uncertainity
Pp = diag([3, 3]);
Rs = diag([0.5^2 0.5^2 0.2^2 1^2]); %process noise
Rp = diag([0.25^2 0.25^2]);
Q = diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise

n = 1500; %number of iteration

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;

sim_noise = 0.05;
gamma = 0.99;
disp('Estimating parameters from kinematic model')
for k = 1:n
    u = [2.5*cos(0.001*k),0.15*sin(0.001*k)];
    %u = [7.5e5,0.1];
    x = kinematic_step(x,params,u,dt,sim_noise);
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
% plot(Xsbar(3,:),'b')
% hold on;
% plot(X(3,:),'r')
% 
figure, 
plot(Xpbar(2,:),'b')

figure,
plot(Xpbar(1,:),'b')
disp('Done')
pause()
% weight_vector = ones(1,size(X,2))*1700;
% length_vector = ones(1,size(X,2))*3;


% figure, 
% plot(Xpbar(2,:),'b')
% hold on;
% plot(weight_vector,'r')
% 
% figure,
% plot(Xpbar(1,:),'b')
% Estimation of Tire Co-efficients using Dynamic model 
clear all
close all
time = [];
current_time = 0;
% loading control input data
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;
pre_computed_params = [1.5037,1698];
orig_params = [1.5,1700];
% dynamic states = vx vy yaw_rate ax ay  
% observing vx yaw_rate ax ay 
x = [1;0;0];
% tire stiffness co-efficients 
param = [20; 5];

% The system model for both EKF's
xsbar = [0.25;0;0]; % initial belief of state
xpbar = [15; 7]; % initial belief of parameter

Ps = diag([1, 1, 0.2]); %initial belief uncertainity
Pp = diag([10 5]);
Rs = 0.1*diag([0.2^2 0.2^2 0.05^2]); %process noise
Rp = diag([3.5^2 1.25^2]);

Q = diag([0.1^2 0.1^2 0.1^2]); %measurement noise
Hs = eye(3);
gamma = 0.95;
n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;

disp('Estimating parameters from dynamic model')

for loop = 1:n
    time = [time;current_time];
    current_time = current_time+dt;
    u = u1(loop,:);
    x = predict_state(x,u,param,orig_params,dt,0.001);
    z = Hs*x + 0.01*randn(3,1); 
    
%     parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
%     state prediction
    xshat = predict_state(xsbar,u,xphat,pre_computed_params,dt,0);
    
    G = get_dynstate_jacobian(xsbar,u,xphat,pre_computed_params,dt);
    Ps = G*Ps*G' + Rs;
    
%     calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
%     update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
%     calculating Kalman gain for parameters
    Hp = get_dynparam_jacobian(xsbar,u,xphat,pre_computed_params,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
%     update step for parameters
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

% disp('Done')
% 
% plot(time(1:end-28),Xsbar(3,30:end),'b')
% hold on;
% plot(time(1:end-28),X(3,30:end),'r')
