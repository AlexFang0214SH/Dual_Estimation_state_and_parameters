% Dual EKF based state and parameter estimation for autonomous vehicle control 

clc
clear all
close all

%% Estimation of COM using Kinematic model - lf, lr
dt= 0.1;

%params = [lf,lr]
%states = [x,y,phi,v];
%u = [acc,steering_angle]

x = [1;2; 0.5; 0.2];
params = [1.5; 1.5];

%The system model for both EKF's
xsbar = [0; 0; 0; 0]; % initial belief
xpbar = [5; 4.75]; % initial belief

Ps = diag([1, 1, 0.2, 2]); %initial belief uncertainity
Pp = diag([3, 5]);
Rs = diag([0.5^2 0.5^2 0.2^2 1^2]); %process noise
Rp = diag([0.25^2 0.25^2]);
Q = 10*diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise

n = 1500; %number of iteration
Hs = eye(4);

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;

sim_noise = 0.1;
gamma = 0.97;
disp('Estimating parameters from kinematic model')
for k = 1:n
    u = [0.25*cos(0.001*k)+0.05*randn(1,1),0.05*sin(0.001*k)+0.01*randn(1,1)];
    %u = [7.5e5,0.1];
    x = kinematic_step(x,params,u,dt,sim_noise);
    z = x + 0.2*randn(4,1); %observing all state in simulation mode
    
    %parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
    %state prediction
    xshat = kinematic_step(xsbar,xphat,u,dt,0);
    
    G = get_state_jacobian(xsbar,xphat,u,dt);
    Ps = G*Ps*G' + Rs;
    
    %calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
    %update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
    %calculating Kalman gain for parameters
    Hp = get_params_jacobian(xsbar,xphat,u,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
    %update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,k+1) = xsbar;
    Xpbar(:,k+1) = xpbar;
    X(:,k+1) = x;
end
figure, 
plot(Xsbar(1,:),'b')
hold on;
plot(X(1,:),'r')

figure,
plot(Xsbar(2,:),'b')
hold on;
plot(X(2,:),'r')


figure,
plot(Xsbar(3,:),'b')
hold on;
plot(X(3,:),'r')

% 
figure, 
plot(Xpbar(2,:),'r')

figure,
plot(Xpbar(1,:),'b')
disp('Done')