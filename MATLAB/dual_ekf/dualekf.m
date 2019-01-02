clc
clear all;
close all;

dt= 0.1;
u = [0.5,1];
%state = [x, y, theta, v]
%parameter = [1/L, k/m]
%Simulation model 
x = [10;20; 0; 0; 0.015; 50];

% The system model for both EKF's
xsbar = [0; 0; 0; 0]; % initial belief 
xpbar = [2.5; 40.5]; % initial belief

Ps = diag([1, 1, 0.2, 5]); %initial belief uncertainity
Pp = 0.1*diag([20, 20]);
Rs = diag([2^2 2^2 0.2^2 10^2]); %process noise
Rp = diag([2^2 2^2]);
Q = 10*diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise
Hs = eye(4);

%H = [1 1 1 1 0 0]; %observation model - directly observing x,y,theta,v for now. 
n = 1000;
%D = [1 1 1 1];
X = zeros(6,n+1);
Xbar = zeros(6,n+1);

X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;

for k = 1:n
    x = step_sim(x,u,dt);
    z = x(1:4,1); %observing all state in simulation mode
    
    %parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
    %state prediction
    xshat = predict_state(xsbar,xphat,u,dt);
    
    G = get_jacobian(xsbar,xphat, u,dt);
    Ps = G*Ps*G' + Rs;
    
    %calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
    %update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
    %calculating Kalman gain for parameters
    Hp = measurement_jac(xsbar,xphat, u,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
    %update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,k+1) = xsbar;
    Xpbar(:,k+1) = xpbar;
    X(:,k+1) = x;
end
% plot(Xbar(1,:),Xbar(2,:),'b')
% hold on;
% plot(X(1,:),X(2,:),'r')
figure, 
plot(Xpbar(1,:),'b')
hold on;
plot(X(5,:),'r')

figure,
plot(Xpbar(2,:),'b')
hold on;
plot(X(6,:),'r')

% figure,plot(Xbar(5,:))
% figure,plot(Xbar(6,:))
