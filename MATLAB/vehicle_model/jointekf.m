clc
clear all
close all

load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;

x = [1;0;0];
param = [1700; 1.5; 1.5; 1.5e5; 4e4];

% The system model for both EKF's
xbar = [0.8;0;0;1680; 1.25; 1.25; 1.48e5; 4e4]; % initial belief of parameter

P = diag([1, 1, 0.2, 50, 1, 1, 500, 500]);
R = diag([0.1^2 0.1^2 0.2^2 10^2 1^2 1^2 50^2 50^2]); %process noise

Q = diag([1^2 1^2 0.1^2]); %measurement noise
H = [eye(3) zeros(3,5)];

n = size(u1,1);
Xbar(:,1) = xbar;
X(:,1) = x;

for loop = 1:n
    u = u1(loop,:);
    x = compute_dx(x,u,param,dt);
    z = x + 0.01*randn(3,1); %observing all state in simulation , added some noise boi)
    
    %predict step
    pstate = compute_dx(xbar(1:3,1),u,xbar(4:8,1),dt);
    xhat = [pstate;xbar(4:8)];
    G = get_jacobian_combined(xbar(1:3,1),u,xbar(4:8,1),dt);
    P = G*P*G' + R;
    
    %calculating Kalman gain 
    K = P*H'*inv(H*P*H' + Q);
    
    %update step
    xbar = xhat + K*(z - H*xhat);
    P = P - K*H*P;
    Xbar(:,loop+1) = xbar;
    X(:,loop+1) = x;
end  
