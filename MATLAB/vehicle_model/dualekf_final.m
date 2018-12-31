clc
clear all
close all

load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;

x = [1;0;0];
param = [1700; 1.5; 1.5; 1.5e5; 4e4];

% The system model for both EKF's
xsbar = [0.8;0;0]; % initial belief of state
xpbar = [1650; 1.25; 1.25; 1.48e5; 4e4]; % initial belief of parameter

Ps = diag([1, 1, 0.2]); %initial belief uncertainity
Pp = diag([50, 1, 1, 500, 500]);
Rs = diag([0.1^2 0.1^2 0.2^2]); %process noise
Rp = diag([10^2 1^2 1^2 50^2 50^2]);

Q = diag([1^2 1^2 0.1^2]); %measurement noise
Hs = eye(3);

n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;
gap = 200;

for loop = 1:n-gap
    u = u1(gap+loop,:);
    x = compute_dx(x,u,param,dt);
    z = x + 0.01*randn(3,1); %observing all state in simulation , added some noise boi)
     
    %parameter prediction
    xphat = xpbar;
    Pp = Pp + Rp;
    
    %state prediction
    xshat = compute_dx(xsbar,u,xphat,dt);
    
    G = get_jacobian(xsbar,u,xphat,dt);
    Ps = G*Ps*G' + Rs;
    
    %calculating Kalman gain for state 
    Ks = Ps*Hs'*inv(Hs*Ps*Hs' + Q);
    
    %update step for state
    xsbar = xshat + Ks*(z - Hs*xshat);
    Ps = Ps - Ks*Hs*Ps;
    
    %calculating Kalman gain for parameters
    Hp = get_param_jacobian(xsbar,u,xphat,dt);
    Kp = Pp*Hp'*inv(Hp*Pp*Hp' + Q);
    
    %update step for parameters
    xpbar = xphat + Kp*(z - Hs*xshat);
    Pp = Pp - Kp*Hp*Pp;
    
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
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