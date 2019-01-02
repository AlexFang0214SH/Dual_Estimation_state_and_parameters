%% Estimation of Tire Co-efficients using Dynamic model using 5 states (not working)

clear all
%loading control input data
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;
pre_computed_params = [1.5,1700];
%dynamic states = vx vy yaw_rate ax ay  
%observing vx yaw_rate ax ay 
x = [1;0;0;0.01;0];
%tire stiffness co-efficients 
param = [20; 5];

% The system model for both EKF's
xsbar = [0.8;0;0;0.05;0]; % initial belief of state
xpbar = [15; 4]; % initial belief of parameter

Ps = diag([1, 1, 0.2, 0.5, 0.5]); %initial belief uncertainity
Pp = diag([10 5]);
Rs = diag([0.2^2 0.2^2 0.05^2 0.1^2 0.1^2]); %process noise
Rp = diag([2.5^2 1.25^2]);

Q = diag([1^2 0.1^2 0.5^2 0.5^2]); %measurement noise
Hs = [1 0 0 0 0;0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1];
%Not measuring lateral velocity due to practical feasibility

n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;

for loop = 1:n
    u = u1(loop,:);
    x = predict_state(x,u,param,pre_computed_params,dt);
    z = Hs*x + 0.01*randn(4,1); 
    
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
