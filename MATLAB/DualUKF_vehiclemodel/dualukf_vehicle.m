% Dual UKF based state and parameter estimation for autonomous vehicle control

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
    
    %% Prediction
    %%Parameter
    %UKF defining tuning variables and weights
    L_p=numel(xpbar);
    alpha_p=1e-2; %tune here.. for parameters
    ki_p=1;
    beta_p=2;
    lambda_p=alpha_p^2*(L_p+ki_p)-L_p;
    c_p=L_p+lambda_p;
    Wm_p=[lambda_p/c_p 0.5/c_p+zeros(1,2*L_p)];           %weights for means
    Wc_p=Wm_p;
    Wc_p(1)=Wc_p(1)+(1-alpha_p^2+beta_p);               %weights for covariance
    c_p=sqrt(c_p);
    
    %UKF defining sigma points Xt-1
    Xp=sigmas(xpbar,Pp,c_p);                            %sigma points around x
    
    %parameter prediction using unscented transform
    xphat = zeros(L_p,1);
    for l=1:size(Xp,2)
        Yp(:,l)=Xp(:,l); %The map of parameter is an identity
        xphat=xphat+Wm_p(l)*Yp(:,l);
    end
    Y1p=Yp-xphat(:,ones(1,size(Xp,2)));
    Pphat=Y1p*diag(Wc_p)*Y1p'+Rp;
    
    %%State
    %UKF defining tuning variables and weights
    L_s=numel(xsbar);
    alpha_s=1e-3;
    ki_s=0;
    beta_s=2;
    lambda_s=alpha_s^2*(L_s+ki_s)-L_s;
    c_s=L_s+lambda_s;
    Wm_s=[lambda_s/c_s 0.5/c_s+zeros(1,2*L_s)];           %weights for means
    Wc_s=Wm_s;
    Wc_s(1)=Wc_s(1)+(1-alpha_s^2+beta_s);               %weights for covariance
    c_s=sqrt(c_s);
    
    %UKF defining sigma points Xt-1
    Xs=sigmas(xsbar,Ps,c_s);                            %sigma points around x
    xsbar_prev = xsbar;
    %parameter prediction using unscented transform
    xshat = zeros(L_s,1);
    for l=1:size(Xs,2)%check if it is 2N or not
        Ys(:,l)=step_sim(Xs(:,l),xphat,u,dt,0); %The map of state
        xshat=xshat+Wm_s(l)*Ys(:,l);
    end
    Y1s=Ys-xshat(:,ones(1,size(Xs,2)));
    Pshat=Y1s*diag(Wc_s)*Y1s'+Rs;
    
    
    
    %% Correction
    %%UKF on measurement for state
    
    %measurement prediction using unscented transform
    zhat = xshat;
    Yz = Ys; %Map of states and measurement is identity in this case, will change if we observe less states (add observability)
    Y1z=Yz-zhat(:,ones(1,size(Yz,2)));
    Pzhat=Y1z*diag(Wc_s)*Y1z'+Q; %measurement noise
    
    %%State correction
    Psz=Y1s*diag(Wc_s)*Y1z'; %line number 10 of the algo
    
    Ks=Psz*inv(Pzhat);
    xsbar=xshat+Ks*(z-zhat);                              %state update
    Ps=Pshat-Ks*Psz';                                %covariance update
    
    %Hope it is correct.. !!
    
    %%UKF on measurement for parameters
    zphat = zeros(L_s,1);
    for l=1:size(Xp,2)%check if it is 2N or not
        Yzp(:,l) = step_sim(xsbar_prev,Xp(:,l),u,dt,0); %Notice that I am generating expected measurement for each sigma point of parameter distribution
        zphat=zphat+Wm_p(l)*Yzp(:,l);
    end
    Y1zp=Yzp-zphat(:,ones(1,size(Xp,2)));
    Pzhat_2=Y1zp*diag(Wc_p)*Y1zp'+Q;
    Ppz = Y1p*diag(Wc_p)*Y1zp';  %critical point, generating cross covariance
    
    Kp=Ppz*inv(Pzhat_2);
    xpbar=xphat+Kp*(z-zphat);                              %state update
    Pp=Pphat-Kp*Ppz';                                %covariance update
    
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
pre_computed_params = [1.4537,1680];
orig_params = [1.5,1700];
%dynamic states = vx vy yaw_rate ax ay
%observing vx yaw_rate ax ay
x = [0.01;0;0];
%tire stiffness co-efficients
param = [20; 5];

% The system model for both EKF's
xsbar = [0.05;0;0]; % initial belief of state
xpbar = [30; 2]; % initial belief of parameter

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
    %% Prediction
    %%Parameter
    %UKF defining tuning variables and weights
    L_p=numel(xpbar);
    alpha_p=1e-3; %tune here.. for parameters
    ki_p=1;
    beta_p=2;
    lambda_p=alpha_p^2*(L_p+ki_p)-L_p;
    c_p=L_p+lambda_p;
    Wm_p=[lambda_p/c_p 0.5/c_p+zeros(1,2*L_p)];           %weights for means
    Wc_p=Wm_p;
    Wc_p(1)=Wc_p(1)+(1-alpha_p^2+beta_p);               %weights for covariance
    c_p=sqrt(c_p);
    
    %UKF defining sigma points Xt-1
    Xp=sigmas(xpbar,Pp,c_p);                            %sigma points around x
    
    %parameter prediction using unscented transform
    xphat = zeros(L_p,1);
    for l=1:size(Xp,2)
        Yp(:,l)=Xp(:,l); %The map of parameter is an identity
        xphat=xphat+Wm_p(l)*Yp(:,l);
    end
    Y1p=Yp-xphat(:,ones(1,size(Xp,2)));
    Pphat=Y1p*diag(Wc_p)*Y1p'+Rp;
    
    %%State
    %UKF defining tuning variables and weights
    L_s=numel(xsbar);
    alpha_s=1e-3;
    ki_s=0;
    beta_s=2;
    lambda_s=alpha_s^2*(L_s+ki_s)-L_s;
    c_s=L_s+lambda_s;
    Wm_s=[lambda_s/c_s 0.5/c_s+zeros(1,2*L_s)];           %weights for means
    Wc_s=Wm_s;
    Wc_s(1)=Wc_s(1)+(1-alpha_s^2+beta_s);               %weights for covariance
    c_s=sqrt(c_s);
    
    %UKF defining sigma points Xt-1
    Xs=sigmas(xsbar,Ps,c_s);                            %sigma points around x
    xsbar_prev = xsbar;
    %parameter prediction using unscented transform
    xshat = zeros(L_s,1);
    for l=1:size(Xs,2)%check if it is 2N or not
        Ys(:,l)=predict_state(Xs(:,l),u,xphat,pre_computed_params,dt); %The map of state
        xshat=xshat+Wm_s(l)*Ys(:,l);
    end
    Y1s=Ys-xshat(:,ones(1,size(Xs,2)));
    Pshat=Y1s*diag(Wc_s)*Y1s'+Rs;
    
    
    
    %% Correction
    %%UKF on measurement for state
    
    %measurement prediction using unscented transform
    zhat = xshat;
    Yz = Ys; %Map of states and measurement is identity in this case, will change if we observe less states (add observability)
    Y1z=Yz-zhat(:,ones(1,size(Yz,2)));
    Pzhat=Y1z*diag(Wc_s)*Y1z'+Q; %measurement noise
    
    %%State correction
    Psz=Y1s*diag(Wc_s)*Y1z'; %line number 10 of the algo
    
    Ks=Psz*inv(Pzhat);
    xsbar=xshat+Ks*(z-zhat);                              %state update
    Ps=Pshat-Ks*Psz';                                %covariance update
    
    %Hope it is correct.. !!
    
    %%UKF on measurement for parameters
    zphat = zeros(L_s,1);
    for l=1:size(Xp,2)%check if it is 2N or not
        Yzp(:,l) = predict_state(xsbar_prev,u,Xp(:,l),pre_computed_params,dt); %Notice that I am generating expected measurement for each sigma point of parameter distribution
        zphat=zphat+Wm_p(l)*Yzp(:,l);
    end
    Y1zp=Yzp-zphat(:,ones(1,size(Xp,2)));
    Pzhat_2=Y1zp*diag(Wc_p)*Y1zp'+Q;
    Ppz = Y1p*diag(Wc_p)*Y1zp';  %critical point, generating cross covariance
    
    Kp=Ppz*inv(Pzhat_2);
    xpbar=xphat+Kp*(z-zphat);                              %state update
    Pp=Pphat-Kp*Ppz';                                %covariance update
    
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
end
figure,
plot(Xpbar(2,:),'b')

figure,
plot(Xpbar(1,:),'b')
disp('Done')