clc
clear all
close all

load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

dt= 0.1;

x = [1;0;0];
param = [1700; 1.5; 1.5; 1.5e5; 4e4];

% The system model for both UKF's
xsbar = [0.8;0;0]; % initial belief of state
xpbar = [1700; 1.5; 1.5; 1.5e5; 4e4]; % initial belief of parameter

Ps = diag([1, 1, 0.2]); %initial belief uncertainity
Pp = diag([20, 1, 1, 50, 50]);
Rs = diag([0.1^2 0.1^2 0.2^2]); %process noise
Rp = diag([10^2 1^2 1^2 1^2 1^2]);

Q = diag([1^2 1^2 0.1^2]); %measurement noise


n = size(u1,1);
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;
X(:,1) = x;
gap = 0;
for loop = 1:n-gap
    loop
    u = u1(gap+loop,:);
    x = compute_dx(x,u,param,dt);
    z = x + 0.01*randn(3,1); %observing all state in simulation , added some noise boi)
    
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
    for k=1:size(Xp,2)%check if it is 2N or not
        Yp(:,k)=Xp(:,k); %The map of parameter is an identity
        xphat=xphat+Wm_p(k)*Yp(:,k);
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
    for k=1:size(Xs,2)%check if it is 2N or not
        %compute_dx(x(:,loop),u,param,dt);
        Ys(:,k)=compute_dx(Xs(:,k),u,xphat,dt); %The map of state 
        xshat=xshat+Wm_s(k)*Ys(:,k);
    end
    Y1s=Ys-xshat(:,ones(1,size(Xs,2)));
    Pshat=Y1s*diag(Wc_s)*Y1s'+Rs;
    
    %% Correction
    %%UKF on measurement for state
    %we can generate another set of sigma points, but I think it is not needed if we are observing the states directly
    
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
    for k=1:size(Xp,2)%check if it is 2N or not
        %compute_dx(x(:,loop),u,param,dt);
        Yzp(:,k) = compute_dx(xsbar_prev,u,Xp(:,k),dt); %Notice that I am generating expected measurement for each sigma point of parameter distribution 
        zphat=zphat+Wm_p(k)*Yzp(:,k);
    end
    Y1zp=Yzp-zphat(:,ones(1,size(Xp,2)));
    Pzhat_2=Y1zp*diag(Wc_p)*Y1zp'+Q;
    Ppz = Y1p*diag(Wc_p)*Y1zp';  %critical point, generating cross covariance
    
    Kp=Ppz*inv(Pzhat_2);
    xpbar=xphat+Kp*(z-zphat);                              %state update
    Pp=Pphat-Kp*Ppz';                                %covariance update
  
    %% Storing
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
end

% figure, 
% plot(Xpbar(2,:),'b')
% hold on;
% plot(X(6,:),'r')
% 
% figure,
% plot(Xpbar(1,:),'b')
% hold on;
% plot(X(5,:),'r')

% figure, 
% plot(Xsbar(2,:),'b')
% hold on;
% plot(X(2,:),'r')
% 
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

