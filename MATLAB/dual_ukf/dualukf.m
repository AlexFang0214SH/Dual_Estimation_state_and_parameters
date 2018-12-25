clc
clear all
close all

dt= 0.1;
u = [0.5,0.01];
%state = [x, y, theta, v]
%parameter = [1/L, k/m]
%Simulation model
x = [10;20; 0; 0; 0.015; 50];

% The system model for both UKF's
xsbar = [0; 0; 0; 0]; % initial belief
xpbar = [2.5; 40.5]; % initial belief

Ps = diag([1, 1, 0.2, 5]); %initial belief uncertainity
Pp = diag([20, 20]);
Rs = diag([2^2 2^2 0.2^2 10^2]); %process noise
Rp = diag([2^2 2^2]);
Q = 10*diag([1^2 1^2 0.1^2 0.1^2]); %measurement noise


n = 1000;
X = zeros(6,n+1);
X(:,1) = x;
Xsbar(:,1) = xsbar;
Xpbar(:,1) = xpbar;

for loop = 1:n
    x = step_sim(x,u,dt);
    z = x(1:4,1)+0.1*randn(4,1); %observing all state in simulation , added some noise boi)
    
    %% Prediction
    %%Parameter
    %UKF defining tuning variables and weights
    L_p=numel(xpbar);
    alpha_p=1e-3; %tune here.. for parameters
    ki_p=0;
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
    
    %parameter prediction using unscented transform
    xshat = zeros(L_s,1);
    for k=1:size(Xs,2)%check if it is 2N or not
        Ys(:,k)=predict_state(Xs(:,k),xphat,u,dt); %The map of state 
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
        Yzp(:,k) = predict_state(xsbar,Xp(:,k),u,dt); %Notice that I am generating expected measurement for each sigma point of parameter distribution 
        zphat=zphat+Wm_p(k)*Yzp(:,k);
    end
    Y1zp=Yzp-zphat(:,ones(1,size(Xp,2)));
    Ppz = Y1p*diag(Wc_p)*Y1zp';  %critical point, generating cross covariance
    
    Kp=Ppz*inv(Pzhat);
    xpbar=xphat+Kp*(z-zhat);                              %state update
    Pp=Pphat-Kp*Ppz';                                %covariance update
  
    %% Storing
    Xsbar(:,loop+1) = xsbar;
    Xpbar(:,loop+1) = xpbar;
    X(:,loop+1) = x;
end

figure, 
plot(Xpbar(2,:),'b')
hold on;
plot(X(6,:),'r')

figure,
plot(Xpbar(1,:),'b')
hold on;
plot(X(5,:),'r')