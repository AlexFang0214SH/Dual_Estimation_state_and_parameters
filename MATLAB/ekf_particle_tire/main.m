
%loading control input data
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));

%   Noise configuration
a = 0.1;
R_SIM = diag([0.01, 0.01, 0.01])*a*4;
Q_SIM = diag([0.01, 0.01, 0.01])*a*10;


R_MODEL = diag([0.01, 0.01, 0.01])*a*4;
Q_MODEL = diag([0.02, 0.02, 0.015])*a*10;



Qt = diag([0.01, 0.01, 0.01])*10;



%   Define simulation model
cycle_sim = Tires();
cycle_sim.R = R_SIM;
cycle_sim.Q = Q_SIM;


%   Define estimation model
cycle_model = Tires();
cycle_model.R = R_MODEL;
cycle_model.Q = Q_MODEL;


N = 600;
X = zeros(3, N+1);
X(:,1) = [0.01; 0; 0];


P = [1700;1.5; 20e5;5e5];

Xhat = X;
Phat = zeros(4, N+1);
Phat(:, 1) = P;

dif_mag = zeros(1, N+1);

ekf = EKF();

ekf.R = R_MODEL;
ekf.Q = Q_MODEL;
ekf.Qt = Qt;

for i=2:N+1
    x = X(:, i-1);
    
    %u = [(4. - x(4,1) + 2.0 * sin(0.01 * i) ) * 0.05; 0.1 + 0.2 * rand()];
    %u = [0.02; 0.01];
    u = u2(i, :)';
    
    [nX, nY] = cycle_sim.step(x, P, u);
    X(:, i) = nX;
    
    xhat = Xhat(:, i-1);
    [xEst, pEst, ekf] = ekf.run(cycle_model, x, P, u, nY);
    Xhat(:, i) = xEst;
    Phat(:, i) = pEst;
    
    dif_mag(1, i) = ekf.diffusion_mag;
    
end

% figure,
% plot(X(1,:), 'b');
% hold on
% plot(Xhat(1,:), 'r');
% 
% 
% figure,
% plot(X(2,:), 'b');
% hold on
% plot(Xhat(2,:), 'r');
% 
% figure,
% plot(X(3,:), 'b');
% hold on
% plot(Xhat(3,:), 'r');
% 
% 
figure,
plot(Phat(3,:), 'r');

figure,
plot(Phat(4,:), 'r');

% figure,
% plot(dif_mag(1,:), 'r');







