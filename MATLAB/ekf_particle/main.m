%   Noise configuration
a = 0.1;
R_SIM = diag([0.001, 0.001, 0.01 * pi / 180, 0.01])*a;
Q_SIM = diag([0.01, 0.01, 0.1 * pi / 180, 0.01])*a;

R_MODEL = diag([0.001, 0.001, 0.01 * pi / 180, 0.01])*a;
Q_MODEL = diag([0.02, 0.02, 0.2 * pi / 180, 0.015])*a;



P = diag([1,1,0.02,0.01]);

Qt = P;


%   Define simulation model
cycle_sim = Bicycle();
cycle_sim.R = R_SIM;
cycle_sim.Q = Q_SIM;


%   Define estimation model
cycle_model = Bicycle();
cycle_model.R = R_MODEL;
cycle_model.Q = Q_MODEL;


N = 100000;
X = zeros(4, N+1);
X(:,1) = [0; 0; pi/6; 0];


P = [1; 50];

Xhat = X;
Phat = zeros(2, N+1);
Phat(:, 1) = [1;50];


ekf = EKF();

ekf.R = R_MODEL;
ekf.Q = Q_MODEL;
ekf.Qt = Qt;

for i=2:N+1
    x = X(:, i-1);
    
    u = [(4. - x(4,1) + 2.0 * sin(0.01 * i) ) * 0.05; 0.1 + 0.2 * rand()];
    
    [nX, nY] = cycle_sim.step(x, P, u);
    X(:, i) = nX;
    
    xhat = Xhat(:, i-1);
    [xEst, pEst, ekf] = ekf.run(cycle_model, x, P, u, nY);
    Xhat(:, i) = xEst;
    Phat(:, i) = pEst;
    
end

figure,
plot(X(1,:), X(2,:), 'b');
hold on
plot(Xhat(1,:), Xhat(2,:), 'r');

figure,
plot(Phat(1,:), 'r');

figure,
plot(Phat(2,:), 'r');


