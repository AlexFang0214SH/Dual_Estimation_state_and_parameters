function pred_sim = step_sim(state,params,u,dt,noise_level)
pred_sim = zeros(size(state));
pred_sim(1,1) = state(1) + state(4)*cos(state(3))*dt + noise_level*randn(1,1);
pred_sim(2,1) = state(2) + state(4)*sin(state(3))*dt + noise_level*randn(1,1);
pred_sim(3,1) = state(3) + state(4)*tan(u(2))*dt/params(1) + noise_level*randn(1,1);
pred_sim(4,1) = state(4) + u(1)*dt/params(2) + noise_level*randn(1,1);
end