function pred_sim = step_sim(state,u,dt)
pred_sim = zeros(size(state));
pred_sim(1,1) = state(1) + state(4)*cos(state(3))*dt + 0.01*randn(1,1);
pred_sim(2,1) = state(2) + state(4)*sin(state(3))*dt + 0.01*randn(1,1);
pred_sim(3,1) = state(3) + state(4)*tan(u(2))*state(5)*dt + 0.01*randn(1,1);
pred_sim(4,1) = state(4) + state(6)*u(2)*dt + 0.01*randn(1,1);
pred_sim(5,1) = state(5) ;
pred_sim(6,1) = state(6);
end