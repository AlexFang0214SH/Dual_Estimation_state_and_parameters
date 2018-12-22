function state_hat = predict_state(state,param,u,dt)
state_hat = zeros(size(state));
state_hat(1,1) = state(1) + state(4)*cos(state(3))*dt;
state_hat(2,1) = state(2) + state(4)*sin(state(3))*dt;
state_hat(3,1) = state(3) + state(4)*tan(u(2))*param(1)*dt;
state_hat(4,1) = state(4) + param(2)*u(2)*dt;
end