function J = get_state_jacobian(state,params, u,dt)

theta = state(3);
v = state(4);
u_2 = u(2);
L = params(1);

J(1,:) = [ 1, 0, -dt*v*sin(theta),   dt*cos(theta)];
J(2,:) = [ 0, 1,  dt*v*cos(theta),   dt*sin(theta)];
J(3,:) = [ 0, 0,                1, (dt*tan(u_2))/L];
J(4,:) = [ 0, 0,                0,               1];

end