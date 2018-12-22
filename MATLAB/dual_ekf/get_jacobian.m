function G = get_jacobian(state,param, u,dt)
theta = state(3);
v = state(4);
u_1 = u(1);
u_2 = u(2);
p1 = param(1);
p2 = param(2);

G(1,:) = [ 1, 0, -dt*v*sin(theta),  dt*cos(theta)];
G(2,:) = [ 0, 1,  dt*v*cos(theta),  dt*sin(theta)];
G(3,:) = [ 0, 0,                1, dt*p1*tan(u_1)];
G(4,:) = [ 0, 0,                0,              1];
end