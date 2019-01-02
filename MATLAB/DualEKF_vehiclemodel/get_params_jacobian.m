function G = get_params_jacobian(state,params, u,dt)

v = state(4);
u_1 = u(1);
u_2 = u(2);
L = params(1);
m = params(2);

G(1,:) = [                    0,             0];
G(2,:) = [                    0,             0];
G(3,:) = [ -(dt*v*tan(u_2))/L^2,             0];
G(4,:) = [                    0, -(dt*u_1)/m^2];

end