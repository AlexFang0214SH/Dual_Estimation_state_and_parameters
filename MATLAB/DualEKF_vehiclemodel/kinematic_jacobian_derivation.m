clc
clear all
close all

syms x y theta v m L u_1 u_2 dt

eqn_1 = x + v*cos(theta)*dt;
eqn_2 = y + v*sin(theta)*dt;
eqn_3 = theta + v*tan(u_2)*dt/L;
eqn_4 = v + u_1*dt/m;

J(1,:) = jacobian(eqn_1,[x, y, theta, v]);
J(2,:) = jacobian(eqn_2,[x, y, theta, v]);
J(3,:) = jacobian(eqn_3,[x, y, theta, v]);
J(4,:) = jacobian(eqn_4,[x, y, theta, v]);


G(1,:) = jacobian(eqn_1,[L, m]);
G(2,:) = jacobian(eqn_2,[L, m]);
G(3,:) = jacobian(eqn_3,[L, m]);
G(4,:) = jacobian(eqn_4,[L, m]);

