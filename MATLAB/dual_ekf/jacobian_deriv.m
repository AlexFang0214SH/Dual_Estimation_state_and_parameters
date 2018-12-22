clc
clear all
close all

syms x y theta v p1 p2 u_1 u_2 dt

eqn_1 = x + v*cos(theta)*dt;
eqn_2 = y + v*sin(theta)*dt;
eqn_3 = theta + v*tan(u_1)*dt*p1;
eqn_4 = v + p2*u_2*dt;

J(1,:) = jacobian(eqn_1,[p1, p2]);
J(2,:) = jacobian(eqn_2,[p1, p2]);
J(3,:) = jacobian(eqn_3,[p1, p2]);
J(4,:) = jacobian(eqn_4,[p1, p2]);

%J(5,:) = jacobian(eqn_5,[x,y,theta,v,p1,p2]);
%J(6,:) = jacobian(eqn_6,[x,y,theta,v,p1,p2]);

J