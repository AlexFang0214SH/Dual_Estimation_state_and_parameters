clc
clear all
close all

dt = 0.01;
slip_FL = 0.0015;
slip_FR = 0.0015;
slip_RR = 0;
slip_RL = 0;
steerAngle = 0.001;
u = [slip_FL slip_FR slip_RL slip_RR steerAngle];

x = [0,0,0];
%       m  = p(1);   /* Vehicle mass.                    */
%       a  = p(2);   /* Distance from front axle to COG. */
%       b  = p(3);   /* Distance from rear axle to COG.  */
%       Cx = p(4);   /* Longitudinal tire stiffness.     */
%       Cy = p(5);   /* Lateral tire stiffness.          */
%       CA = p(6);   /* Air resistance coefficient.      */
param = [1700 1.5 1.5 1.5e5 4e4 0.5];

% pos_x = 0;
% pos_y = 0;
pos = [0 0 0];
N = 200;
for i=1:N
    plot(pos(1),pos(2),'*r');    
    delta = compute_dx(x(i,:),u,param,dt);
    x(i+1,:) = x(i,:) + delta*dt;
%     slipAngle_F = u(5) - atan2(x(i,2)+param(2)*x(i,3),x(i,1));
%     slipAngle_R = atan2(param(3)*x(i,3)-x(i,2),x(i,1));
%     slip(i,:) = [slipAngle_F, slipAngle_R];
    pos = pos + x(i,:);
    hold on;
    %drawnow; 
    %pause(0.5);
end