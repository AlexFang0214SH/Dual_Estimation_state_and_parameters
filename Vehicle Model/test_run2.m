clc
clear all
close all

dt = 0.1;
% slip_FL = 0.0015;
% slip_FR = 0.0015;
% slip_RR = 0;
% slip_RL = 0;
% steerAngle = -0.01;
% u = [slip_FL slip_FR slip_RL slip_RR steerAngle];
load(fullfile(matlabroot, 'toolbox', 'ident', 'iddemos', 'data', 'vehicledata'));
u = u1;

x = [0,0,0];
%       m  = p(1);   /* Vehicle mass.                    */
%       a  = p(2);   /* Distance from front axle to COG. */
%       b  = p(3);   /* Distance from rear axle to COG.  */
%       Cx = p(4);   /* Longitudinal tire stiffness.     */
%       Cy = p(5);   /* Lateral tire stiffness.          */
%       CA = p(6);   /* Air resistance coefficient.      */
param = [1700 1.5 1.5 1.5e5 4e4 0.5];

pos = [0 0 0];
N = size(u1,1);

for i=1:N
    plot(pos(1),pos(2),'*r');    
    delta = compute_dx(x(i,:),u(i,:),param);
    x(i+1,:) = x(i,:) + delta*dt ;

    pos = pos + x(i,:)*dt;
    hold on;
    %drawnow; 
    %pause(0.5);
end