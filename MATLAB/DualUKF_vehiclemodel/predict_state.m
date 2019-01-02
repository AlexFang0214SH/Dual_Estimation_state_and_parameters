%   /* State equations. */
function x_next = predict_state(x, u, p, cp, dt)
%       /* Retrieve model parameters. */
m = cp(2);                      %       m  = p(1);   /* Vehicle mass.                    */
distCog = cp(1)/2;              %       a,b = p(3);   /* Distance of COG (Length/2).  */
tyreCoeff_x = p(1)*10000;       %       Cx = p(4);   /* Longitudinal tire stiffness.     */
tyreCoeff_y = p(2)*10000;       %       Cy = p(5);   /* Lateral tire stiffness.          */
Jz = m*((distCog+distCog)/2)^2;

vel_x = x(1,1);             %       /* x(1): Longitudinal vehicle velocity. */
vel_y = x(2,1);             %       /* x(2): Lateral vehicle velocity. */
yawRate = x(3,1);           %       /* x(3): Yaw rate. */

slipWheel_fl = u(1);      %       u1(t) = s_FL(t)     Slip of Front Left tire (ratio).
slipWheel_fr = u(2);      %       u2(t) = s_FR(t)     Slip of Front Right tire (ratio).
slipWheel_rl = u(3);      %       u3(t) = s_RL(t)     Slip of Rear Left tire (ratio).
slipWheel_rr = u(4);      %       u4(t) = s_RR(t)     Slip of Rear Right tire (ratio).
steerAngle = u(5);        %       u5(t) = delta(t)    Steering angle (rad)

slipAngle_f = steerAngle - atan2((vel_y + distCog*yawRate),vel_x);
slipAngle_r = -atan2((vel_y - distCog*yawRate),vel_x);
%     slipAngle_f = min(0.2,max(-0.2,slipAngle_f));
%     slipAngle_r = min(0.2,max(-0.2,slipAngle_r));
slip = [slipAngle_f slipAngle_r];

Fy_f = tyreCoeff_y*slipAngle_f;
Fy_r = tyreCoeff_y*slipAngle_r;
Fx_f = tyreCoeff_x*(slipWheel_fl+slipWheel_fr);
Fx_r = tyreCoeff_x*(slipWheel_rl+slipWheel_rr);

torqueLong = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle))*distCog - Fy_r*distCog;
acc_x = (Fx_f*cos(steerAngle) - Fy_f*sin(steerAngle) + Fx_r)/m;
acc_y = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle) + Fy_r)/m;
cause = [acc_x acc_y torqueLong];

dvX_dt = (acc_x + vel_y*yawRate);
dvY_dt = (acc_y - vel_x*yawRate);
dyawRate_dt = (torqueLong/Jz);

x_next = x + [dvX_dt; dvY_dt; dyawRate_dt]*dt;

end
