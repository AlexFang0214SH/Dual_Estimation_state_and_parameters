function G = get_dynparam_jacobian(x, u, p, cp, dt)

m = cp(2);                      %       m  = p(1);   /* Vehicle mass.                    */
distCog = cp(1)/2;              %       a,b = p(3);   /* Distance of COG (Length/2).  */

vel_x = x(1,1);                 %       /* x(1): Longitudinal vehicle velocity. */
vel_y = x(2,1);                 %       /* x(2): Lateral vehicle velocity. */
yawRate = x(3,1);               %       /* x(3): Yaw rate. */

slipWheel_fl = u(1);      %       u1(t) = s_FL(t)     Slip of Front Left tire (ratio).
slipWheel_fr = u(2);      %       u2(t) = s_FR(t)     Slip of Front Right tire (ratio).
steerAngle = u(5);        %       u5(t) = delta(t)    Steering angle (rad)

G(1,:) = [           (dt*cos(steerAngle)*(slipWheel_fl + slipWheel_fr))/m,          -(dt*sin(steerAngle)*(steerAngle - atan2(vel_y + distCog*yawRate, vel_x)))/m];
G(2,:) = [           (dt*sin(steerAngle)*(slipWheel_fl + slipWheel_fr))/m,           (dt*cos(steerAngle)*(steerAngle - atan2(vel_y + distCog*yawRate, vel_x)))/m];
G(3,:) = [ (dt*sin(steerAngle)*(slipWheel_fl + slipWheel_fr))/(distCog*m), (dt*cos(steerAngle)*(steerAngle - atan2(vel_y + distCog*yawRate, vel_x)))/(distCog*m)];


end