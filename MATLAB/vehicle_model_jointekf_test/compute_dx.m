%   /* State equations. */
function x_next = compute_dx(x, u, p, dt)
%       /* Retrieve model parameters. */
    m = p(1);                 %       m  = p(1);   /* Vehicle mass.                    */
    distCog_f = p(2);         %       a  = p(2);   /* Distance from front axle to COG. */
    distCog_r = p(3);         %       b  = p(3);   /* Distance from rear axle to COG.  */
    tyreCoeff_x = p(4);       %       Cx = p(4);   /* Longitudinal tire stiffness.     */
    tyreCoeff_y = p(5);       %       Cy = p(5);   /* Lateral tire stiffness.          */
    Jz = m*((distCog_f+distCog_r)/2)^2;
    
    vel_x = x(1,1);             %       /* x(1): Longitudinal vehicle velocity. */
    vel_y = x(2,1);             %       /* x(2): Lateral vehicle velocity. */
    yawRate = x(3,1);           %       /* x(3): Yaw rate. */
    
    slipWheel_fl = u(1);      %       u1(t) = s_FL(t)     Slip of Front Left tire (ratio).
    slipWheel_fr = u(2);      %       u2(t) = s_FR(t)     Slip of Front Right tire (ratio).
    slipWheel_rl = u(3);      %       u3(t) = s_RL(t)     Slip of Rear Left tire (ratio).
    slipWheel_rr = u(4);      %       u4(t) = s_RR(t)     Slip of Rear Right tire (ratio).
    steerAngle = u(5);        %       u5(t) = delta(t)    Steering angle (rad)

    slipAngle_f = steerAngle - atan2((vel_y + distCog_f*yawRate),vel_x);
    slipAngle_r = -atan2((vel_y - distCog_r*yawRate),vel_x);
%     slipAngle_f = min(0.2,max(-0.2,slipAngle_f));
%     slipAngle_r = min(0.2,max(-0.2,slipAngle_r));
    slip = [slipAngle_f slipAngle_r];
    
    Fy_f = tyreCoeff_y*slipAngle_f;
    Fy_r = tyreCoeff_y*slipAngle_r;
    Fx_f = tyreCoeff_x*(slipWheel_fl+slipWheel_fr);
    Fx_r = tyreCoeff_x*(slipWheel_rl+slipWheel_rr);
    
    torqueLong = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle))*distCog_f - Fy_r*distCog_r;
    acc_x = (Fx_f*cos(steerAngle) - Fy_f*sin(steerAngle) + Fx_r)/m;
    acc_y = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle) + Fy_r)/m;
    cause = [acc_x acc_y torqueLong];
    
    dvX_dt = (acc_x + vel_y*yawRate);
    dvY_dt = (acc_y - vel_x*yawRate);
    dyawRate_dt = (torqueLong/Jz);
    
    x_next = x + [dvX_dt; dvY_dt; dyawRate_dt]*dt;
    %dragForces = p(6)*(x(1)^2); CA*(vel_x^2)
  
%   dx(1) = x(2)*x(3)+(1/p(1))*(p(4)*(u(1)+u(2))*cos(u(5))-2*p(5)*(slipAngle_F_adj)*sin(u(5))+p(4)*(u(3)+u(4))-p(6)*(x(1)^2));
%   dx(2) = (-x(1)*x(3))+(1/p(1))*(p(4)*(u(1)+u(2))*sin(u(5))+2*p(5)*(slipAngle_F_adj)*cos(u(5))+2*p(5)*slipAngle_R_adj);
%   dx(3) = (1/((((p(2)+p(3))/2)^2)*p(1)))*(p(2)*(p(4)*(u(1)+u(2))*sin(u(5))+2*p(5)*(slipAngle_F_adj)*cos(u(5)))-2*p(3)*p(5)*slipAngle_R_adj)
%           



end
