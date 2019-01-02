function J = get_dynstate_jacobian(x, u, p, cp, dt)

m = cp(2);                      %       m  = p(1);   /* Vehicle mass.                    */
distCog = cp(1)/2;              %       a,b = p(3);   /* Distance of COG (Length/2).  */
tyreCoeff_y = p(2)*10000;       %       Cy = p(5);   /* Lateral tire stiffness.          */

vel_x = x(1,1);                 %       /* x(1): Longitudinal vehicle velocity. */
vel_y = x(2,1);                 %       /* x(2): Lateral vehicle velocity. */
yawRate = x(3,1);               %       /* x(3): Yaw rate. */
steerAngle = u(5);        %       u5(t) = delta(t)    Steering angle (rad)

J(1,:) = [          1 - (dt*tyreCoeff_y*sin(steerAngle)*(vel_y + distCog*yawRate))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2)), dt*(yawRate + (tyreCoeff_y*vel_x*sin(steerAngle))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2))),  dt*(vel_y + (distCog*tyreCoeff_y*vel_x*sin(steerAngle))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2)))];
J(2,:) = [ -dt*(yawRate - (tyreCoeff_y*cos(steerAngle)*(vel_y + distCog*yawRate))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2))),         1 - (dt*tyreCoeff_y*vel_x*cos(steerAngle))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2)), -dt*(vel_x + (distCog*tyreCoeff_y*vel_x*cos(steerAngle))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2)))];
J(3,:) = [      (dt*tyreCoeff_y*cos(steerAngle)*(vel_y + distCog*yawRate))/(distCog*m*((vel_y + distCog*yawRate)^2 + vel_x^2)),    -(dt*tyreCoeff_y*vel_x*cos(steerAngle))/(distCog*m*((vel_y + distCog*yawRate)^2 + vel_x^2)),                1 - (dt*tyreCoeff_y*vel_x*cos(steerAngle))/(m*((vel_y + distCog*yawRate)^2 + vel_x^2))];

end