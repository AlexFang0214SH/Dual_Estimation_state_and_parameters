function [dvX_dt,dvY_dt,dyawRate_dt] = model_dynamics(u, x, param, dt)
    
    %parameters distCog_f,distCog_r, tyreCoeff_f, tyreCoeff_r, m, Jz
    %states vel_x, vel_y, yawRate
    %inputs steerAngle, slipWheel_f, slipWheel_r
    
    slipWheel_f = u(1);
    slipWheel_r = u(2);
    steerAngle = u(3);
    
    vel_x = x(1);
    vel_y = x(2);
    yawRate = x(3);
    
    distCog_f = param(1);
    distCog_r = param(2);
    tyreCoeff_f = param(3);
    tyreCoeff_r = param(4);
    m = param(5);
    Jz = param(6);
    
    slipAngle_f = (vel_y + distCog_f*yawRate)/vel_x - steerAngle;
    slipAngle_r = (vel_y - distCog_r*yawRate)/vel_x;
    
    Fy_f = -tyreCoeff_f*slipAngle_f;
    Fy_r = -tyreCoeff_r*slipAngle_r;
    Fx_f = -tyreCoeff_f*slipWheel_f;
    Fx_r = -tyreCoeff_f*slipWheel_r;
    
    torqueLong = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle))*distCog_f - Fy_r*distCog_r;
    acc_x = (Fx_f*cos(steerAngle) - Fy_f*sin(steerAngle) + Fx_r)/m;
    acc_y = (Fy_f*cos(steerAngle) + Fx_f*sin(steerAngle) + Fy_r)/m;
    
    dvX_dt = (acc_x + vel_y*yawRate)*dt;
    dvY_dt = (acc_y - vel_x*yawRate)*dt;
    dyawRate_dt = (torqueLong/Jz)*dt;

end
