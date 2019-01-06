classdef Tires
    properties
        dt = 0.1;
        nx = 3;
        ny = 3;
        R
        Q
        
        C = [   1 0 0;
                0 1 0;
                0 0 1; ];
            
        D = [   0 0 0;
                0 0 0;
                0 0 0;];
    end
    
    methods    
        function [xNext, y] = step(self, x, p, u)
            [xNext, y] = predict(self, x, p, u);
            
            process_noise = self.R * randn(self.nx, 1);
            measurement_noise = self.Q * randn(self.ny, 1);
            
            xNext = xNext + process_noise;
            y = y + measurement_noise;
            
        end
    
        function [xNext, y] = predict(self, x, p, u)
            %   state = [Vx, Vy, yaw].Transpose
            %   y = [Vx, Vy, yaw].Transpose
            %   parameter = [m, a, b, Cx, Cy].T
            %   u = [pwm, steer].Transpose

            m = p(1);                 %       m  = p(1);   /* Vehicle mass.                    */
            a = p(2);         %       a  = p(2);   /* Distance from front axle to COG. */
            b = p(2);         %       b  = p(3);   /* Distance from rear axle to COG.  */
            Cx = p(3);       %       Cx = p(4);   /* Longitudinal tire stiffness.     */
            Cy = p(4);       %       Cy = p(5);   /* Lateral tire stiffness.          */
            Jz = m*((a+b)/2)^2;
            
            vel_x = x(1,1);             %       /* x(1): Longitudinal vehicle velocity. */
            vel_y = x(2,1);             %       /* x(2): Lateral vehicle velocity. */
            yawRate = x(3,1);           %       /* x(3): Yaw rate. */
            
            slipWheel_fl = u(1);      %       u1(t) = s_FL(t)     Slip of Front Left tire (ratio).
            slipWheel_fr = u(2);      %       u2(t) = s_FR(t)     Slip of Front Right tire (ratio).
            slipWheel_rl = u(3);      %       u3(t) = s_RL(t)     Slip of Rear Left tire (ratio).
            slipWheel_rr = u(4);      %       u4(t) = s_RR(t)     Slip of Rear Right tire (ratio).
            steer = u(5);        %       u5(t) = delta(t)    Steering angle (rad)
            
            slipAngle_f = steer - atan2((vel_y + a*yawRate),vel_x);
            slipAngle_r = -atan2((vel_y - b*yawRate),vel_x);
            %     slipAngle_f = min(0.2,max(-0.2,slipAngle_f));
            %     slipAngle_r = min(0.2,max(-0.2,slipAngle_r));
            slip = [slipAngle_f slipAngle_r];
            
            Fy_f = Cy*slipAngle_f;
            Fy_r = Cy*slipAngle_r;
            Fx_f = Cx*(slipWheel_fl+slipWheel_fr);
            Fx_r = Cx*(slipWheel_rl+slipWheel_rr);
            
            torqueLong = (Fy_f*cos(steer) + Fx_f*sin(steer))*a - Fy_r*b;
            acc_x = (Fx_f*cos(steer) - Fy_f*sin(steer) + Fx_r)/m;
            acc_y = (Fy_f*cos(steer) + Fx_f*sin(steer) + Fy_r)/m;
            cause = [acc_x acc_y torqueLong];
            
            dvX_dt = (acc_x + vel_y*yawRate);
            dvY_dt = (acc_y - vel_x*yawRate);
            dyawRate_dt = (torqueLong/Jz);
            
            xNext = x + [dvX_dt; dvY_dt; dyawRate_dt]*self.dt;
    %dragForces = p(6)*(x(1)^2); CA*(vel_x^2)
  
%   dx(1) = x(2)*x(3)+(1/p(1))*(p(4)*(u(1)+u(2))*cos(u(5))-2*p(5)*(slipAngle_F_adj)*sin(u(5))+p(4)*(u(3)+u(4))-p(6)*(x(1)^2));
%   dx(2) = (-x(1)*x(3))+(1/p(1))*(p(4)*(u(1)+u(2))*sin(u(5))+2*p(5)*(slipAngle_F_adj)*cos(u(5))+2*p(5)*slipAngle_R_adj);
%   dx(3) = (1/((((p(2)+p(3))/2)^2)*p(1)))*(p(2)*(p(4)*(u(1)+u(2))*sin(u(5))+2*p(5)*(slipAngle_F_adj)*cos(u(5)))-2*p(3)*p(5)*slipAngle_R_adj)
%           
            y = xNext;

            
            
            
        end
    
        
        
        
    end
    
end
        