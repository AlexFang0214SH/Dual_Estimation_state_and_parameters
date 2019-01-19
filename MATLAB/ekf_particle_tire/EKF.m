
classdef EKF
    properties
        dt=0.1;
        R
        Q
        Qt
        
        count = 0;
        
        N = 1000;
        PARTICLES
        weight
        parameter_range = [1700.; 1.; 2e5; 8e4];
        lower_range = [1000.0;1.0; 1e4; 1e4];
        E = [0.; 0.0; 0.1e5; 0.01e5];
        diffusion_mag = 1;
        
        C = 1;
        gamma = 0.9;
        
        sigma
        alpha = 0.998;
        
    end 
    
    methods
        function self = EKF(self)
            self.PARTICLES = rand(4, self.N).*self.parameter_range + repmat(self.lower_range, [1, self.N]);
            self.PARTICLES(1:2, :) = repmat([1700; 1.5], [1, self.N]);
            self.weight = zeros(1, self.N);
            
            self.sigma = diag([1 1 1])*500;
        end
    
        function [G, H] = compute_jacobian(self, dt, x, p, u)
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
            steer = u(5);               %       u5(t) = delta(t)    Steering angle (rad)
            
            
            G(1,:) = [                                                                                                                                          1 - (dt*Cy*sin(steer)*(vel_y + a*yawRate))/(m*((vel_y + a*yawRate)^2 + vel_x^2)),                                                                                                           dt*(yawRate + (Cy*vel_x*sin(steer))/(m*((vel_y + a*yawRate)^2 + vel_x^2))),                                                                                                           dt*(vel_y + (a*Cy*vel_x*sin(steer))/(m*((vel_y + a*yawRate)^2 + vel_x^2)))];
            G(2,:) = [                                           -dt*(yawRate - ((Cy*(vel_y - b*yawRate))/((vel_y - b*yawRate)^2 + vel_x^2) + (Cy*cos(steer)*(vel_y + a*yawRate))/((vel_y + a*yawRate)^2 + vel_x^2))/m),                                                 1 - (dt*((Cy*vel_x)/((vel_y - b*yawRate)^2 + vel_x^2) + (Cy*vel_x*cos(steer))/((vel_y + a*yawRate)^2 + vel_x^2)))/m,                                -dt*(vel_x - ((b*Cy*vel_x)/((vel_y - b*yawRate)^2 + vel_x^2) - (a*Cy*vel_x*cos(steer))/((vel_y + a*yawRate)^2 + vel_x^2))/m)];
            G(3,:) = [ -(dt*((b*Cy*(vel_y - b*yawRate))/((vel_y - b*yawRate)^2 + vel_x^2) - (a*Cy*cos(steer)*(vel_y + a*yawRate))/((vel_y + a*yawRate)^2 + vel_x^2)))/(m*(a/2 + b/2)^2), (dt*((b*Cy*vel_x)/((vel_y - b*yawRate)^2 + vel_x^2) - (a*Cy*vel_x*cos(steer))/((vel_y + a*yawRate)^2 + vel_x^2)))/(m*(a/2 + b/2)^2), 1 - (dt*((b^2*Cy*vel_x)/((vel_y - b*yawRate)^2 + vel_x^2) + (a^2*Cy*vel_x*cos(steer))/((vel_y + a*yawRate)^2 + vel_x^2)))/(m*(a/2 + b/2)^2)];

    
            H = eye(3);
        end
    
        function g = gaussian(self, dx, sigma)
            g = exp(-dx'*inv(sigma)*dx/2)/det(sigma)^0.5;
        end
    
        function self= resampling(self)
            
            [w, idx] = sort(self.weight, 'descend');
            
            % Preserve top 10
            top = self.PARTICLES(:, idx(1,1:10));
            self.PARTICLES = repmat(top, [1, self.N/10]);
               
            
%             [q, i] = max(self.weight(1, 1:self.N));
%             winner = self.PARTICLES(:, i);
%             self.PARTICLES = repmat(winner, [1, self.N]);
            
            
            %   Best candidate resampling
%             [q, i] = max(self.weight(1, 1:self.N/2));
%             winner = self.PARTICLES(:, i);
%             self.PARTICLES(:,  1:self.N/2) = repmat(winner, [1, self.N/2]);
%                 
%             [q, i] = max(self.weight(1, self.N/2: self.N));
%             winner = self.PARTICLES(:, i);
%             self.PARTICLES(:,  self.N/2+1: self.N) = repmat(winner, [1, self.N/2]);
%        
            diffusion = randn(4, self.N).*self.E*self.diffusion_mag;
            diffusion(:, 1:200) = zeros(4, 200);
            self.PARTICLES = self.PARTICLES + diffusion;
            self.weight = repmat(w(1:10), [1, self.N/10])*self.gamma*self.C;
        
        end
    
        function [xEst, pEst, self] = run(self, cycle, x, parameter, u, z)
            self.count = self.count + 1;

            Z = zeros(3, self.N);
            
            diffusion = 0;
            for i=1:self.N
                [xp, zp] = cycle.step(x, self.PARTICLES(:,i), u);
                
                dz = zp - z;
                Z(:, i) = dz;
                diffusion = diffusion + norm(dz);
                
                self.weight(1,i) = self.weight(1,i) + self.gaussian(dz, self.sigma);
            end
            self.diffusion_mag = self.diffusion_mag*self.alpha + (1 - self.alpha)*1*diffusion/self.C/self.N;
            self.diffusion_mag
            
            if self.count == 1
                scatter(Z(1, :), Z(2,:));
                hold on
            end
            if self.count == 50
                scatter(Z(1, :), Z(2,:), [], 'd');
                hold on
            end
            if self.count == 150
                scatter(Z(1, :), Z(2,:), [], 'x');
                hold on
            end
            if self.count == 250
                scatter(Z(1, :), Z(2,:), [], 's');
                hold on
                scatter(0,0, [], 'r', "filled")
            end
            
            
            self.weight = self.weight/sum(self.weight);
            
            
            
            self = self.resampling();
            pEst = mean(self.PARTICLES(:, 1:3), 2);
            
            %pEst = self.PARTICLES(:, 1);
            
            
            %   EKF with new parameter estimate
            
            [xp, zp] = cycle.predict(x, pEst, u);
            
            dz = z - zp;
            
            [G, H] = self.compute_jacobian(cycle.dt, x, pEst, u);
            
            Qtp = G*self.Qt* G' + self.R;
            
            S = H*Qtp*H' + self.Q;
            S_inv = pinv(S);
            K = Qtp * H' * S_inv ;
            
            xEst = xp + K * dz;
            self.Qt = (eye(3) - K* H)*self.Qt;

        end
        
    end 
    
end 