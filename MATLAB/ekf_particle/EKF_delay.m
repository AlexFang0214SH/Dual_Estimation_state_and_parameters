
classdef EKF_delay
    properties
        dt=0.01;
        R
        Q
        Qt
        
        N = 10;
        PARTICLES
        weight
        parameter_range = [4; 100];
        E = [0.01; 0.2];
        evolution_cov = ones(2,1);
        
        C = 30;
        count = 0;
        cov = zeros(2,1);
        
        sigma
        
    end 
    
    methods
        function self = EKF_delay(self)
            self.PARTICLES = rand(2, self.N).*self.parameter_range;
            self.weight = zeros(1, self.N) + 1.0/self.N;
            
            self.sigma = diag([1 1 100*self.dt*pi/180 10.0*self.dt]);
        end
    
        function [G, H] = compute_jacobian(self, dt, x, parameter, u)
            pwm = u(1,1);
            steer = u(2,1);
            
            v = x(4, 1);
            theta = x(3, 1);
            
            
            L_inv = parameter(1,1);
            K_m = parameter(2, 1);
            
            
            G = [   1, 0, -dt*v*sin(theta), dt*cos(theta);
                    0, 1, dt*v*cos(theta), dt*sin(theta);
                    0, 0, 1, dt*L_inv*tan(steer);
                    0, 0, 0, 1];
            
            H = [   1, 0, 0, 0;
                    0, 1, 0, 0;
                    0, 0, 1, 0;
                    0, 0, 0, 1];
        end
    
        function g = gaussian(self, dx, sigma)
            g = exp(-dx'*inv(sigma)*dx/2)/det(sigma)^0.5;
        end
    
        function self= resampling(self)
            %   Best candidate resampling
            [q, i] = max(self.weight);
            winner = self.PARTICLES(:, i);
            self.PARTICLES = repmat(winner, [1, self.N]);
        end
    
        function [xEst, pEst, self] = run(self, cycle, x, parameter, u, z)

            self.count = self.count + 1;
            
            for i=1:self.N
                [xp, zp] = cycle.step(x, self.PARTICLES(:,i), u);
                
                dz = zp - z;
                self.cov = self.cov + dz(3:4);
                
                self.weight(1,i) = self.weight(1,i) + self.gaussian(dz, self.sigma);
            end
            
            w = self.weight/sum(self.weight);
            
            pEst = self.PARTICLES*w';
            
            if mod(self.count, self.C) == 0
                self = self.resampling();
                self.weight(1,i) = 0;
                self.count = 0;
                self.evolution_cov = abs(self.cov/self.C/self.dt);
                self.cov = zeros(2,1);
                
                evolution = randn(2, self.N).*self.E.*self.evolution_cov;
                self.PARTICLES = self.PARTICLES + evolution;
                
            end
            
            
            
            %   EKF with new parameter estimate
            
            [xp, zp] = cycle.step(x, pEst, u);
            
            dz = z - zp;
            
            [G, H] = self.compute_jacobian(cycle.dt, x, pEst, u);
            
            
            Qtp = G*self.Qt* G' + self.R;
            
            
            S = H*Qtp*H' + self.Q;
            S_inv = pinv(S);
            K = Qtp * H' * S_inv ;
            
            xEst = xp + K * dz;
            self.Qt = (eye(4) - K* H)*self.Qt;

        end
        
    end 
    
end 