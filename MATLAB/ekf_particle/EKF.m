
classdef EKF
    properties
        dt=0.1;
        R
        Q
        Qt
        
        N = 500;
        PARTICLES
        weight
        parameter_range = [6; 3000];
        E = [0.02*3; 0.01*1700];
        evolution_cov = ones(2,1);
        
        sigma
        
    end 
    
    methods
        function self = EKF(self)
            self.PARTICLES = rand(2, self.N).*self.parameter_range;
            self.weight = zeros(1, self.N);
            
            self.sigma = diag([1 1 100*self.dt*pi/180 10.0*self.dt])*200;
        end
    
        function [G, H] = compute_jacobian(self, dt, x, parameter, u)
            pwm = u(1,1);
            steer = u(2,1);
            
            v = x(4, 1);
            theta = x(3, 1);
            
            
            L_inv = 1.0/parameter(1,1);
            K_m = 1.0/parameter(2, 1);
            
            
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
%             [q, i] = max(self.weight);
%             winner = self.PARTICLES(:, i);
%             self.PARTICLES = repmat(winner, [1, self.N]);

            % Preserve top 10
            [w, idx] = sort(self.weight, 'descend');
            top = self.PARTICLES(:, idx(1,1:10));
            self.PARTICLES = repmat(top, [1, self.N/10]);
            
%             diffusion = randn(4, self.N).*self.E*self.diffusion_mag;
%             diffusion(:, 1:200) = zeros(4, 200);
%             self.PARTICLES = self.PARTICLES + diffusion;
        end
    
        function [xEst, pEst, self] = run(self, cycle, x, parameter, u, z)

            evolution = randn(2, self.N).*self.E.*self.evolution_cov;
            self.PARTICLES = self.PARTICLES + evolution;
            
            cov = zeros(2,1);
            self.weight = zeros(1, self.N);
            for i=1:self.N
                [xp, zp] = cycle.predict(x, self.PARTICLES(:,i), u);
                
                dz = zp - z;
                cov = cov + dz(3:4);
                
                self.weight(1,i) = self.weight(1,i) + self.gaussian(dz, self.sigma);
            end

            self.evolution_cov = self.evolution_cov*0.997 + 0.003*abs(cov/self.N/self.dt);
            self.weight = self.weight/sum(self.weight);
            
            
            self = self.resampling();
            
            pEst = mean(self.PARTICLES(:, 1:3), 2);
            
            
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