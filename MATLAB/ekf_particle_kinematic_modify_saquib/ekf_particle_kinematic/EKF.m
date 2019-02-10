
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
        parameter_range = [1.0; 1.0];
        lower_range = [1.0; 1.0;];
        E = [0.1; 0.1];
        diffusion_mag = 1;
        
        C = 1;
        gamma = 0.0;
        
        sigma
        alpha = 0.98;
        
    end 
    
    methods
        function self = EKF(self)
            self.PARTICLES = rand(2, self.N).*self.parameter_range + repmat(self.lower_range, [1, self.N]);
            %self.PARTICLES(1:3, :) = repmat([1700;1.5; 1.5;], [1, self.N]);
            self.weight = zeros(1, self.N);
            
            self.sigma = diag([1 1 1 1])*10;
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
            diffusion = randn(2, self.N).*self.E*self.diffusion_mag;
            diffusion(:, 1:200) = zeros(2, 200);
            self.PARTICLES = self.PARTICLES + diffusion;
            self.weight = repmat(w(1:10), [1, self.N/10])*self.gamma*self.C;
        
        end
    
        function [xEst, pEst, self] = run(self, cycle, x, parameter, u, z)
            self.count = self.count + 1;

            Z = zeros(4, self.N);
            
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
            
%             if self.count == 1
%                 scatter(Z(1, :), Z(2,:));
%                 hold on
%             end
%             if self.count == 50
%                 scatter(Z(1, :), Z(2,:), [], 'd');
%                 hold on
%             end
%             if self.count == 150
%                 scatter(Z(1, :), Z(2,:), [], 'x');
%                 hold on
%             end
%             if self.count == 250
%                 scatter(Z(1, :), Z(2,:), [], 's');
%                 hold on
%                 scatter(0,0, [], 'r', "filled")
%             end
            
            
            self.weight = self.weight/sum(self.weight);
            
            
            
            
            pEst = self.PARTICLES*self.weight'
            
            self = self.resampling();
            %pEst = mean(self.PARTICLES(:, 1:3), 2);
            
            %pEst = self.PARTICLES(:, 1);
            
            
            %   EKF with new parameter estimate
            
            [xp, zp] = cycle.predict(x, pEst, u);
            
            dz = z - zp;
            
            [G, H] = cycle.compute_jacobian(cycle.dt, x, pEst, u);
            
            Qtp = G*self.Qt* G' + self.R;
            
            S = H*Qtp*H' + self.Q;
            S_inv = pinv(S);
            K = Qtp * H' * S_inv ;
            
            xEst = xp + K * dz;
            self.Qt = (eye(4) - K* H)*self.Qt;

        end
        
    end 
    
end 