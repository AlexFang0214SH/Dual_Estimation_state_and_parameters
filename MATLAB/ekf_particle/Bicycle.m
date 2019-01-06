classdef Bicycle
    properties
        dt = 0.01;
        nx = 4;
        ny = 4;
        R
        Q
        
        C = [   1 0 0 0;
                0 1 0 0;
                0 0 1 0;
                0 0 0 1; ];
            
        D = [   0 0 0 0;
                0 0 0 0;
                0 0 0 0;
                0 0 0 0; ];
    end
    
    methods    
        function [xNext, y] = predict(self, x, parameter, u)
            %   state = [x, y, theta, v].Transpose
            %   y = [x, y, theta, v].Transpose
            %   parameter = [L, m].T
            %   u = [pwm, steer].Transpose
            
            pwm = u(1,1);
            steer = u(2,1);
            
            theta = x(3,1);
            v = x(4,1);
            L_inv = 1.0/parameter(1,1);
            k_by_m = 1.0/parameter(2,1);
            
            H = [v*cos(theta); v*sin(theta); v*L_inv*tan(steer); k_by_m*pwm];
            
            
            xNext = x + H*self.dt ;
            
            y = self.C * xNext + self.D*H;
        end
    
        function [xNext, y] = step(self, x, parameter, u)
            [xNext, y] = self.predict(x, parameter, u);
            
            noise_process = self.R * randn(self.nx, 1);
            noise_measurement = self.Q * randn(self.ny, 1);
            
            xNext = xNext + noise_process ;
            
            y = y + noise_measurement;
        end
    
        
        
        
    end
    
end
        