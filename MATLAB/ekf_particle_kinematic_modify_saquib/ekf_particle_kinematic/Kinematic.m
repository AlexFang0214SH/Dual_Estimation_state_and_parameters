classdef Kinematic
    properties
        dt = 0.1;
        nx = 4;
        ny = 4;
        R
        Q
        
    end
    
    methods    
        function [xNext, y] = step(self, x, p, u)
            [xNext, y] = predict(self, x, p, u);
            
            process_noise = self.R * randn(self.nx, 1);
            measurement_noise = self.Q * randn(self.ny, 1);
            
            xNext = xNext + process_noise;
            y = y + measurement_noise;
            
        end
        
        function [G, H] = compute_jacobian(self, dt, x, p, u)

            a = p(1);         %       a  = p(2);   /* Distance from front axle to COG. */
            b = p(2);         %       b  = p(3);   /* Distance from rear axle to COG.  */
            
            acc = u(1);
            steer = u(2);
            
            beta = atan((b*tan(steer))/(a+b));
            
            G = eye(4) + [0 0 -x(4)*sin(x(3)+beta) cos(x(3)+beta); 0 0 x(4)*cos(x(3)+beta) sin(x(3)+beta); 0 0 0 sin(beta)/b; 0 0 0 0]*self.dt;
            H = eye(4);
        end
    
        function [xNext, y] = predict(self, x, p, u)
            %   state = [x,y,phi,v].Transpose
            %   y = [x,y,phi,v].Transpose
            %   parameter = [a, b].T
            %   u = [acc, steer].Transpose


            a = p(1);         %       a  = p(2);   /* Distance from front axle to COG. */
            b = p(2);         %       b  = p(3);   /* Distance from rear axle to COG.  */
            
            acc = u(1);
            steer = u(2);
            
            %params = [lf,lr]
            %states = [x,y,phi,v];
            %u = [acc,steering_angle]
            %equations based on the paper Kong et al, Kinematic and Dynamic Vehicle
            %Models for Autonomous Driving control design
            
            beta = atan((b*tan(steer))/(a+b));
            xNext = x + [x(4)*cos(x(3)+beta); x(4)*sin(x(3)+beta); x(4)*sin(beta)/b; acc]*self.dt;
            
            y = xNext;

            
        end
    
    end
    
end
        