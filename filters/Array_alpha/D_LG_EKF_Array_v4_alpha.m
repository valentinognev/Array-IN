classdef D_LG_EKF_Array_v4_alpha
    %D_LG_EKF_CLASSIC Summary of this class goes here
    %   Detailed explanation goes here
%  alpha = [omega_dot; s]
    properties
        T                       % Sampling time 
        g                       % gravity in navigation frame
        % Equation (35) Inertial navigationusing an Inertial sensor array
        inds_R = 1:3            % Rotation matrix R_nb, from body to navigation frame
        inds_omega = 4:6        % Angular velocity in body frame
        inds_p = []             % Position in navigation frame
        inds_v = []             % Velocity in navigation frame        
        inds_b_alpha = []       % Bias in alpha (omega_dot)
        % inds_b_s = []         % (Not in use, see Eqn (25)) Bias in specific force
        inds_b_g = []           % Bias in gyroscope 

        % Equation (36) Inertial navigationusing an Inertial sensor array
        inds_w_alpha = []       % alpha (omega_dot) white noise (that drives the bias terms for the accelerometers and gyroscopes)
        % inds_w_s = []         % (Not in use, see Eqn (25)) specific force white noise
        inds_w_b_alpha = []     % Process noise bias alpha
        %inds_w_b_s = []        % Process noise bias specific force        
        inds_w_b_g = []         % Process noise bias gyroscope

        Nx
        Nw
        N_a
        A
        A_omega_dot
        A_s
        r                       % accelerometer sensor position on body coordinates
        T2_R
        constant_b_alpha = []
        Q
        R_pos                   % Covariance matrix for position estimation
        R_rot                   
        R_gyro                  % Covariance matrix for gyro estimation
    end    
    methods
        function obj = D_LG_EKF_Array_v4_alpha(settings)
            %D_LG_EKF_CLASSIC Construct an instance of this class
            %   Detailed explanation goes here
            % Assume one gyro
            
            obj.T = settings.T;
            obj.g = settings.g;
            obj.N_a = size(settings.r,2);

            obj.r = settings.r;
            obj.A = compute_A_non_center(obj.r);
            obj.A_omega_dot = obj.A(1:3,:);
            obj.A_s = obj.A(4:6,:);
            
            if isfield(settings, "set_T2_R_zero") && settings.("set_T2_R_zero")
                obj.T2_R = 0;
            else
                obj.T2_R = obj.T^2;
            end

            % R, omega, defeault
            Nx = 6;                       
            obj.inds_w_alpha = 1:6;  % process noise indexes
            Nw = 6;            
            
            if isfield(settings,"propagate_position") && settings.propagate_position
                obj.inds_p = (1:3) + Nx;
                Nx = Nx + 3;
            end
            if isfield(settings,"propagate_velocity") && settings.propagate_velocity
                obj.inds_v = (1:3) + Nx;
                Nx = Nx + 3;
            end
            
            if isfield(settings,"propagate_bias_alpha") && settings.propagate_bias_alpha
                obj.inds_b_alpha = (1:6) + Nx;
                Nx = Nx + 6;
                                
                obj.inds_w_b_alpha = (1:6) + Nw;
                Nw = Nw + 6;
            end
            
            if isfield(settings,"propagate_bias_gyro") && settings.propagate_bias_gyro
                obj.inds_b_g = (1:3) + Nx;
                Nx = Nx + 3;
                
                obj.inds_w_b_g = (1:3) + Nw;                
                Nw = Nw + 3;
            end
            
            obj.Nx = Nx;
            obj.Nw = Nw;

            % -----------------------------------------------------------------
            % Other constants
            % -----------------------------------------------------------------
           
            if isfield(settings,"constant_b_alpha")
                obj.constant_b_alpha = settings.constant_b_alpha;
            end
            
            if isfield(settings,"constant_b_alpha") && ~isempty(obj.inds_b_alpha)
                error("Cannot have constant b_alpha and have in the state vector")
            end
            
            % -------------------------------------------------------------
            % Fill Q
            % -------------------------------------------------------------
            Q = zeros(Nw,Nw);
            if isfield(settings,"Q_alpha")
                Q(obj.inds_w_alpha, obj.inds_w_alpha) = settings.Q_alpha;
            elseif isfield(settings,"Q_acc")
                Q(obj.inds_w_alpha, obj.inds_w_alpha) = obj.A*settings.Q_acc*obj.A';
            else
                error("No covariance for alpha")
            end

            if ~isempty(obj.inds_w_b_alpha) 
                if isfield(settings,"Q_bias_alpha")
                    Q(obj.inds_w_b_alpha,obj.inds_w_b_alpha) = settings.Q_bias_alpha;
                elseif isfield(settings,"Q_bias_acc")
                    Q(obj.inds_w_b_alpha,obj.inds_w_b_alpha) = obj.A*settings.Q_bias_acc*obj.A';
                else
                    error("No covariance for bias alpha")
                end
            end

            if ~isempty(obj.inds_w_b_g)
                Q(obj.inds_w_b_g, obj.inds_w_b_g) = settings.Q_bias_gyro;
            end   
            obj.Q = Q;
            
            % -----------------------------------------------------------------
            % Update covariance
            % -----------------------------------------------------------------
           
            if isfield(settings,"R_pos")
                obj.R_pos = settings.R_pos;
                assert(issymmetric(obj.R_pos))
                assert(all(size(obj.R_pos) == [3,3]))
            end
            if isfield(settings,"R_gyro")
                obj.R_gyro = settings.R_gyro;
                assert(issymmetric(obj.R_gyro))
                assert(all(size(obj.R_gyro) == [3,3]))
            end
            if isfield(settings,"R_rot")
                obj.R_rot = settings.R_rot;
                assert(issymmetric(obj.R_rot))
                assert(all(size(obj.R_rot) == [3,3]))
            end
        end        
        function u = get_input(obj, sensorData)
            assert(size(sensorData.acc_measurements,1) == 3*obj.N_a);
            u = sensorData.acc_measurements;
        end
        function Q = get_Q(obj)
            Q = obj.Q;         
        end
        function res = propagate(obj, Rnb, x, y, w)
            % Rnb - Lie group rotation matrix - Rotation between body frame to navigation frame
            % x - position
            %     x(1:3);     % angular velocity in body frame
            %     x(4:6);     % position in navigation frame
            %     x(7:9);     % velocity in navigation frame
            %     x(10:12)    % bias in angular acceleration measurement
            %     x(13:15)    % bias in linear acceleration measurement
            %     x(16:18)    % bias in gyro measurement

            % y - sensor acceleration data in body frame

            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
                            
            assert(length(x) == obj.Nx - 3);
            assert(length(w) == obj.Nw);

            
            omega = x(obj.inds_omega-3);     % angular velocity
            
            if ~isempty(obj.inds_v)
                v = x(obj.inds_v-3);             % velocity
            end
            
            if ~isempty(obj.inds_b_alpha)            % alpha is omega dot
                b_alpha = x(obj.inds_b_alpha - 3);   % bias of angular acceleration
            elseif ~isempty(obj.constant_b_alpha)
                b_alpha = obj.constant_b_alpha;
            else
                b_alpha = zeros(6,1);
            end
            
            % Process noise
            w_alpha = w(obj.inds_w_alpha);     
            
            if ~isempty(obj.inds_w_b_alpha)
                w_b_alpha = w(obj.inds_w_b_alpha);     % omega dot(angular acceleration) white noise bias process 
            else
                w_b_alpha = zeros(6,1);   
            end
                        
            if ~isempty(obj.inds_w_b_g)
                w_b_g = w(obj.inds_w_b_g);   % gyro white noise bias process
            else
                w_b_g = zeros(3,1);   
            end            
                        
            if 0
                y_a = y;
                h = reshape(HatSO3(omega)^2*obj.r,[],1); % h = h(omega, T_a,r)
                alpha = obj.A*(y_a - h) + b_alpha + w_alpha;
            else
                omegaPrev = omega;
                [bdotdot, omega_dot_b, omegaNew]=cardou12(y, omega, obj.r, obj.T);
                alpha(1:3,1) = omega_dot_b;
                alpha(4:6,1) = bdotdot;
                alpha = alpha + b_alpha + w_alpha;
            end
            omega_dot = alpha(1:3); % Angular acceleration in body coordinates           
            s = alpha(4:6);         % Specific force (acceleration) in body coordinates       

            % Navigation acceleration
            v_dot = obj.g + Rnb*s;    % Acceleration in navigation coordinates 
  
            Omega = zeros(obj.Nx,1);
            % R
            Omega(obj.inds_R) = omega*obj.T + omega_dot*obj.T2_R/2;
            % omega            
            Omega(obj.inds_omega) = omega_dot*obj.T;

            % p
            if ~isempty(obj.inds_p)
                Omega(obj.inds_p) = v*obj.T + v_dot*obj.T^2/2;
            end
            % v
            if ~isempty(obj.inds_v)
                Omega(obj.inds_v) = v_dot*obj.T;
            end
            
            if ~isempty(obj.inds_b_g)
                Omega(obj.inds_b_g) = w_b_g;
            end
            
            if ~isempty(obj.inds_b_alpha)
                Omega(obj.inds_b_alpha) = w_b_alpha;
            end
                        
            % -----------------------------------------------------------------
            % Jacobian in x            
            d_h_d_omega = compute_d_h_d_omega(obj, omega);
            
            d_omega_dot_d_omega = -obj.A_omega_dot*d_h_d_omega;
            
            d_v_dot_d_R = -Rnb*HatSO3(s);
            d_v_dot_d_omega = -Rnb*obj.A_s*d_h_d_omega;
            
            
            % Fill Jacobian
            dOmega_de = zeros(obj.Nx, obj.Nx);
            
            % -----------------------------------------------------------------
            % R equation
            % Nothing in R
            
            dOmega_de(obj.inds_R,obj.inds_omega) = eye(3)*obj.T + d_omega_dot_d_omega*obj.T2_R/2; % omega
            
            % nothing in p
            % nothing in v
            if ~isempty(obj.inds_b_alpha)
                d_omega_dot_d_b_alpha = [eye(3) zeros(3)];
                dOmega_de(obj.inds_R,obj.inds_b_alpha) = d_omega_dot_d_b_alpha*obj.T2_R/2;
            end
            
            % -----------------------------------------------------------------
            % w equation
            % nothing in R            
            dOmega_de(obj.inds_omega,obj.inds_omega) = d_omega_dot_d_omega*obj.T; % omega                
                                                                                  % nothing in p
                                                                                  % nothing in v
            if ~isempty(obj.inds_b_alpha) 
                dOmega_de(obj.inds_omega,obj.inds_b_alpha) = d_omega_dot_d_b_alpha*obj.T;
            end
                            
            % -----------------------------------------------------------------
            % p equation
            if ~isempty(obj.inds_p)
                dOmega_de(obj.inds_p,obj.inds_R) = d_v_dot_d_R*obj.T^2/2;         % R
                dOmega_de(obj.inds_p,obj.inds_omega) = d_v_dot_d_omega*obj.T^2/2; % omega
                
                if ~isempty(obj.inds_v)
                    dOmega_de(obj.inds_p,obj.inds_v) = eye(3)*obj.T;                      % v
                end
                if ~isempty(obj.inds_b_alpha)
                    d_v_dot_d_b_alpha = [zeros(3) Rnb];
                    dOmega_de(obj.inds_p,obj.inds_b_alpha) = d_v_dot_d_b_alpha*obj.T^2/2;            % b_alpha
                end
            end

            % -----------------------------------------------------------------
            % v equation
            if ~isempty(obj.inds_v)
                dOmega_de(obj.inds_v, obj.inds_R) = d_v_dot_d_R*obj.T;          % R
                dOmega_de(obj.inds_v, obj.inds_omega) = d_v_dot_d_omega*obj.T;  % omega
                
                if ~isempty(obj.inds_b_alpha)
                    dOmega_de(obj.inds_v,obj.inds_b_alpha) = d_v_dot_d_b_alpha*obj.T;         % b_alpha
                end
            end
            
            % -----------------------------------------------------------------
            % Jacobian in w            

            dOmega_dw = zeros(obj.Nx, obj.Nw);
                        
            d_omega_dot_d_w_alpha = [eye(3) zeros(3)];
            d_v_dot_d_w_alpha = [zeros(3) Rnb];
            
            % -----------------------------------------------------------------------------
            % R equation
            dOmega_dw(obj.inds_R, obj.inds_w_alpha) = d_omega_dot_d_w_alpha*obj.T2_R/2;% w_alpha
            % nothing in w_b_alpha
            % nothing in w_b_g
            
            % -----------------------------------------------------------------------------
            % omega equation            
            dOmega_dw(obj.inds_omega, obj.inds_w_alpha) = d_omega_dot_d_w_alpha*obj.T;  % w_alpha
            % nothing in w_b_alpha
            % nothing in w_b_g
            
            % -----------------------------------------------------------------------------
            % p equation
            if ~isempty(obj.inds_p)
                dOmega_dw(obj.inds_p, obj.inds_w_alpha) = d_v_dot_d_w_alpha*obj.T^2/2;  % w_alpha
            end
            % nothing in w_b_alpha
            % nothing in w_b_g
            
            % -----------------------------------------------------------------------------
            % v equation
            if ~isempty(obj.inds_v)
                dOmega_dw(obj.inds_v,obj.inds_w_alpha) = d_v_dot_d_w_alpha*obj.T;      % w_alpha
            end
            % nothing in w_b_alpha
            % nothing in w_b_g
            
            % -----------------------------------------------------------------------------
            % b_s equation
            
            if ~isempty(obj.inds_b_alpha)
                dOmega_dw(obj.inds_b_alpha, obj.inds_w_b_alpha) = eye(6);
            end
            
            if ~isempty(obj.inds_b_g)
                dOmega_dw(obj.inds_b_g,obj.inds_w_b_g) = eye(3);
            end            
            
            res = struct;
            res.Omega = Omega;
            res.dOmega_de = dOmega_de;
            res.dOmega_dw = dOmega_dw;
            res.v_dot = v_dot;
            res.omega_dot = omega_dot;
            res.s = s;
        end
        function d_h_d_omega = compute_d_h_d_omega(obj, w)

            row1 = 1:3:3*obj.N_a;
            row2 = 2:3:3*obj.N_a;
            row3 = 3:3:3*obj.N_a;
            d_h_d_omega = zeros(3*obj.N_a,3);
            
            r1 = obj.r(1,:);
            r2 = obj.r(2,:);
            r3 = obj.r(3,:);
            
            r1w1 = w(1).*r1;
            r1w2 = r1.*w(2);
            r1w3 = r1.*w(3);
            
            r2w1 = r2.*w(1);
            r2w2 = w(2).*r2;
            r2w3 = r2.*w(3);
            
            r3w1 = r3.*w(1);
            r3w2 = r3.*w(2);
            r3w3 = w(3).*r3;
            
            d_h_d_omega(row1,1) = r2w2 + r3w3;
            d_h_d_omega(row2,1) = r1w2 - 2*r2w1;
            d_h_d_omega(row3,1) = r1w3 - 2*r3w1;
            
            d_h_d_omega(row1,2) = r2w1 - 2*r1w2;
            d_h_d_omega(row2,2) = r1w1 + r3w3;
            d_h_d_omega(row3,2) = r2w3 - 2*r3w2;
            
            d_h_d_omega(row1,3) = r3w1 - 2*r1w3;
            d_h_d_omega(row2,3) = r3w2 - 2*r2w3;
            d_h_d_omega(row3,3) = r1w1 + r2w2;
            
        end
        function [e,H,Q] = position_update(obj, p_obs, ~, x_in)
            
            assert(~isempty(obj.R_pos))
            Q = obj.R_pos;
            
            H = zeros(3, obj.Nx);
            H(:,obj.inds_p) = eye(3); % p
                        
            
            p_pred = x_in(obj.inds_p-3);
            e =  p_obs - p_pred;
        end
        function [e, H, Q] = gyroscope_update(obj, u_g, ~, x_in)
            
            assert(~isempty(obj.R_gyro))
            Q = obj.R_gyro;
            
            H = zeros(3, obj.Nx);
            H(:,obj.inds_omega) = eye(3);
            
            omega = x_in(obj.inds_omega - 3);

            if ~isempty(obj.inds_b_g)
                b_g = x_in(obj.inds_b_g - 3);
                H(:,obj.inds_b_g) = eye(3);
            else
                b_g = zeros(3,1);                
            end

            e =  u_g - omega - b_g;
        end
        function [e,H,Q] = rotation_update(obj, R_obs, R_pred, ~)
       
            assert(~isempty(obj.R_rot))
            Q = obj.R_rot;
            
            % H rotation
            H = zeros(3, obj.Nx);
            H(:,1:3) = eye(3); % R    

            e = logSO3(invSO3(R_pred)*R_obs);
        end
        function initOut = get_initial_conditions(obj, initIn)
            
            
            initOut = struct;
            use_full = (isfield(initIn,"x") || isfield(initIn,"P") || isfield(initIn,"R"));
            use_partial = (isfield(initIn,"mean") || isfield(initIn,"cov"));
            if use_full && use_partial
                error("Both (R0,x0,P0) and (mean,cov) defined, use only one")
            elseif isfield(initIn,"x") && isfield(initIn,"P")
                initOut.R0 = initIn.R;
                initOut.x0 = initIn.x;
                initOut.P0 = initIn.P;
            elseif isfield(initIn,"mean") && isfield(initIn,"cov")
                m = initIn.mean;
                
                assert(all(size(m.R) == [3,3]));
                initOut.R0 = m.R;
                x0 = zeros(obj.Nx-3,1);
                x0(obj.inds_omega - 3) = m.omega;
                
                if ~isempty(obj.inds_p)
                    x0(obj.inds_p - 3) = m.p;
                end
                
                if ~isempty(obj.inds_v)
                    x0(obj.inds_v - 3) = m.v;
                end
                                
                if ~isempty(obj.inds_b_alpha)
                    if isfield(m,"b_alpha")
                        x0(obj.inds_b_alpha - 3) = m.b_alpha;
                    elseif isfield(m,"b_a")
                        x0(obj.inds_b_alpha - 3) = -obj.A*m.b_a;
                    else
                        error("No mean value alpha bias")
                    end
                end

                initOut.x0 = x0;
                
                c = initIn.cov;
                P0 = zeros(obj.Nx, obj.Nx);
                P0(obj.inds_R, obj.inds_R) = c.R;
                P0(obj.inds_omega, obj.inds_omega) = c.omega;
                
                if ~isempty(obj.inds_p)
                    P0(obj.inds_p, obj.inds_p) = c.p;
                end
                if ~isempty(obj.inds_v)
                    P0(obj.inds_v, obj.inds_v) = c.v;
                end
                if ~isempty(obj.inds_b_g)
                    P0(obj.inds_b_g, obj.inds_b_g) = c.b_g;
                end
                if ~isempty(obj.inds_b_alpha)
                    if isfield(m,"b_alpha")
                        P0(obj.inds_b_alpha,obj.inds_b_alpha) = c.b_alpha;
                    elseif isfield(m,"b_a")
                        P0(obj.inds_b_alpha,obj.inds_b_alpha) = obj.A*c.b_a*obj.A';
                    else
                        error("No cov value alpha bias")
                    end
                end
                
                initOut.P0 = P0;
            else
                error("No initial conditions")
            end
        end
        function Sout = extract_variables(obj,Sin)
            Sout = struct;
            Sout.mean.R = Sin.R;
            Sout.mean.omega = Sin.x(obj.inds_omega - 3, :);
            
            if ~isempty(obj.inds_p)
                Sout.mean.p = Sin.x(obj.inds_p - 3, :);
            end
            if ~isempty(obj.inds_v)
                Sout.mean.v = Sin.x(obj.inds_v - 3, :);
            end                        
            if ~isempty(obj.inds_b_g)
                Sout.mean.b_g = Sin.x(obj.inds_b_g - 3, :);
            end
            if ~isempty(obj.inds_b_alpha)
                Sout.mean.b_alpha = Sin.x(obj.inds_b_alpha - 3, :);
                Sout.mean.b_omega_dot = Sout.mean.b_alpha(1:3,:);
                Sout.mean.b_s = Sout.mean.b_alpha(4:6,:);
            end            
            
            Sout.std.R = Sin.std(obj.inds_R, :);
            Sout.std.omega = Sin.std(obj.inds_omega, :);
            
            if ~isempty(obj.inds_p)
                Sout.std.p = Sin.std(obj.inds_p, :);
            end
            if ~isempty(obj.inds_v)
                Sout.std.v = Sin.std(obj.inds_v, :);
            end
            
            if ~isempty(obj.inds_b_g)
                Sout.std.b_g = Sin.std(obj.inds_b_g, :);
            end
            if ~isempty(obj.inds_b_alpha)
                Sout.std.b_alpha = Sin.std(obj.inds_b_alpha, :);
                Sout.std.b_omega_dot = Sout.std.b_alpha(1:3,:);
                Sout.std.b_s = Sout.std.b_alpha(4:6,:);
            end            
                     
        end
        function print_info(obj, S_init, SensorData, settings)
            
            fprintf("Sampling time: %.2e [s]\n", settings.T)
            fprintf("Sampling freq: %.1f [Hz]\n", 1/settings.T)
            fprintf("T2 for rotation: %.1e [s]\n", obj.T2_R)
            fprintf("Propagate rotation: %d\n", ~isempty(obj.inds_R))
            fprintf("Propagate omega: %d\n", ~isempty(obj.inds_omega))            
            fprintf("Propagate position: %d\n", ~isempty(obj.inds_p))
            fprintf("Propagate velocity: %d\n", ~isempty(obj.inds_v))
            fprintf("Propagate bias alpha: %d\n", ~isempty(obj.inds_b_alpha))
            fprintf("Propagate bias gyro: %d\n", ~isempty(obj.inds_b_g))

            %------------------------------------------------------------------            
            fprintf("Propagation data:\n")
            fprintf("\tAccelerometer data:\n")
            try
                Q_sqrt = chol(obj.Q(obj.inds_w_alpha, obj.inds_w_alpha));
                fprintf("\t\twhite noise chol(Q_alpha) [mixed units]\n")
                fprintf("\t\t%10.1e%10.1e%10.1e%10.1e%10.1e%10.1e\n", Q_sqrt')
            catch
                Q_alpha = (obj.Q(obj.inds_w_alpha, obj.inds_w_alpha));
                fprintf("\t\twhite noise Q_alpha [mixed units]^2\n")
                fprintf("\t\t%10.1e%10.1e%10.1e%10.1e%10.1e%10.1e\n", Q_alpha')
            end
            
            if ~isempty(obj.inds_b_alpha)
                try
                    Q_sqrt = chol(obj.Q(obj.inds_w_b_alpha, obj.inds_w_b_alpha));
                    fprintf("\t\tchol(Q) bias alpha noise: [mixed units]\n")
                    fprintf("\t\t%10.1e%10.1e%10.1e%10.1e%10.1e%10.1e\n", Q_sqrt')
                    fprintf("\t\t[mixed units]\n")
                catch
                    fprintf("\t\tQ bias alpha noise: [mixed units]\n")
                    fprintf("\t\t%10.1e%10.1e%10.1e%10.1e%10.1e%10.1e\n", obj.Q(obj.inds_w_b_alpha, obj.inds_w_b_alpha)')
                    fprintf("\t\t[mixed units]\n")
                end
            end
            
            
            if isfield(settings,"constant_b_alpha")
                fprintf("\t\tConstant b_alpha:          [")
                fprintf("%.1e ", settings.constant_b_alpha)
                fprintf("] []\n")
            end
            
            fprintf("\t\tStart of accelerometer triad measurements:\n")
            fprintf("\t\t%10.1f %10.1f %10.1f\n", SensorData.acc_measurements(1:3,1:3)')

            fprintf("\t\tStart of accelerometer positions:\n")
            fprintf("\t\t%10.1e %10.1e %10.1e\n", obj.r(1:3,1:3)')
            
            fprintf("\t\tMean accelerometer positions: [")
            fprintf("%.1e ", mean(obj.r,2))
            fprintf("] [m]\n")

            fprintf("\tGyroscope data:\n")
            if ~isempty(obj.inds_b_g)
                try
                    Q_sqrt = rad2deg(chol(obj.Q(obj.inds_w_b_g, obj.inds_w_b_g)));
                    fprintf("\t\tchol(Q) bias gyro noise: [deg/s]\n")
                    fprintf("\t\t%10.1e%10.1e%10.1e\n", Q_sqrt')
                    fprintf("\t\t[deg/s]\n")
                catch
                    fprintf("\t\tQ bias gyro noise: [deg/s]^2\n")
                    fprintf("\t\t%10.1e%10.1e%10.1e\n", rad2deg(rad2deg(obj.Q(obj.inds_w_b_g, obj.inds_w_b_g)))')
                    fprintf("\t\t[deg/s]^2\n")
                end
            end
 
            fprintf("Initial conditions:\n")
            fprintf("\tMean:\n")
            e_deg = rad2deg(my_rotm2eul(S_init.mean.R));
            fprintf("\t\tRotation: roll: %.1f, pitch: %.1f, yaw: %.1f [deg]\n", e_deg(1), e_deg(2), e_deg(3) )
            
            if ~isempty(obj.inds_omega)
                fprintf("\t\tAngular Velocity:  [")
                fprintf("%.1f ", rad2deg((S_init.mean.omega)))
                fprintf("] [deg/s]\n")
            end
            if ~isempty(obj.inds_p)
                fprintf("\t\tPosition:          [")
                fprintf("%.1f ", (S_init.mean.p))
                fprintf("] [m]\n")
            end
            
            if ~isempty(obj.inds_v)
                fprintf("\t\tVelocity:          [")
                fprintf("%.1f ", (S_init.mean.v))
                fprintf("] [m/s]\n")
            end
                       
            if ~isempty(obj.inds_b_alpha)
                fprintf("\t\tBias omega dot:    [")
                fprintf("%.1f ", rad2deg((S_init.mean.b_alpha(1:3))))
                fprintf("] [deg/s^2]\n")
                fprintf("\t\tBias s:        [")
                fprintf("%.1f ", (S_init.mean.b_alpha(4:6)))
                fprintf("] [m/s^2]\n")
            end
            
            fprintf("\tStd:\n")
            
            fprintf("\t\tRotation:          [")
            fprintf("%.1f ", rad2deg(S_init.std.R))
            fprintf("] [deg]\n")
            
            if ~isempty(obj.inds_omega)
                fprintf("\t\tAngular Velocity:  [")
                fprintf("%.1f ", rad2deg(S_init.std.omega))
                fprintf("] [deg/s]\n")
            end
            
            if ~isempty(obj.inds_p)
                fprintf("\t\tPosition:          [")
                fprintf("%.1f ", (S_init.std.p))
                fprintf("] [m]\n")
            end
            if ~isempty(obj.inds_v)
                fprintf("\t\tVelocity:          [")
                fprintf("%.1f ", (S_init.std.v))
                fprintf("] [m/s]\n")
            end
                                    
            if ~isempty(obj.inds_b_alpha)
                fprintf("\t\tBias omega dot:    [")
                fprintf("%.1f ", rad2deg((S_init.std.b_alpha(1:3))))
                fprintf("] [deg/s^2]\n")
                fprintf("\t\tBias s:        [")
                fprintf("%.1f ", (S_init.std.b_alpha(4:6)))
                fprintf("] [m/s^2]\n")
            end
        end
    end
end

