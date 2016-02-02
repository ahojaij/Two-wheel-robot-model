% Extended Kalman filter example
clear;clc;

% Discrete time step
dt = 0.2;

% Initial State
x0 = [0 0 0 0 0 0]';

% Prior
mu = [22 -1.8 3.5]'; % mean (mu)
S = 1*eye(6);% covariance (Sigma)

% Discrete motion model
% Ad = [ 0 0 -(r*w1+r*w2)/2*sin(x3)*dt ; 
%     0 0 (r*w1+r*w2)/2*cos(x3)*dt; 
%     0 0 0];

R = ([1 0 0 0 0 0 
    0 0.1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 0.1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 0.1].^2)./1000;

[RE, Re] = eig (R);

% Measurement model defined below
Q = [1 0 0
    0 1 0
    0 0 3];

% Simulation Initializations
Tf = 10;
T = 0:dt:Tf;
% n = length(Ad(1,:));
x = zeros(6,length(T));
x(:,1) = x0;
% m = length(Q(:,1));
% y = zeros(m,length(T));
% mup_S = zeros(n,length(T));
% mu_S = zeros(n,length(T));
w1 = 10;
w2 = 15;
r = 0.25;
l = 0.3;
figure;
hold on;

%% Main loop
for t=2:length(T)
    %% Simulation
    % Select a motion disturbance
    e = RE*sqrt(Re)*randn(6,1);
    % Update state
%     x(:,t) = x(:,t-1) + Ad*x(:,t-1) + e;
    x(:,t) = [
        x(1,t-1) + dt*x(2,t-1)
        (w1+w2)*r/2*cos(x(5,t-1))
        x(3,t-1) + dt*x(4,t-1)
        (w1+w2)*r/2*sin(x(5,t-1))
        x(5,t-1) + dt*x(6,t-1)
        (r*(w1-w2)/(2*l))
        ]+e;
    
    x1 = x(1,t);
    x2 = x(3,t);
    plot(x(1,t), x(3,t), 'b*');
    
% 
%     % Take measurement
%     % Select a motion disturbance
     d = sqrt(Q)*randn(3,1);
%     % Determine measurement
%     
%     
     y(:,t) = [x(1,t); x(3,t);  99.7 - x(5,t)] + d;
% 
% 
%     %% Extended Kalman Filter Estimation
%     % Prediction update
     mup = [
        x(1,t-1) + dt*x(2,t-1)
        (w1+w2)*r/2*cos(x(5,t-1))
        x(3,t-1) + dt*x(4,t-1)
        (w1+w2)*r/2*sin(x(5,t-1))
        x(5,t-1) + dt*x(6,t-1)
        (r*(w1-w2)/(2*l))
        ];
    
    Ad = [
        1 dt 0 0 0 0
        0 0 0 0 -(w1+w2)/2*r*sin(x(5,t-1)) 0
        0 0 1 dt 0 0
        0 0 0 0 (w1+w2)/2*r*cos(x(5,t-1)) 0
        0 0 0 0 1 dt
        0 0 0 0 0 0];

     Sp = Ad*S*Ad' + R;  
          
% 
%     % Linearization
     Ht = [1 0 0 0 0 0
         0 0 1 0 0 0 
         0 0 0 0 -1 0];
     
% 
%     % Measurement update
     K = Sp*Ht'*inv(Ht*Sp*Ht'+Q);
     h = [mup(1)
         mup(3)
         99.7 - mup(5)];
     
     mu = mup + K*(y(:,t)-h);
     S = (eye(6)-K*Ht)*Sp;
     plot(mu(1), mu(3), 'r*');
% 
%     % Store results
%     mup_S(:,t) = mup;
%     mu_S(:,t) = mu;
%     K_S(:,t) = K;
% 
% 
%     %% Plot results
%     figure(1);clf; hold on;
%     plot(0,0,'bx', 'MarkerSize', 6, 'LineWidth', 2)
%     plot([20 -1],[0 0],'b--')
%     plot(x(1,2:t),x(3,2:t), 'ro--')
%     plot(mu_S(1,2:t),mu_S(3,2:t), 'bx--')
%     mu_pos = [mu(1) mu(3)];
%     S_pos = [S(1,1) S(1,3); S(3,1) S(3,3)];
%     error_ellipse(S_pos,mu_pos,0.75);
%     error_ellipse(S_pos,mu_pos,0.95);
%     title('True state and belief')
%     axis([-1 20 -1 10])
%     F(t-1) = getframe;
end

