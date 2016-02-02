clear all
close all
clc

Ex = 1;
Exx = 0.1;
Ey = 1;
Eyy = 0.1;
Eth = 1*pi/180;
Ew = 0.1*pi/180;
EE = 1;
EN = 1;
EM = 3*pi/180;
dt = 0.1;
m = 1;
r = 0.25;
l = 0.3;
I = 0.2;
T1 = 0;
T2 = 0;
II = eye(6);
% w = 10;
% w1 = 15;
% w2 = 2*l/r*w + w1;
% v1 = w1*r;
% v2 = w2*r;
% v = (v1 + v2)/2;
vx = 3;
vy = 0;
w = 0;
theta = atan2(vy,vx);
Mu = [0; vx; 0; vy; theta; w];
x(:,1) = Mu;
y(:,1) = [x(1,1); x(3,1); 99.7*pi/180 - x(5,1)];
E = zeros(6)+0.1;
T = 60;
R = 0*dt*[Ex 0 0 0 0 0
    0 Exx 0 0 0 0
    0 0 Ey 0 0 0
    0 0 0 Eyy 0 0
    0 0 0 0 Eth 0
    0 0 0 0 0 Ew];
[RE,Re] = eig(R);
Q = 0*dt*[EE 0 0
    0 EN 0
    0 0 EM];
[REy,Rey] = eig(Q);
check = 1;
t = 1;
% Fixed vehicle parameters
v = 3; % Speed
delta_max = 360*pi/180; % max steering angle
k = 2.5; % Gain
% Simulation time
Tmax = 3;  % End point
% dt =0.001; % Time step
T = 0:dt:Tmax; % Time vector

% Desired line state through 0,0 [theta]
xp = [0]; % Line heading

% Initial conditions in [e psi]
x0 = [1 0]; % translational offset
%x0 = [0 2]; % heading offset
%x0 = [5 2]; % both
% Simulation setup
xd = zeros(length(T)-1,2); % Derivative of state ([edot psidot])
xc = zeros(length(T),2);  % State ([e psi]
xc(1,:) = x0; % Initial condition
delta = zeros(length(T),1); % Steering angles
p = zeros(length(T),2); % Position in x,y plane
p(1,:) = (x0(1))*[sin(xp) cos(xp)]; % Initial position
x(1,1) = p(1,1);
x(3,1) = p(1,2);

% % % delta(t) = k*xc(t,1)/v;
% % % % State derivatives
% % % xd(t,1) = v*sin(xc(t,2));
% % % xd(t,2) = -delta(t);%(v1-v2)/(l);
% % % w = xd(t,2);
% Calculate steering angle
delta(t) = max(-delta_max,min(delta_max,xc(t,2) + atan2(k*xc(t,1),v)));
% State derivatives
xd(t,1) = v*sin(xc(t,2) - delta(t));
xd(t,2) = (delta(t) - 2*xc(t,2))/dt;
w = xd(t,2);
% State update
xc(t+1,1) = xc(t,1)+dt*xd(t,1);
xc(t+1,2) = xc(t,2)+dt*xd(t,2);
% Position update
p(t+1,1) = p(t,1) + dt*v*cos(xc(t,2)-delta(t)-xp);
p(t+1,2) = p(t,2) + dt*v*sin(xc(t,2)-delta(t)-xp);
x(1,t+1) = p(t+1,1);
x(3,t+1) = p(t+1,2);

% while(check ~= 5)
while (t <= (Tmax/dt))
    t = t + 1;
    ex = RE*(Re)*randn(6,1);
    ey = REy*(Rey)*randn(3,1);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %% Trajectory tracking
    
    % Desired line state through 0,0 [theta]
    xp = [0]; % Line heading
    
    % Initial conditions in [e psi]
    x0 = [1 0]; % translational offset
    %x0 = [0 2]; % heading offset
    %x0 = [5 2]; % both
    
    
    
    %     for t=1:length(T)-1
    % % %     % Calculate steering angle
    % % %     delta(t) = atan2(2*k*xc(t,1),v);%xd(t-1,1)/(1+(k*xc(t,1))^2);%k*xc(t,1)/v;
    % % %     % State derivatives
    % % %     xd(t,1) = v*sin(xc(t,2));
    % % %     xd(t,2) = -delta(t);%(v1-v2)/(l);
    % % %     w = xd(t,2);
    % Calculate steering angle
    delta(t) = max(-delta_max,min(delta_max,xc(t,2) + atan2(k*xc(t,1),v)));
    % State derivatives
    xd(t,1) = v*sin(xc(t,2) - delta(t));
    xd(t,2) = (delta(t) - 2*xc(t,2))/dt;
    w = xd(t,2);
    % State update
    xc(t+1,1) = xc(t,1)+dt*xd(t,1);
    xc(t+1,2) = xc(t,2)+dt*xd(t,2);
%     [xc(t+1,1), xc(t+1,2)]
    % Position update
        p(t+1,1) = p(t,1) + dt*v*cos(xc(t,2)-delta(t)-xp);
        p(t+1,2) = p(t,2) + dt*v*sin(xc(t,2)-delta(t)-xp);
%         end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %     if ((x(1,t-1)-x(1,1) >= 20) && check ==1)
    %         vx = 0;
    %         vy = 3;
    %         %         w = 0;
    %         %         theta = atan2(vy,vx);
    %         check = 2;
    %         t2 = t;
    %     elseif ((x(3,t-1)-x(3,t2) >= 5) && check ==2)
    %         vx = -3;
    %         vy = 0;
    %         %         w = 0;
    %         %         theta = atan2(vy,vx);
    %         check = 3;
    %         t2 = t;
    %     elseif ((x(1,t-1)-x(1,t2) <= -20) && check ==3)
    %         vx = 0;
    %         vy = -3;
    %         %         w = 0;
    %         %         theta = atan2(vy,vx);
    %         check = 4;
    %         t2 = t;
    %     elseif ((x(3,t-1)-x(3,t2) <= -5) && check ==4)
    %         vx = 0;
    %         vy = -3;
    %         %         w = 0;
    %         %         theta = atan2(vy,vx);
    %         check = 5;
    %     end
    
    %     x(:,t) = [x(1,t-1) + x(2,t-1)*dt
    %         x(2,t-1) + (-v1/w-l)*(x(6,t-1)^2)*sin(x(5,t-1))*dt + 1/(m*r)*(T1 + T2)*cos(x(5,t-1))*dt
    %         x(3,t-1) + (x(4,t-1))*dt
    %         x(4,t-1) - (-v1/w-l)*((x(6,t-1))^2)*cos(x(5,t-1))*dt + 1/(m*r)*(T1 + T2)*sin(x(5,t-1))*dt
    %         x(5,t-1) + x(6,t-1)*dt
    %         x(6,t-1) + l/(r*I)*(-T1 + T2)*dt] + ex;
%     x(:,t) = [x(1,t-1) + x(2,t-1)*dt
%         v*cos(x(5,t-1))
%         x(3,t-1) + (x(4,t-1))*dt
%         v*sin(x(5,t-1))
%         x(5,t-1) + x(6,t-1)*dt
%         w] + ex;
    x(:,t+1) = [p(t+1,1)
        v*cos(x(5,t))
        p(t+1,2)
        v*sin(x(5,t))
        x(5,t) + x(6,t)*dt
        w] + ex;

    
    y(:,t) = [x(1,t); x(3,t); 99.7*pi/180 - x(5,t)] + ey;
    
    figure(1)
    plot(x(1,t),x(3,t),'*r')
    hold on
%     figure(2)
%     plot(p(t,1),p(t,2),'r');
%     hold on
%     %     plot(y(1,t),y(2,t),'*g')
%     %     plot (Mu(1),Mu(3),'*')
%     %     legend('Motion Model', 'Measurement Model', 'Filtered Image')
%     drawnow
end
figure(1)
xlabel('x [m]')
ylabel('y [m]')
title('Simulating the Motion of a Rebot and Using the EKF to Correct with the Help of the Measurement')
hold on
for t =2:Tmax/dt
    % % % % % % % % % % %
    % Prediction Update
    % % % % % % % % % % %
    G = [1, dt, 0, 0, 0, 0
        0, 1, 0, 0, (-v1/w - l)*(x(6,t-1)^2)*cos(x(5,t-1))*dt - 1/(m*r)*(T1 + T2)*sin(x(5,t-1))*dt, 2*(-v1/w-l)*(x(6,t-1))*sin(x(5,t-1))*dt
        0, 0, 1, dt, 0, 0
        0, 0, 0, 1, (-v1/w - l)*(x(6,t-1)^2)*sin(x(5,t-1))*dt + 1/(m*r)*(T1 + T2)*cos(x(5,t-1))*dt, -2*(-v1/w-l)*(x(6,t-1))*cos(x(5,t-1))*dt
        0, 0, 0, 0, 1, dt
        0, 0, 0, 0, 0, 1];
    
    Mup = [Mu(1) + (Mu(2))*dt
        Mu(2) + (-v1/w-l)*((Mu(6))^2)*sin(Mu(5))*dt + 1/(m*r)*(T1 + T2)*cos(Mu(5))*dt
        Mu(3) + (Mu(4))*dt
        Mu(4) - (-v1/w-l)*((Mu(6))^2)*cos(Mu(5))*dt + 1/(m*r)*(T1 + T2)*sin(Mu(5))*dt
        Mu(5) + (Mu(6))*dt
        Mu(6) + l/(r*I)*(-T1 + T2)*dt];
    
    
    R = [Ex 0 0 0 0 0
        0 Exx 0 0 0 0
        0 0 Ey 0 0 0
        0 0 0 Eyy 0 0
        0 0 0 0 Eth 0
        0 0 0 0 0 Ew];
    
    Ep = G*E*G' + R;
    
    % % % % % % % % % % %
    % Measurement Update
    % % % % % % % % % % %
    H = [1 0 0 0 0 0
        0 0 1 0 0 0
        0 0 0 0 -1 0];
    
    Q = [EE 0 0
        0 EN 0
        0 0 EM];
    
    K = Ep*H'*(H*Ep*H' + Q)^-1;
    h = [Mup(1); Mup(3); 99.7*pi/180 - Mup(5)];
    Mu = Mup + K*(y(:,t) - h);
    E = (II - K*H)*Ep;
    
    
    figure(1)
    plot(x(1,t),x(3,t),'*r')
    plot(y(1,t),y(2,t),'*g')
    plot (Mu(1),Mu(3),'*')
    legend('Motion Model', 'Measurement Model', 'Filtered Image')
    drawnow
end