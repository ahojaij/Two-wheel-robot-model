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
dt = 0.02;
m = 1;
r = 0.25;
l = 0.3;
I = 0.2;
T1 = 0;
T2 = 0;
II = eye(6);
w = 10;
w1 = 15;
w2 = 2*l/r*w + w1;
v1 = w1*r;
v2 = w2*r;
v = (v1 + v2)/2;
theta = 99.7*pi/180;
Mu = [0; v*cos(theta); 0; v*sin(theta); theta; w];
x(:,1) = Mu;
y(:,1) = [x(1,1); x(3,1); 99.7*pi/180 - x(5,1)];
E = zeros(6)+0.1;
T = 60;
t1 = 0;
R = dt/100*[Ex 0 0 0 0 0
    0 Exx 0 0 0 0
    0 0 Ey 0 0 0
    0 0 0 Eyy 0 0
    0 0 0 0 Eth 0
    0 0 0 0 0 Ew];
[RE,Re] = eig(R);
Q = dt/100*[EE 0 0
    0 EN 0
    0 0 EM];
[REy,Rey] = eig(Q);
for t = 2:T/dt
    ex = RE*sqrt(Re)*randn(6,1);
    ey = REy*sqrt(Rey)*randn(3,1);
    
    x(:,t) = [x(1,t-1) + x(2,t-1)*dt
        x(2,t-1) + (-v1/x(6,t-1)-l)*(x(6,t-1)^2)*sin(x(5,t-1))*dt + 1/(m*r)*(T1 + T2)*cos(x(5,t-1))*dt
        x(3,t-1) + (x(4,t-1))*dt
        x(4,t-1) - (-v1/x(6,t-1)-l)*((x(6,t-1))^2)*cos(x(5,t-1))*dt + 1/(m*r)*(T1 + T2)*sin(x(5,t-1))*dt
        x(5,t-1) + x(6,t-1)*dt
        x(6,t-1) + l/(r*I)*(-T1 + T2)*dt] + ex;
    

    if t*dt >= t1
        y(:,t) = [x(1,t); x(3,t); 99.7*pi/180 - x(5,t)] + ey;
        t1 = t1 + 0.06;
    else
        y(:,t) = y(:,t-1);
    end
end
figure(1)
xlabel('x [m]')
ylabel('y [m]')
title('Simulating the Motion of a Rebot and Using the EKF to Correct with the Help of the Measurement')
hold on
for t =2:T/dt
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
    
    Ep = G*E*G' + R;
    
    % % % % % % % % % % %
    % Measurement Update
    % % % % % % % % % % %
    H = [1 0 0 0 0 0
        0 0 1 0 0 0
        0 0 0 0 -1 0];
    
    K = Ep*H'*(H*Ep*H' + Q)^-1;
    h = [Mup(1); Mup(3); 99.7*pi/180 - Mup(5)];
    Mu = Mup + K*(y(:,t) - h);
    E = (II - K*H)*Ep;
    
    
    
    figure(1)
    plot(x(1,t),x(3,t),'*r')
    plot(y(1,t),y(2,t),'*g')
    plot (Mu(1),Mu(3),'*')
    legend('Motion Model', 'Measurement Model', 'Filtered Image')
        mu_pos = [Mu(1) Mu(3)];
        S_pos = [E(1,1) E(1,3); E(3,1) E(3,3)];
    %     error_ellipse(S_pos,mu_pos,0.75);
        error_ellipse(S_pos,mu_pos,0.95);
    drawnow
end