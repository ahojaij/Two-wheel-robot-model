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
dt = 0.01;
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
R = 0.01*[Ex 0 0 0 0 0
    0 Exx 0 0 0 0
    0 0 Ey 0 0 0
    0 0 0 Eyy 0 0
    0 0 0 0 Eth 0
    0 0 0 0 0 Ew];
[RE,Re] = eig(R);
Q = 0.01*[EE 0 0
    0 EN 0
    0 0 EM];
[REy,Rey] = eig(Q);
for t = 2:T/dt
    ex = RE*(Re)*randn(6,1);
    ey = REy*(Rey)*randn(3,1);
    x(:,t) = [x(1,t-1) + (x(2,t-1))*dt
        v*cos(x(5,t-1))
        x(3,t-1) + (x(4,t-1))*dt
        v*sin(x(5,t-1))
        x(5,t-1)+ (x(6,t-1))*dt
        x(6,t-1) + l/(r*I)*(-T1 + T2)*dt] + ex;
    
    y(:,t) = [x(1,t); x(3,t); 99.7*pi/180 - x(5,t)] + ey;
end

for t =2:T/dt
    % % % % % % % % % % %
    % Prediction Update
    % % % % % % % % % % %
    G = [1 dt 0 0 0 0
        0 0 0 0 -v*sin(x(5,t-1)) 0
        0 0 1 dt 0 0
        0 0 0 0 v*cos(x(5,t-1)) 0
        0 0 0 0 1 dt
        0 0 0 0 0 1];
    
    Mup = [Mu(1) + Mu(2)*dt
        v*cos(Mu(5))
        Mu(3) + Mu(4)*dt
        v*sin(Mu(5))
        Mu(5)+ Mu(6)*dt
        Mu(6) + l/(r*I)*(-T1 + T2)*dt] + 0.01*RE*(Re)*randn(6,1);
    
    %     Mup = [x(1,t-1) + x(2,t-1)*dt
    %         v*cos(x(5,t-1))
    %         x(3,t-1) + x(4,t-1)*dt
    %         v*sin(x(5,t-1))
    %         x(5,t-1) + x(6,t-1)*dt
    %         x(6,t-1) + l/(r*I)*(-T1 + T2)*dt] + 0.01*RE*(Re)*randn(6,1);
    
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
    plot (Mu(1),Mu(3),'*')
    plot(y(1,t),y(2,t),'*g')
    plot(x(1,t),x(3,t),'*r')
    drawnow
    hold on;
end