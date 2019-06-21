clear all
close all
clc

%% ODE 45
global u U index ek2 ek1 ek last_t M B Asc Acog1 Acog3 omega_y Case;
u=0;index=1;ek2 = 0;ek1 = 0;ek = 0;last_t = 0;
M = 0.085; B = 0.35; Asc = 0; Acog1 = 0; Acog3 = 0;
omega_y= 2*pi/0.06;Case = 1;

omega_m = 15; zeta_m = 1;
yf = 0; yfdot = 0; yfddot = 0; uf = 0; ufdot = 0; ufddot = 0; 
b0 = 1/0.225; a1 = 0.225/0.055; bm0 = omega_m^2; a0 = 100; am1 = 2*zeta_m*omega_m; am2 = omega_m^2;

Pm = 100;
P = eye(4)*Pm;
Pdot = zeros(4,4);

theta = [0.1;0.1;0.1;0.1]; thetadot = [0;0;0;0];
uk = 0;
alpha = 0.995;

X = zeros(21,1);
% X(7) = 43.8536719069886;
% X(8) = 5026.65770886679;
% X(9) = 2722.37886216010;
% X(10) = 22498.2602437943;
X(7) = 18.1818;
X(8) = 3500;
X(9) = 2709.92;
X(10) = 22500;


X(11) = Pm;
X(15) = Pm;
X(18) = Pm;
X(20) = Pm;

time_step=0.0005; %sampling period
T=4; % simualtion time span
xp = [0;0];
x1 = xp(1);
x2 = xp(2);
X1(1) = xp(1);
X2(1) = xp(2);
U = 0;
Time = 0;
ss = 0;
tolerance = 0.15;
ym=0; ym_dot=0; uc=0.2; omega_m=15; zeta_m=1.0; %reference model
xc=[ym;ym_dot];
% save_ym = 0;

%% Discret system
for time_index = 1:T/time_step    

    time(time_index) = time_step * (time_index - 1);
    time_index;
    Time(end+1) = time_index * time_step; 
       
%% Build Trajector for step tracking
    ym=xc(1); ym_dot=xc(2);
    ym_ddot=-omega_m^2*ym-2*zeta_m*omega_m*ym_dot+omega_m^2*uc;
    xc_dot(1,1)=ym_dot;
    xc_dot(2,1)=ym_ddot;
    xc=xc+xc_dot*time_step;
%     save_ym(end+1)=xc(1); %save the reference output 
    
    
%% Update states
y = x1;
yf = X(1);
yfdot = X(2);
yfddot = X(3);
yfdddot = -(am1+a0) * yfddot - (am2+am1*a0) * yfdot - am2*a0*yf + am2*a0*y;

uf = X(4);
ufdot = X(5);
ufddot = X(6);
ufdddot = -(am1+a0) * ufddot - (am2+am1*a0) * ufdot - am2*a0*uf + am2*a0*u;
theta(1) = X(7);
theta(2) = X(8);
theta(3) = X(9);
theta(4) = X(10);
P(1,1) = X(11);
P(1,2) = X(12);
P(2,1) = X(12);
P(1,3) = X(13);
P(3,1) = X(13);
P(1,4) = X(14);
P(4,1) = X(14);
P(2,2) = X(15);
P(2,3) = X(16);
P(3,2) = X(16);
P(2,4) = X(17);
P(4,2) = X(17);
P(3,3) = X(18);
P(3,4) = X(19);
P(4,3) = X(19);
P(4,4) = X(20);
uk = X(21);

phi = [ufdot;uf;yfdot;yf];
e = am2*a0*y-phi'*theta
thetadot = P*phi*e;
A = P*(phi*phi')*P;
Pdot = alpha*P - A;

ukdot = -theta(2)/theta(1)*uk ...
        + (bm0*a0/theta(1) - bm0*theta(2)/theta(1)^2) * uc ...
        - (theta(4)/theta(1)-theta(3)*theta(2)/theta(1)^2) * y;


Xdot = [yfdot; yfddot; yfdddot; ufdot; ufddot; ufdddot; thetadot(1); thetadot(2); thetadot(3); ...
        thetadot(4); Pdot(1,1); Pdot(1,2); Pdot(1,3); Pdot(1,4); Pdot(2,2); Pdot(2,3); Pdot(2,4); ...
        Pdot(3,3); Pdot(3,4); Pdot(4,4); ukdot];

X = X + Xdot * time_step;    





%% Controller

    u = bm0/theta(1) * uc - theta(3)/theta(1) * y + uk;
    if (abs(u) > 10)
        u = sign(u) * 10;
    end
    U(end+1) = u;
    
    

    
    tspan = [time(time_index) time(time_index)+time_step];
    [t, y] = ode45(@sys, tspan, xp);
    [NN,MM]=size(y);
    xp=y(NN,:); x1=xp(1,1); x2=xp(1,2); 
    
    X1(end+1)=x1; X2(end+1)=x2;
    
    if (mod(time_index*time_step,0.6) == 0)
       NN=floor(time_index*time_step/0.6); uc=uc+(-1)^NN*0.2;
    end
    theta0 = [11.7647;1480.97;2706.66;22500];
    save_u(time_index)=u; %save control input for plotting
    save_y(time_index)=x1; %save output
    save_ym(time_index)=ym; %save the reference output
    save_uc(time_index)=uc; %save the reference command input
    save_theta(:,time_index)=theta; %save parameter estimates
    save_theta0(:,time_index)=theta0; %save parameter estimates
end


% r = 0.1*(1-cos(4*pi*Time));
% r = save_ym;
% 
% subplot(2,2,1)
% plot(Time,X1,Time,r)
% legend('x1','r')
% xlabel('time')
% ylabel('Position (m)')
% title('Position')
% subplot(2,2,2)
% plot(Time,U)
% legend('u')
% xlabel('time')
% ylabel('Voltage')
% title('Voltage Input')
% subplot(2,2,3)
% plot(Time,X2)
% legend('x2')
% xlabel('time')
% ylabel('Speed')
% title('Speed')
% subplot(2,2,4)
% plot(Time,r-X1)
% legend('error')
% xlabel('time')
% ylabel('Error')
% title('Error')
% 
% figure 
% plot(Time,ss)
% title('s')
% xlabel('time (sec)')
TT = time;

figure
hold on
plot(TT,save_uc)
plot(TT,save_ym)
plot(TT,save_y)
hold off
title("Output, Reference and Uc")
xlabel("time (sec)")
legend('Reference Input','Reference','Actual Position')

figure
plot(save_y-save_ym)
title("Tracking error")
xlabel("time (sec)")
ylabel("error in position (m)")

figure
plot(save_u)
title("Control Input")
xlabel("time (sec)")
ylabel("Voltage (volt)")

figure
subplot(2,2,1)
plot(TT,save_theta(1,:))
hold on
plot(TT,save_theta0(1,:))
title("Convergence of r0t")
legend('Estimated r0t','true r0t')
subplot(2,2,2)
plot(TT,save_theta(2,:))
hold on
plot(TT,save_theta0(2,:))
title("Convergence of r1t")
legend('Estimated r1t','true r1t')
subplot(2,2,3)
plot(TT,save_theta(3,:))
hold on
plot(TT,save_theta0(3,:))
title("Convergence of s0t")
legend('Estimated s0t','true s0t')
subplot(2,2,4)
plot(TT,save_theta(4,:))
hold on
plot(TT,save_theta0(4,:))
title("Convergence of s1t")
legend('Estimated s1t','true s1t')


function  dx = sys(t, x)

    global u M B Asc Acog1 Acog3 omega_y;
    S=saturation([-1 1]);
    dx(1) = x(2);
    dx(2) = (u - B * x(2) - Asc * evaluate(S,1000*x(2)) ...
        - Acog1 * sin(omega_y*x(1)) - Acog3 ...
        * sin(3*omega_y*x(1))) / M;
    dx = dx';
    
end





