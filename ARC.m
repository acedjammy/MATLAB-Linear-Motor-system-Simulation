clear
close all
clc

%% Parameter setup
global u U index ek2 ek1 ek last_t M B Asc Acog1 Acog3 omega_y Case d0;
u=0;index=1;ek2 = 0;ek1 = 0;ek = 0;last_t = 0;
M = 0.025; B = 0.1; Asc = 0.1; Acog1 = 0.01; Acog3 = 0.05; d0 = 1;
Mc = 0.055; Bc = 0.225; Ascc = 0.125; Acog1c = 0.03; Acog3c = 0.03; d = 0;
omega_y= 2*pi/0.06;
%% Trajectory Setup
ym=0; ym_dot=0; uc=0.2; omega_m=15; zeta_m=1.0; %reference model
%% Controller Setup
k1=100;k2=20;
phi = zeros(6,1);
theta = zeros(6,1);
Gamma=diag([3.1003 1.035 0.0551 0.0551 244.9648 0.1853]);
theta(1) = Bc; theta(2) = Ascc; theta(3) = Acog1c; 
theta(4) = Acog3c; theta(5) = d; theta(6) = Mc; 
S=saturation([-1 1]);kf=1000; P = 0.06;
theta_min = [0.1;0.1;0.01;0.01;-2;.025];
theta_max = [0.35;0.15;0.05;0.05;2;0.085];
epsilon = 1;
%% Euler Approximation update setup
X = [zeros(2,1);theta];
X_dot = zeros(8,1);
%% Loop Setup
time_step=0.0005; %sampling period
T=4; % simualtion time span
Pt = 0.6;
xp = [0;0];
x1 = xp(1);
x2 = xp(2);
X1(1) = xp(1);
X2(1) = xp(2);
U = 0;
Time = 0;


%% Discret system
for time_index = 1:T/time_step    

    time(time_index) = time_step * (time_index - 1);
    time_index;
    Time(end+1) = time_index * time_step; 
    tt = time_index * time_step
%% Update states step 1
y = x1;
y_dot = x2;
ym = X(1);
ym_dot = X(2);
ym_ddot = -omega_m^2*ym-2*zeta_m*omega_m*ym_dot+omega_m^2*uc;

phi = [-y_dot;-evaluate(S,kf*y_dot);sin(2*pi/P);sin(6*pi/P);-1;-(ym_ddot-k1*(y_dot-ym_dot))];
p=(y_dot-ym_dot)+k1*(y-ym);



% theta = X(3:8,1);
theta_dot = Gamma * phi * p;

Xdot = [ym_dot;ym_ddot;theta_dot];

%% Controller
    h = abs(phi)'*(theta_max-theta_min)+2;
    u = -phi'*theta-k2*p - h^2/(4*epsilon)*p;
    if (abs(u) > 10)
        u = sign(u) * 10;
    end
    U(end+1) = u; 
    d0 = 1+(-1)^(round(10*tt*sin(20*tt)));
    
    tspan = [time(time_index) time(time_index)+time_step];
    [t, y] = ode45(@sys, tspan, xp);
    [NN,MM]=size(y);
    xp=y(NN,:); x1=xp(1,1); x2=xp(1,2); 
    
    X1(end+1)=x1; X2(end+1)=x2;
    
    if (mod(time_index*time_step,Pt) == 0)
       NNu=round(time_index*time_step/Pt); uc=uc+(-1)^NNu*0.2;       
    end
    
    
    
    save_u(time_index)=u; %save control input for plotting
    save_y(time_index)=x1; %save output
    save_ym(time_index)=ym; %save the reference output
    save_uc(time_index)=uc; %save the reference command input
    
    save_Mc(:,time_index)=theta(6); 
    save_Bc(:,time_index)=theta(1); 
    save_Ascc(:,time_index)=theta(2);
    save_Acog1c(:,time_index)=theta(3);
    save_Acog3c(:,time_index)=theta(4);
    save_dc(:,time_index)=theta(5); 
    
    save_M(:,time_index)=M; 
    save_B(:,time_index)=B; 
    save_Asc(:,time_index)=Asc;
    save_Acog1(:,time_index)=Acog1;
    save_Acog3(:,time_index)=Acog3;
    save_d(:,time_index)=d0;
 %% Update states step 2   
    X = X + Xdot * time_step;
    theta_p = X(3:8,1);
    for j=1:6
    if (theta_p(j) < theta_max(j)) & (theta_p(j) > theta_min(j))
    theta(j)=theta_p(j);
    end
    end
end

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
plot(TT,save_y-save_ym)
title("Tracking error")
xlabel("time (sec)")
ylabel("error in position (m)")

figure
plot(TT,save_u)
title("Control Input")
xlabel("time (sec)")
ylabel("Voltage (volt)")

figure
subplot(2,3,1)
plot(TT,save_Mc)
hold on
plot(TT,save_M)
title("Convergence of M")
legend('Estimated M','true M')
subplot(2,3,2)
plot(TT,save_Bc)
hold on
plot(TT,save_B)
title("Convergence of B")
legend('Estimated B','true B')
subplot(2,3,3)
plot(TT,save_Ascc)
hold on
plot(TT,save_Asc)
title("Convergence of Asc")
legend('Estimated Asc','true Asc')
subplot(2,3,4)
plot(TT,save_Acog1c)
hold on
plot(TT,save_Acog1)
title("Convergence of Acog1")
legend('Estimated Acog1','true Acog1')
subplot(2,3,5)
plot(TT,save_Acog3c)
hold on
plot(TT,save_Acog3)
title("Convergence of Acog3")
legend('Estimated Acog3','true Acog3')
subplot(2,3,6)
plot(TT,save_dc)
hold on
plot(TT,save_d)
title("Convergence of d0")
legend('Estimated d0','true d0')



function  dx = sys(t, x)

    global u M B Asc Acog1 Acog3 omega_y d0;
    S=saturation([-1 1]);
    dx(1) = x(2);
    dx(2) = (u - B * x(2) - Asc * evaluate(S,1000*x(2)) ...
        - Acog1 * sin(omega_y*x(1)) - Acog3 ...
        * sin(3*omega_y*x(1)) + d0 ) / M;
    dx = dx';
    
end