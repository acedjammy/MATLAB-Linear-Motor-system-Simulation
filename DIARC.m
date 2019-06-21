clear
close all
clc
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
%% Parameter setup
global u U M B Asc Acog1 Acog3 omega_y dt;
% u=0; M = 0.025; B = 0.1; Asc = 0.1; Acog1 = 0.01; Acog3 = 0.05; d0 = 1;
u=0; M = 0.085; B = 0.35; Asc = 0.15; Acog1 = 0.05; Acog3 = 0.05; d0 = 1;
Mc = 0.055; Bc = 0.225; Ascc = 0.125; Acog1c = 0.03; Acog3c = 0.03; d = 0;
omega_y= 2*pi/0.06;
%% Trajectory Setup
ym=0; ym_dot=0; uc=0.2; omega_m=15; zeta_m=1.0; %reference model
%% Controller Setup
k1=600;k2=10;
phi = zeros(6,1);
theta = zeros(6,1);
Gamma=diag([10.0, 500.0, 100, 100, 500, 0.2]);
theta(1) = Bc; theta(2) = Ascc; theta(3) = Acog1c; 
theta(4) = Acog3c; theta(5) = d; theta(6) = Mc; 
S=saturation([-1 1]);kf=1000; P = 0.06;
theta_min = [0.1;0.1;0.01;0.01;-2;0.025];
theta_max = [0.35;0.15;0.05;0.05;2;0.085];
epsilon = 1;
yf_dot=0; mu1=10; mu2=1; tau_f=0.01; Np=6; dcm=2;gammad=32036;dc=0;
%% Euler Approximation update setup
X = zeros(8,1);
X_dot = zeros(8,1);
%% Discret system
for time_index = 1:T/time_step    
    time(time_index) = time_step * (time_index - 1);
    time_index;
    Time(end+1) = time_index * time_step; 
    tt = time_index * time_step
%% Generate model reference
    y = x1;
    y_dot = x2;
    ym = X(1);
    ym_dot = X(2);
    ym_ddot = -omega_m^2*ym-2*zeta_m*omega_m*ym_dot+omega_m^2*uc;

%% Controller
    phi = [y_dot;evaluate(S,kf*y_dot);-sin(omega_y*y);-sin(3*omega_y*y);-1;(ym_ddot-k1*(y_dot-ym_dot))];
    p=(y_dot-ym_dot)+k1*(y-ym);
    h = abs(phi)'*(theta_max-theta_min)+2;   
    u = phi'*theta -k2*p - h^2/(4*epsilon)*p-dc;
    if (abs(u) > 10)
        u = sign(u) * 10;
    end
    U(end+1) = u; 
%     dt = d0 +(-1)^(round(10*tt*sin(20*tt)));   
    dt = d0;
%% Update dc phi_uf, theta
    if (abs(dc) >= dcm && dc*p>0)
        dc_dot=0;
    else
        dc_dot = gammad*p;
    end
    
    uf = X(3);
    uf_dot = (-uf+u)/tau_f;

    yf_dot = X(4);
    yf_ddot = (-yf_dot+y_dot)/tau_f;

    Ascf = X(5);
    Ascf_dot = (-Ascf+evaluate(S,kf*y_dot))/tau_f;

    Acog1f = X(6);
    Acog1f_dot = (-Acog1f-sin(2*pi/P*y))/tau_f;

    Acog3f = X(7);
    Acog3f_dot = (-Acog3f-sin(6*pi/P*y))/tau_f;

    df = X(8);
    df_dot = (-df-1)/tau_f;

    phi_uf = [yf_dot;Ascf;Acog1f;Acog3f;df;yf_ddot];

    err_u=uf - phi_uf'*theta;
    theta_dot=Gamma*phi_uf*err_u;
    
    theta_p = theta_dot * time_step + theta;
        
    for j=1:6
        if (theta_p(j) < theta_max(j)) && (theta_p(j) > theta_min(j))
            theta(j)=theta_p(j);
        end
    end
    Gamma=(Gamma-(mu2*time_step/(1-mu1*time_step))*Gamma*phi_uf/ ...
    (1+(mu2*time_step/(1-mu1*time_step))*phi_uf'*Gamma*phi_uf)*phi_uf'*Gamma)/(1-mu1*time_step);

    Xdot = [ym_dot;ym_ddot;uf_dot;yf_ddot;Ascf_dot;Acog1f_dot;Acog3f_dot;df_dot];
    X = X + Xdot * time_step;
    dc = dc+dc_dot*time_step;
%% Save history
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
    
%% Simulate discrete
    tspan = [time(time_index) time(time_index)+time_step];
    [t, y] = ode45(@sys, tspan, xp);
    [NN,MM]=size(y);
    xp=y(NN,:); x1=xp(1,1); x2=xp(1,2); 
    
    X1(end+1)=x1; X2(end+1)=x2;
    
    if (mod(time_index*time_step,Pt) == 0)
       NNu=round(time_index*time_step/Pt); uc=uc+(-1)^NNu*0.2;       
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
    global u M B Asc Acog1 Acog3 omega_y dt;
    S=saturation([-1 1]);
    dx(1) = x(2);
    dx(2) = (u - B * x(2) - Asc * evaluate(S,1000*x(2)) ...
        + Acog1 * sin(omega_y*x(1)) + Acog3 ...
        * sin(3*omega_y*x(1)) + dt ) / M;
    dx = dx';
    
end