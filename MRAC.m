clear all
close all
clc

%% ODE 45
global u U index ek2 ek1 ek last_t M B Asc Acog1 Acog3 omega_y Case;
u=0;index=1;ek2 = 0;ek1 = 0;ek = 0;last_t = 0;
M = 0.085; B = 0.35; Asc = 0.15; Acog1 = 0; Acog3 = 0;
Mc = 0.055; Bc = 0.225; Ascc = 0.125;
omega_y= 2*pi/0.06;
Case = 1;

omega_m = 15; zeta_m = 1;
yf = 0; yfdot = 0; yfddot = 0; uf = 0; ufdot = 0; ufddot = 0; 

k1=100;k2=20;
gamma_M=0.4; gamma_B=40.0; gamma_F=1000.0; %adaptation rate

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
X=zeros(5,1);
X(3) = 0.055;
X(4) = 0.225;
X(5) = 0.125;
% save_ym = 0;

%% Discret system
for time_index = 1:T/time_step    

    time(time_index) = time_step * (time_index - 1);
    time_index;
    Time(end+1) = time_index * time_step; 
    t = time_index * time_step
%% Build Trajector for step tracking
%     ym=xc(1); ym_dot=xc(2);
%     ym_ddot=-omega_m^2*ym-2*zeta_m*omega_m*ym_dot+omega_m^2*uc;
%     xc_dot(1,1)=ym_dot;
%     xc_dot(2,1)=ym_ddot;
%     xc=xc+xc_dot*time_step;
%     save_ym(end+1)=xc(1); %save the reference output 
    
    
%% Update states step 1
y = x1;
y_dot = x2;
ym = X(1);
ym_dot = X(2);
Mc = X(3);
Bc = X(4);
Ascc = X(5);
ym_ddot = -omega_m^2*ym-2*zeta_m*omega_m*ym_dot+omega_m^2*uc;
p=(y_dot-ym_dot)+k1*(y-ym);
Mc_dot = -gamma_M*(ym_ddot-k1*(y_dot-ym_dot))*p;
Bc_dot = -gamma_B*y_dot*p;
Ascc_dot = -gamma_F*sign(y_dot)*p;
Xdot = [ym_dot;ym_ddot;Mc_dot;Bc_dot;Ascc_dot];

%% Controller

    u = Bc*x2+Ascc*sign(x2)+Mc*(ym_ddot-k1*(x2-ym_dot))-k2*p;
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
    save_u(time_index)=u; %save control input for plotting
    save_y(time_index)=x1; %save output
    save_ym(time_index)=ym; %save the reference output
    save_uc(time_index)=uc; %save the reference command input
    save_Mc(:,time_index)=Mc; 
    save_Bc(:,time_index)=Bc; 
    save_Ascc(:,time_index)=Ascc;
    save_M(:,time_index)=M; 
    save_B(:,time_index)=B; 
    save_Asc(:,time_index)=Asc;
 %% Update states step 2   
    X = X + Xdot * time_step;
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
subplot(3,1,1)
plot(TT,save_Mc)
hold on
plot(TT,save_M)
title("Convergence of M")
legend('Estimated M','true M')
subplot(3,1,2)
plot(TT,save_Bc)
hold on
plot(TT,save_B)
title("Convergence of B")
legend('Estimated B','true B')
subplot(3,1,3)
plot(TT,save_Ascc)
hold on
plot(TT,save_Asc)
title("Convergence of Asc")
legend('Estimated Asc','true Asc')



function  dx = sys(t, x)

    global u M B Asc Acog1 Acog3 omega_y;
    S=saturation([-1 1]);
    dx(1) = x(2);
    dx(2) = (u - B * x(2) - Asc * evaluate(S,1000*x(2)) ...
        - Acog1 * sin(omega_y*x(1)) - Acog3 ...
        * sin(3*omega_y*x(1))) / M;
    dx = dx';
    
end