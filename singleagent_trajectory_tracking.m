close all
clear
clc
%%
%Global Variables
global r a d mc mw mt Im Ic Iw I M_bar C_bar E_bar 

%Robot Global Variables
r = 0.1;    %radius of wheels in m
a = 0.3;    %distance between center of chasi of wheels to wheels
d = 0;%0.2; %distance between of center of mass and center of chasi
mc = 3;     %mass of the robot without wheels and motors in[kg]
mw = 1;     %mass of wheel in[kg]
Im = mc*a^2;%moment of inertia of robot
Ic = 11.4180;
Iw = 0.0456;
mt = mc + 2*mw;
I = mc*d^2 + Ic + 2*mw*(d^2+a^2) + 2*Im;
%%=================================%%
%Euler-lagrange matrix
M_bar = [Iw + r^2*(mt/4 - (mt*d^2)/(4*a^2) + I/(4*a^2)), r^2*(mt/4 + (mt*d^2)/(4*a^2) - I/(4*a^2));
     r^2*(mt/4 + (mt*d^2)/(4*a^2) - I/(4*a^2)), Iw + r^2*(mt/4 - (mt*d^2)/(4*a^2) + I/(4*a^2))];
inv_M_bar = inv(M_bar);
C_bar = zeros(2);
E_bar = eye(2);
inv_E_bar = inv(E_bar);
%%=================================%%
%Functions
%transform [v w]' to q_dot
Sv = @(theta) [cos(theta) -d*sin(theta); sin(theta) d*cos(theta); 0 1]; 
 %from inertia frame to local frame
rotaion_matrix = @(theta)[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
%map [V omega] to [phi_dotr phi_dotL]
Omega = [1/r a/r; 1/r -a/r]; 
%Constraint Matrix
S = @(theta)[(r/(2*a))*(a*cos(theta) - d*sin(theta)), (r/(2*a))*(a*cos(theta) + d*sin(theta));
    (r/(2*a))*(a*sin(theta) + d*cos(theta)), (r/(2*a))*(a*sin(theta)-d*cos(theta));
    r/(2*a), -r/(2*a);];
%%=================================%%
%% Init 
Phi = [0; 0];   %[phidot_r, phidot_L]
Phi_dot = [0; 0];   %[phi2dot_r, phi2dot_L]
q = [0; 0; 0];
q_dot = [0; 0; 0];
qe = [0 ;0; 0];
V = [0; 0];

%Vectors
q_vec = [];
q_dot_vec = [];
qr_vec = [];
qe_vec = [];
V_vec = [];
Phi_vec = [];
Vr_vec = [];
tau_vec = [];
%%=============================%%
N = 1; %number of agent

%Simulation Time
simul_time = 20;
tsteps = 0.001;
time_sim = 0:tsteps:simul_time;
%%=============================%%
%% Params
%Dynamic Controller param
Kp_dyn = 20;

%Reference
Vref = [1; 0.5];  %[vr, wr].' in (m/s)
vr = Vref(1, :);
wr = Vref(2, :);
r_path = vr/wr;

%Kinematic Control
K = [5; 10; 2]; %control parameters
%%=============================%%
%% Simulation
iteration = 1;
t = 0;
tvec = [];

while t <= simul_time
    
    %% Kinematic Desired
    %Posture Reference
    qr = [r_path*cos(wr*t); r_path*sin(wr*t); wr*t+pi/2];
    
    %Posture Error
    ex = qr(1, :) - q(1, :);
    ey = qr(2, :) - q(2, :);
    etheta = qr(3, :) - q(3, :);
    qe_temp = [ex; ey; etheta];
    qe = rotaion_matrix(q(3,:))*qe_temp;

    %Output of Kinematic Controller
    vc = vr*cos(qe(3,:)) + K(1,:)*qe(1,:);
    wc = wr + vr*K(2, :)*qe(2,:) + vr*K(3, :)*sin(qe(3,:));
    Vref = [vc; wc];

    Phi_ref = Omega*Vref;

    %% Dynamic Control 
    a_dyn = Kp_dyn*(Phi_ref - Phi);
    tau = inv_E_bar*M_bar*a_dyn;
    
    %% Dynamic
    Phi_dot = inv_M_bar*E_bar*tau;
  
    %% Integrator
    th = q(3);
    Phi = Phi + Phi_dot*tsteps;
    q_dot = S(th)*Phi;
    q = q + q_dot*tsteps;
    V = pinv(Sv(th))*q_dot;
    
    %% Vectors
    tvec = [tvec, t];
    V_vec = [V_vec , V];
    Vr_vec = [Vr_vec , Vref];
    Phi_vec = [Phi_vec , Phi];
    q_vec = [q_vec , q];
    q_dot_vec = [q_dot_vec , q_dot]; 
    qr_vec = [qr_vec, qr];
    Vr_vec = [Vr_vec, Vref];
    qe_vec = [qe_vec, qe];
    tau_vec = [tau_vec, tau];
    
    %%Upadte time and iteration
    t = t + tsteps;
    iteration = iteration + 1;

end

%% Plot
%Plot desired and real trajectory
fig1 = figure('Name','Desired and Real trajectory for robot','NumberTitle','off');
hold on
title('Desired and Real trajectory in X-Y plane');
grid on
xlabel({'X [m]'});
ylabel({'Y [m]'});
plot(qr_vec(1, :), qr_vec(2, :), 'b-');
hold on
plot(q_vec(1, :), q_vec(2, :), 'r--');
legend({'Reference path' , 'Real path'}, 'Location', 'northeast')

%Plot X error of robot
fig2 = figure('Name','Error of X position of robot','NumberTitle','off');
hold on
title('Error of X for robot');
grid on
xlabel({'Time [s]'});
ylabel({'Error [m]'});
plot(tvec, qe_vec(1, :), 'b');
grid minor

%Plot Y error of robot
fig3 = figure('Name','Error of Y position of robot','NumberTitle','off');
hold on
title('Error of Y for robot');
grid on
xlabel({'Time [s]'});
ylabel({'Error [m]'});
plot(tvec, qe_vec(2, :), 'b');
grid minor

%Plot Theta error of robot
fig4 = figure('Name','Error of Theta position of robot','NumberTitle','off');
hold on
title('Error of Theta for robot');
grid on
xlabel({'Time [s]'});
ylabel({'Error [m]'});
plot(tvec, qe_vec(3, :), 'b');
grid minor

%Plot Torque of wheels
fig5 = figure('Name','Torques of wheels','NumberTitle','off');
hold on
title('Torques of wheels');
grid on
xlabel({'Time [s]'});
ylabel({'Torque [N.m]'});
plot(tvec, tau_vec(1, :), 'b');
hold on
plot(tvec, tau_vec(2, :), 'r');
legend({'torque of rigth wheel' , 'torque of left wheel'}, 'Location', 'northeast')

