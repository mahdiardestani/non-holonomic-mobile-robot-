close all
clear
clc

%% Global Variables
global r a d mc mw mt Im Ic Iw I M_bar C_bar E_bar

%Robot Variables
r = 0.1; %radius of wheels in cm
a = 0.3; %distance between center of chasi of wheels to wheels
d = 0; %distance between of center of mass and center of chasi
mc = 3; %mass of the robot without wheels and motors in[kg]
mw = 1; %mass of wheel in[kg]
mt = mc + 2*mw;
Im = mc*a^2;    %moment of inertia of robot
Ic = 11.4180;
Iw = 0.0456;
I = mc*d^2 + Ic + 2*mw*(d^2+a^2) + 2*Im;
%%=========================================%%

%% Euler-lagrange matrix
M_bar = [Iw + r^2*(mt/4 - (mt*d^2)/(4*a^2) + I/(4*a^2)), r^2*(mt/4 + (mt*d^2)/(4*a^2) - I/(4*a^2));
     r^2*(mt/4 + (mt*d^2)/(4*a^2) - I/(4*a^2)), Iw + r^2*(mt/4 - (mt*d^2)/(4*a^2) + I/(4*a^2))];
C_bar = zeros(2);
E_bar = eye(2);
%Jacobian Matrix: transform [v, w].' to [phi_r, phi_l]
Omega = [1/r a/r; 1/r -a/r];

%%=========================================%%
%% Vectors
qr_vec = [];
qe_vec = [];
q_vec = [];
N = 1; %number of agent
q = zeros(3, N);
q_vec(:, 1) = [q_vec, q];
Vref_vec = [];
tau_vec = [];
Phi_vec = [];
Phi = [0; 0];
Phi_vec = [Phi_vec, Phi];

%%=========================================%%
Vref = [1; 1];  %[vr, wr].' in (m/s)
vr = Vref(1, :);
wr = Vref(2, :);

%%=========================================%%
%% Simulation Time
simul_time = 20;
tsteps = 0.001;
time_sim = 0:tsteps:simul_time;

%%=========================================%%
%% PD Dynamic Controller
kp_dyn = 20;

%Kinematic Controller
% K = [10; 10; 5]; %control parameters
K = [5; 10; 2];

%% Simulation time
iteration = 1;
t = tsteps;
tvec = [];

%% Simulation

while t <= simul_time
    
    tvec = [tvec, t];

    %% Reference Trajectory
    [qr] = reference_trajectory(t, vr, wr);
    qr_vec = [qr_vec, qr];
    
    %% Kinematic Control
    [Vref, qe] = kinematic_controller(qr, q_vec(:, iteration), K, vr, wr);
    qe_vec = [qe_vec, qe];


    %% Dynamic Controller
    [Phi_dot, tau] = dynamic_controller(Vref, Phi, Omega, kp_dyn);
    Vref_vec = [Vref_vec, Vref];
    tau_vec = [tau_vec, tau];
    
    %% Agent
    [q, q_dot, Phi] = agent(q_vec(:, iteration), Phi_vec(:, iteration), Phi_dot, tsteps);
    q_vec = [q_vec, q];
    Phi_vec = [Phi_vec, Phi];
    
    %% Update iteration
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

%% Functions
%reference_trajectory
function [qr] = reference_trajectory(t, vr, wr)

    r_path = vr/wr;
    qr = [r_path*cos(wr*t); r_path*sin(wr*t); wr*t+pi/2];
    
end

%Kinematic Controller
function [Vref, qe] = kinematic_controller(qr, q, K, vr, wr)

    theta = q(3, :);
    rotaion_matrix = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1]; %from inertia to local
    
    xe = qr(1, :) - q(1, :);
    ye = qr(2, :) - q(2, :);
    thetae = qr(3, :) - q(3, :);
    qe_temp = [xe; ye; thetae];
    qe = rotaion_matrix*qe_temp;
    vc = vr*cos(qe(3,:)) + K(1,:)*qe(1,:);
    wc = wr + vr*K(2, :)*qe(2,:) + vr*K(3, :)*sin(qe(3,:));
    Vref = [vc; wc];
    
end

%Dynamic Controller
function [Phi_dot, tau] = dynamic_controller(Vref, Phi, Omega, Kp_dyn)
    
    global  M_bar E_bar
    Phi_ref = Omega*Vref;
    a_dyn = Kp_dyn*(Phi_ref - Phi);
    tau = inv(E_bar)*M_bar*a_dyn;
    Phi_dot = inv(M_bar)*E_bar*tau;
    
end

%Agent
function [q, q_dot, Phi] = agent(q, Phi, Phi_dot, tsteps)

    global r a d
    theta = q(3);
    S = @(theta)[(r/(2*a))*(a*cos(theta) - d*sin(theta)), (r/(2*a))*(a*cos(theta) + d*sin(theta));
    (r/(2*a))*(a*sin(theta) + d*cos(theta)), (r/(2*a))*(a*sin(theta)-d*cos(theta));
    r/(2*a), -r/(2*a);];

    Phi = Phi + Phi_dot*tsteps;
    q_dot = S(theta)*Phi;
    q = q + q_dot*tsteps;

end


