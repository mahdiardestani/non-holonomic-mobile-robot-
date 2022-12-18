%Calculate Euler_lagrange matrixs
clc;
clear;
close all;
%%
syms a r theta(t) t mt I Iw real

%Constraint matrix with c != p
% S = [(r/(2*a))*(a*cos(theta) - d*sin(theta)), (r/(2*a))*(a*cos(theta) + d*sin(theta));
%     (r/(2*a))*(a*sin(theta) + d*cos(theta)), (r/(2*a))*(a*sin(theta)-d*cos(theta));
%     r/(2*a), -r/(2*a);
%     1 , 0;
%     0 , 1];
d=0;
S = [(r/(2*a))*(a*cos(theta) - d*sin(theta)), (r/(2*a))*(a*cos(theta) + d*sin(theta));
    (r/(2*a))*(a*sin(theta) + d*cos(theta)), (r/(2*a))*(a*sin(theta)-d*cos(theta));
    r/(2*a), -r/(2*a);];

%Constraint matrix with c = p
% S = 1/2*[r*cos(theta(t)) r*cos(theta(t));
%         r*sin(theta(t)) r*sin(theta(t));
%         r/a -r/a;
%         2 0;
%         0 2];

%Inertia matrix
% H = [mt 0 mt*d*sin(theta(t)) 0 0;
%     0 mt -mt*d*cos(theta(t)) 0 0;
%     mt*d*sin(theta(t)) -mt*d*cos(theta(t)) I 0 0;
%     0 0 0 Iw 0;
%     0 0 0 0 Iw];

H = [mt 0 mt*d*sin(theta(t));
    0 mt -mt*d*cos(theta(t));
    mt*d*sin(theta(t)) -mt*d*cos(theta(t)) I;];

H_bar = S.'*H*S;
H_dot = diff(H,t);
simplify(H_bar);

% C = [0 0 mt*d*diff(theta(t),t)*cos(theta(t)) 0 0;
%     0 0 mt*d*diff(theta(t),t)*sin(theta(t)) 0 0;
%     0 0 0 0 0;
%     0 0 0 0 0;
%     0 0 0 0 0];
C = [0 0 mt*d*diff(theta(t),t)*cos(theta(t));
    0 0 mt*d*diff(theta(t),t)*sin(theta(t));
    0 0 0;];

S_dot = diff(S, t);
C_bar = S.'*H*S_dot + S.'*C*S;
simplify(C_bar);

% E = [0 0;
%     0 0;
%     0 0;
%     1 0;
%     0 1];

E = 1/r*[cos(theta) cos(theta);sin(theta) sin(theta);a -a];

E_bar = S.'*E;
simplify(E_bar);

%Skew symmetric
H_bardot = S_dot.'*H*S + S.'*H_dot*S + S.'*H*S_dot;
simplify(H_bardot);