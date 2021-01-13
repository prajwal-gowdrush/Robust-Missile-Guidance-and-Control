function Xdot=ode_augmented_plant_ctrllr_observer(~,X)
global Ap Bp Cp A_comp B1_comp B2_comp C_comp D1_comp D2_comp
% Augmented State Vector(X):
% X =[X_comp ; Xp];
% WHERE 
% Compensator State Vector, X_comp = [eIAz; eIAz_hat; alpha_hat; q_hat];
% Plant state vector, Xp=[alpha; q; dele; dele_dot]

% Define pointers to state variables
%--------------------------------------------------------------------------
%%Pointers to states
%Compensator States:
sel_eIAz      =1;
sel_eIAz_hat  =2;
sel_alpha_hat =3;
sel_q_hat     =4; 
%Plant States:
sel_alpha     =5;
sel_q         =6; 
sel_dele      =7; 
sel_dele_dot  =8;

X_comp=X(sel_eIAz:sel_q_hat);%Compensator State
Xp=X(sel_alpha:sel_dele_dot);%Plant State
Yp=Cp*Xp; %Plant Output, Note that Dp=0 

A_Z_cmd=32.17;% Reference Command that A_z has to track.

%Compensator Dynamics:
X_comp_dot  =  A_comp*X_comp + B1_comp*Yp + B2_comp*A_Z_cmd;
dele_command = C_comp*X_comp + D1_comp*Yp + D2_comp*A_Z_cmd;

% Plant Dynamics
Xp_dot= Ap*Xp + Bp*dele_command;

% Augmented Dynamics:
Xdot = [X_comp_dot;Xp_dot];