function Xdot=ode_augmented_pronav_plant_ctrllr(~,X)
global Ap Bp Cp Ac Bc1 Bc2 Cc Dc1 Dc2 Np nT 
% Augmented State Vector(X):
% X =[X_PN; Xc ; Xp];
% WHERE 
% ProNav_State_vector = [beta; RT1; RT2; RM1; RM2; VT1; VT2; VM1; VM2];
% Controller State,Xc is the integral error of A_z 
% Plant state vector,Xp=[alpha; q; dele; dele_dot]

% Define pointers to state variables
%--------------------------------------------------------------------------
% Pointers to states
%ProNav States:
sel_beta = 1;
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;
%Controller State:
sel_Xc   = 10;
%Plant States:
sel_alpha= 11;
sel_q    = 12; 
sel_dele = 13; 
sel_dele_dot=14;

% target and missile velocity magnitudes
VT = sqrt( X(sel_VT1)^2 + X(sel_VT2)^2 );

% relative positions and velocities
RTM1 = X(sel_RT1) - X(sel_RM1);
RTM2 = X(sel_RT2) - X(sel_RM2);
VTM1 = X(sel_VT1) - X(sel_VM1);
VTM2 = X(sel_VT2) - X(sel_VM2);
 
% relative distance
RTM = sqrt(RTM1^2 + RTM2^2);

% line of sight angle and time derivative
lambda     = atan2( RTM2, RTM1 );
lambda_dot = (RTM1*VTM2 - RTM2*VTM1)/RTM^2;

% closing velocity
VC = -(RTM1*VTM1 + RTM2*VTM2)/RTM;

% command acceleration (True Proportional Navigation)
A_Z_cmd = Np*VC*lambda_dot;% Reference Command that A_z has to track.

Xc=X(sel_Xc);  %Controller State
Xp=X(sel_alpha:sel_dele_dot);%Plant State
Yp=Cp*Xp; %Plant Output, Note that Dp=0 
A_z_actual = Yp(1);

%% Dynamics: 

%ProNav Dynamics:
X_PNdot(1,1) =  nT/VT;
X_PNdot(2,1) = -VT*cos(X(sel_beta));
X_PNdot(3,1) =  VT*sin(X(sel_beta));
X_PNdot(4,1) =  X(sel_VM1);
X_PNdot(5,1) =  X(sel_VM2);
X_PNdot(6,1) =  nT*sin(X(sel_beta));
X_PNdot(7,1) =  nT*cos(X(sel_beta));
X_PNdot(8,1) = -A_z_actual*sin(lambda);
X_PNdot(9,1) =  A_z_actual*cos(lambda);

%Controller Dynamics:
Xcdot = Ac*Xc + Bc1*Yp + Bc2*A_Z_cmd;
dele_command = Cc*Xc + Dc1*Yp + Dc2*A_Z_cmd;

% % Plant Dynamics (Actuator Incorporated)
Xpdot= Ap*Xp + Bp*dele_command;

% Augmented Dynamics:
Xdot = [X_PNdot;Xcdot;Xpdot];
