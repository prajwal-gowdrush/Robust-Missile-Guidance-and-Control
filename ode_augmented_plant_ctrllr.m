function Xdot=ode_augmented_plant_ctrllr(~,X)
global Ap Bp Cp Ac Bc1 Bc2 Cc Dc1 Dc2
% Augmented State Vector Y (Includes Controller State)
% X =[Xc ; alpha; q; dele; dele_dot];
% where x_c is the integral error of A_z 
% augmented to the plant state vector Xp=[alpha; q; dele; dele_dot]

Xc=X(1);  %Controller State
Xp=X(2:5);%Plant State
Yp=Cp*Xp; %Plant Output, Note that Dp=0 

r=32.17;% Reference Command that A_z has to track.

%Controller Dynamics:
Xcdot = Ac*Xc + Bc1*Yp + Bc2*r;
dele_command = Cc*Xc + Dc1*Yp + Dc2*r;

% Plant Dynamics
Xpdot= Ap*Xp + Bp*dele_command;

% Augmented Dynamics:
Xdot = [Xcdot;Xpdot];