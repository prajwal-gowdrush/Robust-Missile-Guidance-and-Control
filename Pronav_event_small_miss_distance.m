function [value,isterminal,direction] = Pronav_event_small_miss_distance(t,y)
global HE_rad Np nT
% Define pointers to state variables
%--------------------------------------------------------------------------
% Pointers to states
sel_beta = 1;
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;
    % target and missile velocity magnitudes
    VT = sqrt( y(sel_VT1).^2 + y(sel_VT2).^2 );

    % relative positions and velocities
    RTM1 = y(sel_RT1) - y(sel_RM1);
    RTM2 = y(sel_RT2) - y(sel_RM2);
    VTM1 = y(sel_VT1) - y(sel_VM1);
    VTM2 = y(sel_VT2) - y(sel_VM2);

    % relative distance
    RTM = sqrt(RTM1.^2 + RTM2.^2);

    % line of sight angle and time derivative
    lambda     = atan2( RTM2, RTM1 );
    lambda_dot = (RTM1.*VTM2 - RTM2.*VTM1)./RTM.^2;

    % closing velocity
    VC = -(RTM1.*VTM1 + RTM2.*VTM2)./RTM;

    % Compute acc commands
    nc = Np*VC.*lambda_dot;
    
    
%Value that takes on 0 instructing the integrator to stop 
value = (nc/32.2 - (-6.27))>0;
isterminal = 1;  % Halt integration 
direction = 0;   