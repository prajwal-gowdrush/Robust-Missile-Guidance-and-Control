function [value,isterminal,direction] = obs_event_small_miss_distance(t,y)

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
%Compensator States:
sel_eIAz      =10;
sel_eIAz_hat  =11;
sel_alpha_hat =12;
sel_q_hat     =13; 
%Plant States:
sel_alpha     =14;
sel_q         =15; 
sel_dele      =16; 
sel_dele_dot  =17;
% relative positions 
RTM1 = y(sel_RT1) - y(sel_RM1);
RTM2 = y(sel_RT2) - y(sel_RM2);

% relative distance
RTM = sqrt(RTM1^2 + RTM2^2);
dele = y(sel_dele);
dele_dot = y(sel_dele_dot);
eIAz = y(sel_eIAz);

stop_good_enough_miss_distance = 0;% (RTM<0.5);
stop_dele_maxed_out = (t>7)&&(abs(180/pi*dele)>30); 
stop_dele_dot_maxed_out = (t>7)&&(abs(180/pi*dele_dot)>30);
stop_eIAz_blowing_up = (t>7)&&(abs(eIAz)>50);

%Value that takes on 0 instructing the integrator to stop if any of the above conditions are met: 
value = ~(stop_good_enough_miss_distance + stop_dele_maxed_out + stop_dele_dot_maxed_out + stop_eIAz_blowing_up);
isterminal = 1;  % Halt integration 
direction = -1;