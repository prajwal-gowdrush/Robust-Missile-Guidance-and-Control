%% Mainscript Project: Controller Design
clc; clear; close all;
global Ap Bp Cp Ac Bc1 Bc2 Cc Dc1 Dc2 HE_rad Np nT
V0 = 2813.32; %(Trimmed Airspeed in ft/s)(= 2.5 Mach)

% Parameters (ProNav)
nT = 3*32.2;%Evasive 3g acceleration normal to it's velocity
HE_rad = -20*pi/180;
Np = 3.4; %Typically ranges between 3 to 5

%Linearized Dynamics of the Missile(Pursuer) (MRAAM):
A = [-1.57, 0 , 0, 1, 0; 0, -0.5, 0.17, 0, -1; -21.13, -2876.7, -2.10, -0.14, -0.05; -82.92, -11.22, -0.01, -0.57, 0; -0.19, -11.86, -0.01, 0, -0.57];
B = [0, -0.1, 0; -0.07, 0, 0.11; -1234.7, -30.49, -1803.2; -4.82, -119.65, -7; 14.84, 0.27, -150.58];
open_loop_poles_A = eig(A);
%% Controller-Design
%% Extracted Short-Period Dynamics:
% State = [alpha,q]; , Control=delta_e;
A_sp = A([1,4],[1,4]);
B_sp = B([1,4],2);
Z_alpha = A_sp(1,1)* V0;% Aerodynamic coefficient
Z_dele =  B_sp(1,1)*V0;% Aerodynamic coefficient
C_sp = [Z_alpha 0 ; 1 0; 0 1]; 
D_sp = [Z_dele; 0 ; 0];
sys_sp = ss(A_sp,B_sp,C_sp,D_sp);

%% Actuator Dynamics:
w_n = 35*2*pi; %Natural Frequency (in radians per second)
z_damp = 0.71; %Damping Factor
A_act = [0, 1; -w_n^2, -2*z_damp*w_n];
B_act = [0; w_n^2];
C_act = [1,0];
D_act = 0;
sys_act = ss(A_act,B_act,C_act,D_act);

%% Extended system OR Control design model:
C_reg = [Z_alpha 0]; 
D_reg = Z_dele; 
A_tilda = [0 C_reg; zeros(size(A_sp,1),1) A_sp];
B_tilda = [D_reg; B_sp];
B_command = [-1; zeros(size(B_sp))];

%% Plant Dynamics:(incorporates the actuator model)
sys_plant = series(sys_act,sys_sp);
[Ap,Bp,Cp,Dp]=ssdata(sys_plant);

%% Index of Design Metrics for Controller-Design:
%1=Gain Margin at Plant Input
%2=Phase Margin at Plant Input
%3=Maximum Elevator displacement (in degrees)
%4=Maximum Elevator Rate(in degrees per second)
%5=(1/3*w_n - Wcp) (should be greater than zero)
    % , where Wcp = Gain Cross-Over Frequency of Lu , w_n = actuator natural frequency
%6=Undershoot percentage
%7=Overshoot percentage
%8=Max Sensitivity over 1e-1 to 1e3 rad/s
%9=Max Co-Sensitivity over 1e-1 to 1e3 rad/s
%10=Minimum Closed-Loop Eigenvalues
%11=Miss Distance

%% LQR Penalty choices with corresponding time and frequency responses.
qq =logspace(-5,-2,50);
R = 1;
qqlen = length(qq);

%%INITIALIZING
does_design_meet_specs = false(qqlen,10); % Initializing a logical array
DesignMetrics = zeros(qqlen,10);%Initializing Design Metrics Matrix
percentage_satisfaction_specs = zeros(qqlen,10);%Initializing Percentage Satisfaction Matrix
Gain_CrossOverFreqs = zeros(qqlen,1);
closed_loop_poles = zeros(qqlen,5);
% Will be used later to recognize which designs meet the specifications
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Full-State Feedback - Controller Design - Running through test values of qq');
for h=1:qqlen
    
Q = [qq(h), 0, 0; 0, 0, 0; 0, 0, 0];
[K_lqr,~,~] = lqr(A_tilda,B_tilda,Q,R,0);

%%Controller Dynamics:(for each LQR penalty choice)
Ac=0;  Bc1=[1 0 0];  Bc2= -1;
Cc = -K_lqr(1); Dc1 = [0,-K_lqr(1,2:3)]; Dc2=0;
sys_ctrllr = ss(Ac,Bc1,Cc,Dc1);
sys_ctrllr_tocomputeloopgain = ss(Ac,Bc1,-Cc,-Dc1);

%%Closed Loop System
yp2yclosed = [1,0,0];%Picking up only A_z as the plant output.
[Aclosed, Bclosed, Cclosed, Dclosed]=closedloop_from_plantplusctrler(Ap,Bp,Cp,Dp,Ac,Bc1,Bc2,Cc, Dc1, Dc2, yp2yclosed);
sys_closedloop = ss(Aclosed, Bclosed, Cclosed, Dclosed);
closed_loop_poles(h,:) = eig(Aclosed)';

%Step Response related design metrics
StepRes= stepinfo(sys_closedloop);

%Loop Gain Model at Plant Input:
Lu_ss=series(sys_plant,sys_ctrllr_tocomputeloopgain); 

%Loop Gain Model at Plant Output:
Ly_ss=series(sys_ctrllr_tocomputeloopgain,sys_plant); 

%Loop Gain Lu : Frequency response
% |w| = frequencies at which frequency response of Lu is evaluated
w=1i*logspace(-1,3,1000);
Lu = freqresp(Lu_ss,w); 
Lu=squeeze(Lu(1,:,:)); %Loop Gain Model Frequency Response at Plant Input
Ly = freqresp(Ly_ss,w);
Ly_Az_Az = squeeze(Ly(1,1,:));%A_z channel Ly frequency response
Ly_q_q = squeeze(Ly(3,3,:)); % q channel Ly frequency response
S_Az_Az = 1./(1 + Ly_Az_Az);  %Sensitivity in the A_z channel
T_Az_Az = Ly_Az_Az./(1 + Ly_Az_Az);%Co-Sensitivity in the A_z channel
S_q_q = 1./(1 + Ly_q_q); %Sensitivity in the q channel
T_q_q = Ly_q_q./(1 + Ly_q_q); %Co-Sensitivity in the q channel
MaxSen_Az = max(abs(S_Az_Az));
MaxCoSen_Az = max(abs(T_Az_Az));
MaxSen_q = max(abs(S_q_q));
MaxCoSen_q = max(abs(T_q_q));

MaxSen = max(MaxSen_Az,MaxSen_q);
MaxCoSen = max(MaxCoSen_Az,MaxCoSen_q);


RD = 1+Lu; %Return Differences at frequencies w
[~,MinRD_index] = min(abs(RD));
[~,PM,~,Wcp] = margin(Lu_ss);
%As can be seen in the Nyquist plot there are multiple phase cross-over
%frequencies and therefore multiple Gain Margins: The positive Gain Margin
%is of interest to us in this particular control system design. 
%So, I've ignored the GM and Wcg that the margin MATLAB function gives us.
%We'll have to compute them ourselves in a different way:
Lu_selected=Lu(real(Lu)<-0.05);
Lu_selected=Lu_selected(real(Lu_selected)>-1);
%We're limiting our search to Lu's with real parts between -1 and -0.05.
%(This is by inspection of Nyquist plots)
Lu_real_GM = interp1(imag(Lu_selected),real(Lu_selected),0);
GM = -1/Lu_real_GM;

%% ProNav+Autopilot Simulations to help us find the Max Elevator Displacement and Max Elevator Rate

%TimeSpan:
tspan=[0,11];

% Define pointers to state variables
%--------------------------------------------------------------------------
% Pointers to states
%ProNav State Pointers:
sel_beta = 1; 
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;
%Controller State(eI_A_Z) Pointer :
sel_Xc   = 10; 
%Plant State Pointers:
sel_alpha= 11;
sel_q    = 12; 
sel_dele = 13; 
sel_dele_dot=14;

% ProNav Initial Conditions
    beta_rad = 0;
    RT1      = 20000;%[ft]
    RT2      = 40000;%[ft]
    RM1      = 0;%[ft]
    RM2      = 40000;%[ft]
    VM       = 2813.32; %[ft/s]
    VT       = 300;%[ft/s]
    VT1      = -VT*cos(beta_rad);%[ft/s]
    VT2      =  VT*sin(beta_rad);%[ft/s]
    
% relative positions and velocities
    RTM1 = RT1 - RM1;
    RTM2 = RT2 - RM2;

    % relative distance
    RTM = sqrt(RTM1^2 + RTM2^2);

    % line of sight angle and time derivative
    lambda = atan2( RTM2, RTM1 );

    % missile lead angle
    L = asin( VT*sin( beta_rad + lambda )/VM );

    % missile velocity components
    VM1  = VM*cos(lambda + L + HE_rad);
    VM2  = VM*sin(lambda + L + HE_rad);
    
% Initial ProNav State vector
X_PN0 = [beta_rad, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]';

%Initial Controller State
Xc0 = 0;

%Initial Plant State
Xp0=[0;0;0;0];
    
%%Augmented Initial State Vector
X0 =[X_PN0;Xc0;Xp0];
options = odeset('AbsTol',1e-6,'RelTol',1e-6, 'Events',@event_small_miss_distance);
[time,X_Sol,time_event,X_Sol_event,ie]=ode45(@ode_augmented_pronav_plant_ctrllr,tspan,X0,options);

%Extracting ProNav States from the ode45 solution matrix:
  % target and missile velocity magnitudes
    VT = sqrt( X_Sol(:,sel_VT1).^2 + X_Sol(:,sel_VT2).^2 );

    % relative positions and velocities
    RTM1 = X_Sol(:,sel_RT1) - X_Sol(:,sel_RM1);
    RTM2 = X_Sol(:,sel_RT2) - X_Sol(:,sel_RM2);
    VTM1 = X_Sol(:,sel_VT1) - X_Sol(:,sel_VM1);
    VTM2 = X_Sol(:,sel_VT2) - X_Sol(:,sel_VM2);

    % relative distance
    RTM = sqrt(RTM1.^2 + RTM2.^2);
    MissDistance = min(RTM);
    
    % line of sight angle and time derivative
    lambda     = atan2( RTM2, RTM1 );
    lambda_dot = (RTM1.*VTM2 - RTM2.*VTM1)./RTM.^2;

    % closing velocity
    VC = -(RTM1.*VTM1 + RTM2.*VTM2)./RTM;

    % Compute acceleration commands
    A_z_cmd = Np*VC.*lambda_dot;

%Extracting Controller and Plant States from the ode45 solution matrix: 
Xc = X_Sol(:,sel_Xc);
alpha = X_Sol(:,sel_alpha);
q = X_Sol(:,sel_q);
dele = X_Sol(:,sel_dele);
dele_dot = X_Sol(:,sel_dele_dot);
A_z_actual = Z_alpha*alpha + Z_dele*dele; 
MaxElevatorDisplacement_degrees = (180/pi)*max(abs(dele));
MaxElevatorRate_degrees_per_second = (180/pi)*max(abs(dele_dot));

%% Plots
% A_z command tracking during engagement:
subplot(1,2,1)
ax1=gca;
% hold off;
plot(ax1,time,A_z_actual./32.17,'b','LineWidth',1.2); hold on
plot(ax1,time,A_z_cmd./32.17,'r--','LineWidth',1);
xlabel(ax1,'Time'); str1=ylabel(ax1,'$A_{z_{cmd}}$,$A_{z_{actual}}$ (in G''s)');
set(str1,'Interpreter','latex');
str2=title(ax1,'Tracking of $A_{z_{cmd}}$ by Missile''s $A_{z_{actual}}$'); grid(ax1,'on');
set(str2,'Interpreter','latex');
ax1.GridLineStyle = ':'; ax1.GridAlpha = 0.3; ax1.FontSize = 16; ax1.LineWidth = 1.4;
xlim([0,11]); ylim([-10 20]); axis square;
leg=legend('$A_{z_{actual}}$','$A_{z_{cmd}}$'); set(leg,'Interpreter','latex');

subplot(1,2,2)
ax2=gca;
plot(ax2,real(Lu),imag(Lu),'b','LineWidth',1); hold on;
%Plotting the return Difference Vector
plot([-1 real(Lu(MinRD_index))],[0 imag(Lu(MinRD_index))], 'r','LineWidth',0.5);
hold on; plot(-1,0,'r*','LineWidth',1.2);
xlabel(ax2,'Real Axis'); ylabel(ax2,'Imaginary Axis');
title(ax2,'A portion of Nyquist Plot of L_{u_{LQR}}'); grid(ax2,'on');
ax2.GridLineStyle = ':'; ax2.GridAlpha = 0.3; ax2.FontSize = 16; ax2.LineWidth = 1.4;
xlim([-4 4]); ylim([-4 4]); axis square;
circle(-1,0,1);
pause(0.01);

%%DESIGN METRIC VALUES FOR CURRENT qq(h)
DesignMetrics(h,1)=20*log10(GM);%in Decibels
DesignMetrics(h,2)=PM;
DesignMetrics(h,3)=MaxElevatorDisplacement_degrees;
DesignMetrics(h,4)=MaxElevatorRate_degrees_per_second;
DesignMetrics(h,5)=Wcp;% in radians per second
DesignMetrics(h,6)=StepRes.Undershoot;
DesignMetrics(h,7)=StepRes.Overshoot;
DesignMetrics(h,8)= 20*log10(MaxSen);   %in Decibels
DesignMetrics(h,9)= 20*log10(MaxCoSen); %in Decibels
DesignMetrics(h,10)=min(real(eig(Aclosed)));
DesignMetrics(h,11)=MissDistance;
Gain_CrossOverFreqs(h)= Wcp;
%% Recognizing which design specifications are met:
%Bounds for Barely Meeting the specifications:
LB1 = 6;
LB2 = 35;
UB3 = 35;
UB4 = 350;
UB5 = (1/3*w_n);
UB6 = 10;
UB7 = 10;
UB8 = 6;
UB9 = 6;
LB10 = -500;

does_design_meet_specs(h,1) =DesignMetrics(h,1)>=LB1;
does_design_meet_specs(h,2) =DesignMetrics(h,2)>=LB2;
does_design_meet_specs(h,3) =DesignMetrics(h,3)<=UB3;
does_design_meet_specs(h,4) =DesignMetrics(h,4)<=UB4;
does_design_meet_specs(h,5) =DesignMetrics(h,5)<UB5;
does_design_meet_specs(h,6) =DesignMetrics(h,6)<UB6 ;
does_design_meet_specs(h,7) =DesignMetrics(h,7)<UB7;
does_design_meet_specs(h,8) =DesignMetrics(h,8)<UB8;
does_design_meet_specs(h,9) =DesignMetrics(h,9)<UB9;
does_design_meet_specs(h,10) =DesignMetrics(h,10)>LB10;

percentage_satisfaction_specs(h,1)= 100*(DesignMetrics(h,1)-LB1)/abs(LB1);
percentage_satisfaction_specs(h,2)= 100*(DesignMetrics(h,2)-LB2)/abs(LB2);
percentage_satisfaction_specs(h,3)= 100*(-DesignMetrics(h,3)+UB3)/abs(UB3);
percentage_satisfaction_specs(h,4)= 100*(-DesignMetrics(h,4)+UB4)/abs(UB4);
percentage_satisfaction_specs(h,5)=100*(-DesignMetrics(h,5)+UB5)/abs(UB5);
percentage_satisfaction_specs(h,6)= 100*(-DesignMetrics(h,6)+UB6)/abs(UB6);
percentage_satisfaction_specs(h,7)= 100*(-DesignMetrics(h,7)+UB7)/abs(UB7);
percentage_satisfaction_specs(h,8)= 100*(-DesignMetrics(h,8)+UB8)/abs(UB8);
percentage_satisfaction_specs(h,9)= 100*(-DesignMetrics(h,9)+UB9)/abs(UB9);
percentage_satisfaction_specs(h,10)=100*(DesignMetrics(h,10)-LB10)/abs(LB10);

end

%To make sure we don't need heavy Loop Transfer Recovery to meet design
%specifications during observer design:
min_satisfaction_percent_2qualify = 10; 

min_percentage_satisfaction_specs = min(percentage_satisfaction_specs')';
qq_qualifying = qq(min_percentage_satisfaction_specs>min_satisfaction_percent_2qualify);
MissDistances_qualifying = DesignMetrics(min_percentage_satisfaction_specs>min_satisfaction_percent_2qualify,11);

%Picking the Qualifying Design with minimum miss distance:
minMissDistance=min(MissDistances_qualifying);
winning_design_index = find(DesignMetrics(:,11)==minMissDistance); 

%So, the best LQR penalty for the integral error among the ones tried out
%is:
qq_design=qq(winning_design_index);

%% Design Charts :  
% The code has already picked out the best penalty matrix for us in terms of minimizing the rise time while meeting the design specifications. 
% Just for the sake of visual confirmation. We plot some design charts
% below:
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(Gain_CrossOverFreqs,DesignMetrics(:,1),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,1),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[LB1,LB1],'g--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); ylabel('GM at Plant Input [dB]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('GM(dB)','Design Point','Design above line'); set(leg,'Interpreter','latex');

subplot(2,2,2);
plot(Gain_CrossOverFreqs,DesignMetrics(:,2),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,2),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[LB2,LB2],'g--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); ylabel('PM at Plant Input [degrees]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('PM(degrees)','Design Point','Design above line'); set(leg,'Interpreter','latex');

subplot(2,2,3);
plot(Gain_CrossOverFreqs,DesignMetrics(:,3),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,3),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB3,UB3],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); str=ylabel('Max $\delta_e$(degrees)'); grid on;
set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('$\delta_e$(degrees)','Design Point','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,4);
plot(Gain_CrossOverFreqs,DesignMetrics(:,4),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,4),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB4,UB4],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); str=ylabel('Max $\dot{\delta_e}$(in degrees/sec)'); grid on;
set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('$\dot{\delta_e}$','Design Point','Design below line'); set(leg,'Interpreter','latex');

%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(Gain_CrossOverFreqs,DesignMetrics(:,5),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,5),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB5,UB5],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies [rad/s]'); str=ylabel('$W_{cp}$ [rad/s] '); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('$W_{cp}$ [rad/s]','Design Point','$\frac{1}{3}\omega_n$ [rad/s] '); set(leg,'Interpreter','latex');

subplot(2,2,2);
plot(Gain_CrossOverFreqs,DesignMetrics(:,6),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,6),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB6,UB6],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); ylabel('Percentage Undershoot'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('% Undershoot','Design Point','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,3);
plot(Gain_CrossOverFreqs,DesignMetrics(:,7),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,7),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB7,UB7],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); ylabel('Percentage Overshoot'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
ylim([0,12]);
leg=legend('% Overshoot','Design Point','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,4);
plot(Gain_CrossOverFreqs,DesignMetrics(:,8),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,8),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB8,UB8],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); str=ylabel('Max Sensitivity [dB]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
ylim([0,8]);
legend('Max S in A_z and q channels [dB]','Design Point','Design below line');

%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(Gain_CrossOverFreqs,DesignMetrics(:,9),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,9),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[UB9,UB9],'k--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); str=ylabel('Max Co-Sensitivity [dB]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
ylim([0,8]);
legend('Max T in A_z and q channels [dB]','Design Point','Design below line');

subplot(2,2,2);
plot(Gain_CrossOverFreqs,DesignMetrics(:,10),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,10),'ro','LineWidth',1.5);
plot([Gain_CrossOverFreqs(1) Gain_CrossOverFreqs(end)],[LB10,LB10],'g--','LineWidth',1);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); str=ylabel('Minimum Re$(\lambda_{cl})$'); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('Minimum Re($\lambda_{cl}$','Design Point','Design above line'); set(leg,'Interpreter','latex');

subplot(2,2,3);
plot(Gain_CrossOverFreqs,DesignMetrics(:,11),'b','LineWidth',1.5);hold on;
plot(Gain_CrossOverFreqs(winning_design_index),DesignMetrics(winning_design_index,11),'ro','LineWidth',1.5);
xlabel('Gain Cross-Over Frequencies (in rad/s)'); str=ylabel('Miss Distance [ft]'); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([Gain_CrossOverFreqs(1),Gain_CrossOverFreqs(end)]);
leg=legend('Miss Distance [ft]','Design Point'); set(leg,'Interpreter','latex');

%% Recomputing time and frequency responses for the design penalty matrix

% Frequency responses for a random penalty matrix
Q = [qq_design, 0, 0; 0, 0, 0; 0, 0, 0];
R = 1;
[K_lqr,~,~] = lqr(A_tilda,B_tilda,Q,R,0);

%%Controller Dynamics:(for each LQR penalty choice)
Ac=0;  Bc1=[1 0 0];  Bc2= -1;
Cc = -K_lqr(1); Dc1 = [0,-K_lqr(1,2:3)]; Dc2=0;
sys_ctrllr = ss(Ac,Bc1,Cc,Dc1);
sys_ctrllr_tocomputeloopgain = ss(Ac,Bc1,-Cc,-Dc1);

%%Closed Loop System
yp2yclosed = [1,0,0];%Picking up only A_z as the plant output.
[Aclosed, Bclosed, Cclosed, Dclosed]=closedloop_from_plantplusctrler(Ap,Bp,Cp,Dp,Ac,Bc1,Bc2,Cc, Dc1, Dc2, yp2yclosed);
sys_closedloop = ss(Aclosed, Bclosed, Cclosed, Dclosed);

%Damping related Design Metrics
[~,zeta_closed] = damp(sys_closedloop);

%Step Response related design metrics
StepRes= stepinfo(sys_closedloop);

%Loop Gain Model at Plant Input:
Lu_ss=series(sys_plant,sys_ctrllr_tocomputeloopgain); 

%Loop Gain Lu : Frequency response
% |w| = frequencies at which frequency response of Lu is evaluated
w_no0=1i*logspace(0.75,3,1000);
angle = 0:0.01:pi/2 ;
w_go_around_origin = abs(w_no0(1))*(cos(angle)+1i*sin(angle));
w = [w_go_around_origin,w_no0] ;
% w=1i*logspace(-1,3,1000);
Lu = freqresp(Lu_ss,w); Lu=squeeze(Lu(1,:,:)); %Loop Gain Model Frequency Response at Plant Input
RD = 1+Lu; %Return Differences at frequencies w
[MinRD,MinRD_index] = min(abs(RD));
[~,PM,~,Wcp] = margin(Lu_ss);
Lu_selected=Lu(real(Lu)<-0.05);
Lu_selected=Lu_selected(real(Lu_selected)>-1);
Lu_real_GM = interp1(imag(Lu_selected),real(Lu_selected),0);
GM = -1/Lu_real_GM;

%Plotting the return Difference Vector
figure('units','normalized','outerposition',[0 0 1 1]);
plot(real(Lu),imag(Lu),'b','LineWidth',1.5); hold on;
plot(real(Lu),-imag(Lu),'b--','LineWidth',1.5);
ax=gca;
plot([-1 real(Lu(MinRD_index))], [0 imag(Lu(MinRD_index))], 'r','LineWidth',1);
hold on; plot(-1,0,'r*','LineWidth',1.2);
xlabel('Real Axis'); ylabel('Imaginary Axis');
title('LQR Design-Nyquist Plot of L_u'); grid('on');
ax2.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
xlim([-25 15]); ylim([-30 30]); axis square;
circle(-1,0,MinRD);
leg=legend('$L_{u_{LQR}}$ for positive $\omega$','$L_{u_{LQR}}$ for negative $\omega$','Return-Difference-Vector','Critical Point','Return-Difference-Disk'); set(leg,'Interpreter','latex');

%% Time Response Simulations for the design point
%TimeSpan:
tspan=[0,11];

% Define pointers to state variables
%--------------------------------------------------------------------------
% Pointers to states
%ProNav State Pointers:
sel_beta = 1; 
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;
%Controller State(eI_A_Z) Pointer :
sel_Xc   = 10; 
%Plant State Pointers:
sel_alpha= 11;
sel_q    = 12; 
sel_dele = 13; 
sel_dele_dot=14;

% ProNav Initial Conditions
    beta_rad = 0;
    RT1      = 20000;%[ft]
    RT2      = 40000;%[ft]
    RM1      = 0;%[ft]
    RM2      = 40000;%[ft]
    VM       = 2813.32; %[ft/s]
    VT       = 300;%[ft/s]
    VT1      = -VT*cos(beta_rad);%[ft/s]
    VT2      =  VT*sin(beta_rad);%[ft/s]
    
% relative positions and velocities
    RTM1 = RT1 - RM1;
    RTM2 = RT2 - RM2;

    % relative distance
    RTM = sqrt(RTM1^2 + RTM2^2);

    % line of sight angle and time derivative
    lambda = atan2( RTM2, RTM1 );

    % missile lead angle
    L = asin( VT*sin( beta_rad + lambda )/VM );

    % missile velocity components
    VM1  = VM*cos(lambda + L + HE_rad);
    VM2  = VM*sin(lambda + L + HE_rad);
    
% Initial ProNav State vector
X_PN0 = [beta_rad, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]';

%Initial Controller State
Xc0 = 0;

%Initial Plant State
Xp0=[0;0;0;0];
    
%%Augmented Initial State Vector
X0 =[X_PN0;Xc0;Xp0];
options = odeset('AbsTol',1e-6,'RelTol',1e-6, 'Events',@event_small_miss_distance);
[time,X_Sol,time_event,X_Sol_event,ie]=ode45(@ode_augmented_pronav_plant_ctrllr,tspan,X0,options);

%Extracting ProNav States from the ode45 solution matrix:
  % target and missile velocity magnitudes
    VT = sqrt( X_Sol(:,sel_VT1).^2 + X_Sol(:,sel_VT2).^2 );

    % relative positions and velocities
    RTM1 = X_Sol(:,sel_RT1) - X_Sol(:,sel_RM1);
    RTM2 = X_Sol(:,sel_RT2) - X_Sol(:,sel_RM2);
    VTM1 = X_Sol(:,sel_VT1) - X_Sol(:,sel_VM1);
    VTM2 = X_Sol(:,sel_VT2) - X_Sol(:,sel_VM2);

    % relative distance
    RTM = sqrt(RTM1.^2 + RTM2.^2);
    MissDistance = RTM(end);

    %Time to collision:
    time2collision = time_event;
    
    % line of sight angle and time derivative
    lambda     = atan2( RTM2, RTM1 );
    lambda_dot = (RTM1.*VTM2 - RTM2.*VTM1)./RTM.^2;

    % closing velocity
    VC = -(RTM1.*VTM1 + RTM2.*VTM2)./RTM;

    % Compute acceleration commands
    A_z_cmd = Np*VC.*lambda_dot;

%Extracting Controller and Plant States from the ode45 solution matrix: 
Xc = X_Sol(:,sel_Xc);
alpha = X_Sol(:,sel_alpha);
q = X_Sol(:,sel_q);
dele = X_Sol(:,sel_dele);
dele_dot = X_Sol(:,sel_dele_dot);
A_z_actual = Z_alpha*alpha + Z_dele*dele; 
MaxElevatorDisplacement_degrees = (180/pi)*max(abs(dele));
MaxElevatorRate_degrees_per_second = (180/pi)*max(abs(dele_dot));

%% Plots
% Engagement: 
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Full-State Feedback - Chosen Design');
ax=gca;
plot(X_Sol(:,sel_RM1), X_Sol(:,sel_RM2),'b','LineWidth',1.5);hold on 
plot(X_Sol(:,sel_RT1), X_Sol(:,sel_RT2),'r','LineWidth',1.5); 
xlabel('Downrange [ft]'); str1=ylabel('Crossrange [ft]');
set(str1,'Interpreter','latex');
str2=title('Missile and Drone Paths'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; axis square;
leg=legend('Missile Path','Drone Path'); set(leg,'Interpreter','latex');

% A_z command tracking during engagement:
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Full-State Feedback - Chosen Design');
ax1=gca;
hold off;
plot(ax1,time,A_z_actual./32.17,'b','LineWidth',1.2); hold on
plot(ax1,time,A_z_cmd./32.17,'r--','LineWidth',1);
xlabel(ax1,'Time'); str1=ylabel(ax1,'$A_{z_{cmd}}$,$A_{z_{actual}}$ (in G''s)');
set(str1,'Interpreter','latex');
str2=title(ax1,'Tracking of $A_{z_{cmd}}$ by Missile''s $A_{z_{actual}}$'); grid(ax1,'on');
set(str2,'Interpreter','latex');
ax1.GridLineStyle = ':'; ax1.GridAlpha = 0.3; ax1.FontSize = 16; ax1.LineWidth = 1.4;
xlim([0,11]); ylim([-10 20]); axis square;
leg=legend('$A_{z_{actual}}$','$A_{z_{cmd}}$'); set(leg,'Interpreter','latex');

%Elevator Deflection
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Full-State Feedback - Chosen Design');
subplot(1,2,1)
ax=gca;
plot(time,180/pi*dele,'b','LineWidth',1.5);hold on;
xlabel('Time'); str1=ylabel('$\delta_e$ in degrees');
set(str1,'Interpreter','latex');
str2=title('Elevator Displacement Vs Time'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; axis square;
leg=legend('$\delta_e$ in degrees'); set(leg,'Interpreter','latex');

%Elevator Rates
subplot(1,2,2)
ax=gca;
plot(time,180/pi*dele_dot,'b','LineWidth',1.5);hold on;
xlabel('Time'); str1=ylabel('$\dot{\delta_e}$ in degrees/sec');
set(str1,'Interpreter','latex');
str2=title('Elevator Rate Vs Time'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; axis square;
leg=legend('$\dot{\delta_e}$ in degrees'); set(leg,'Interpreter','latex');

%Saving Design Parameters from LQR for Observer Design:
qq_LQR = qq_design;
Wcp_LQR = Gain_CrossOverFreqs(winning_design_index);
save LQRdesign qq_LQR Wcp_LQR