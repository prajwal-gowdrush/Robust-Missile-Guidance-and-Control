%% Mainscript Project: Observer Design
clc; clear; close all;
global Ap Bp Cp_OF A_comp B1_comp B2_comp C_comp D1_comp D2_comp HE_rad Np nT
V0 = 2813.32; %(Trimmed Airspeed in ft/s)(= 2.5 Mach)

%%FROM LQR DESIGN:
load LQRdesign
Q_LQR = [qq_LQR, 0, 0; 0, 0, 0; 0, 0, 0]; %LQR Q Penalty Matrix
R_LQR = 1; %LQR R Penalty Matrix

% Parameters (ProNav)
nT = 3*32.2;%Evasive 3g acceleration normal to it's velocity
HE_rad = -20*pi/180;
Np = 3.4; %Typically ranges between 3 to 5

%Linearized Dynamics of the Missile(Pursuer) (MRAAM):
A = [-1.57, 0 , 0, 1, 0; 0, -0.5, 0.17, 0, -1; -21.13, -2876.7, -2.10, -0.14, -0.05; -82.92, -11.22, -0.01, -0.57, 0; -0.19, -11.86, -0.01, 0, -0.57];
B = [0, -0.1, 0; -0.07, 0, 0.11; -1234.7, -30.49, -1803.2; -4.82, -119.65, -7; 14.84, 0.27, -150.58];

%% Extracted Short-Period Dynamics:
% State = [alpha,q]; , Control=delta_e;
A_sp = A([1,4],[1,4]);
B_sp = B([1,4],2);
Z_alpha = A_sp(1,1)* V0;% Aerodynamic coefficient
Z_dele =  B_sp(1,1)*V0;% Aerodynamic coefficient
C_sp_OF = [Z_alpha 0 ; 0 1]; %Note that alpha is not available as SP output
D_sp_OF = [Z_dele; 0]; %Note that alpha is not available as SP output
% OF stands for Output-Feedback
sys_sp = ss(A_sp,B_sp,C_sp_OF,D_sp_OF);

%% Actuator Dynamics:
w_n = 35*2*pi; %Natural Frequency (in radians per second)
z_damp = 0.71; %Damping Factor
A_act = [0, 1; -w_n^2, -2*z_damp*w_n];
B_act = [0; w_n^2];
C_act = [1,0];
D_act = 0;
sys_act = ss(A_act,B_act,C_act,D_act);

%% Control design model:
C_reg_ctrl = [Z_alpha 0]; 
D_reg_ctrl = Z_dele; 
A_tilda_ctrl = [0 C_reg_ctrl; zeros(size(A_sp,1),1) A_sp];
B_tilda_ctrl = [D_reg_ctrl; B_sp];
B_command_ctrl = [-1; zeros(size(B_sp))];
%LQR Gain Matrix:
[K_lqr,~,~] = lqr(A_tilda_ctrl,B_tilda_ctrl,Q_LQR,R_LQR,0);

%% Observer Design Model:
C_reg_obs = [Z_alpha 0]; 
D_reg_obs = Z_dele; 
A_tilda_obs = [0 C_reg_obs; zeros(size(A_sp,1),1) A_sp];
B_tilda_obs = [D_reg_obs; B_sp];
B_command_obs = [-1; zeros(size(B_sp))];
C_tilda_obs = [1,0,0; 0,0,1];

%% Plant Dynamics:(incorporates the actuator model)
sys_plant = series(sys_act,sys_sp);
[Ap,Bp,Cp_OF,Dp_OF]=ssdata(sys_plant);

%% Index of Design Metrics for Observer-Design:
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
%12=abs(Wcp_Comp - Wcp_LQR)

%% LQ Penalty choices for Observer with corresponding time and frequency responses.
Qo = diag([1,1,1]);
% Qo = diag([1e2,1e2,1e4]);
%Qe is used for determining the gain matrix L:
rho = logspace(0,-2,30); %LTR Parameter
Re = diag([1,1]);
rholen = length(rho);

%%INITIALIZING
does_design_meet_specs = false(rholen,11); % Initializing a logical array
DesignMetrics = zeros(rholen,12);%Initializing Design Metrics Matrix
percentage_satisfaction_specs = zeros(rholen,11);%Initializing Percentage Satisfaction Matrix
Gain_CrossOverFreqs=zeros(rholen,1);
closed_loop_poles = zeros(rholen,8);
% Will be used later to recognize which designs meet the specifications
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Output Feedback - Observer Design - Running through test values of rho');
for h=1:rholen
 Qe = Qo + 1/rho(h)*(B_tilda_obs)*(B_tilda_obs');
[L_dash,~,~] = lqr(A_tilda_obs',C_tilda_obs',Qe,Re,0);
 L = L_dash';
 
%%COMPENSATOR DYNAMICS:
% Compensator Inputs: y=[A_z;q] and r=Az_cmd
%Compensator Output: dele_command
%Compensator States: [eI_A_z; eI_A_z_hat; alpha_hat; q_hat]
%This is obtained by putting together 
%1)an integrator block
%2)Observer dynamics and
%3) Controller Dynamics
%Integrator Matrices:
Ai=0; Bi1=[1,0]; Bi2= -1;
Ci=[1;0]; Di1=[0,0;0,1]; Di2=0;
%Observer Matrix:
A_o_cl = (A_tilda_obs - B_tilda_obs*K_lqr - L*C_tilda_obs);
%Matrices for Compensator Dynamics:
A_comp  = [Ai,zeros(1,size(A_o_cl,2)) ; L*Ci,A_o_cl];
B1_comp = [Bi1; L*Di1];
B2_comp = [Bi2; B_command_obs];
C_comp  = [0,-K_lqr];
D1_comp = [0,0];
D2_comp = 0;
sys_comp = ss(A_comp,B1_comp,C_comp,D1_comp);
sys_comp_tocomputeloopgain = ss(A_comp,B1_comp,-C_comp,-D1_comp);

%%Closed Loop System
yp2yclosed = [1,0];%Picking up only A_z as the plant output.
[Aclosed, Bclosed, Cclosed, Dclosed]=closedloop_from_plantplusctrler(Ap,Bp,Cp_OF,Dp_OF,A_comp,B1_comp,B2_comp,C_comp,D1_comp,D2_comp,yp2yclosed);
sys_closedloop = ss(Aclosed, Bclosed, Cclosed, Dclosed);
closed_loop_poles(h,:) = eig(Aclosed)';

%Step Response related design metrics
StepRes= stepinfo(sys_closedloop);

%Loop Gain Model at Plant Input:
Lu_ss=series(sys_plant,sys_comp_tocomputeloopgain); 

%Loop Gain Model at Plant Output:
Ly_ss=series(sys_comp_tocomputeloopgain,sys_plant); 

%Loop Gain Lu : Frequency response
% |w| = frequencies at which frequency response of Lu is evaluated
w=1i*logspace(-1,3,1000);
Lu = freqresp(Lu_ss,w); 
Lu=squeeze(Lu(1,:,:)); %Loop Gain Model Frequency Response at Plant Input
Ly = freqresp(Ly_ss,w);
Ly_Az_Az = squeeze(Ly(1,1,:));%A_z channel Ly frequency response
Ly_q_q = squeeze(Ly(2,2,:)); % q channel Ly frequency response
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

%Initial Compensator State
X_comp0 = [0;0.01;0.01;0.01];

%Initial Plant State
Xp0=[0;0;0;0];
    
%%Augmented Initial State Vector
X0 =[X_PN0;X_comp0;Xp0];
options = odeset('AbsTol',1e-6,'RelTol',1e-6, 'Events',@obs_event_small_miss_distance);
[time,X_Sol,time_event,X_Sol_event,ie]=ode45(@ode_augmented_pronav_plant_ctrllr_observer,tspan,X0,options);

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

%Extracting Compensator and Plant States from the ode45 solution matrix: 
%Compensator States:
eIAz = X_Sol(:,sel_eIAz);
eIAz_hat = X_Sol(:,sel_eIAz_hat);
alpha_hat = X_Sol(:,sel_alpha_hat);
q_hat = X_Sol(:,sel_q_hat);
%Plant States:
alpha = X_Sol(:,sel_alpha);
q = X_Sol(:,sel_q);
dele = X_Sol(:,sel_dele);
dele_dot = X_Sol(:,sel_dele_dot);
%Calculating A_z_actual, Max elevator displacement and rate:
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
xlabel(ax1,'Time (seconds)'); str1=ylabel(ax1,'$A_{z_{cmd}}$,$A_{z_{actual}}$ (in G''s)');
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
title(ax2,'A portion of the Nyquist Plot of L_{u_{Comp}}'); grid(ax2,'on');
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
DesignMetrics(h,11)=abs(Wcp-Wcp_LQR)/(2*pi);%in Hz
DesignMetrics(h,12)=MissDistance;
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
UB11 = 0.25;

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
does_design_meet_specs(h,11) =DesignMetrics(h,11)<UB11;

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
percentage_satisfaction_specs(h,11)=100*(-DesignMetrics(h,11)+UB11)/abs(UB11);
end

min_satisfaction_percent_2qualify = 1;
min_percentage_satisfaction_specs = min(percentage_satisfaction_specs')';
rho_qualifying = rho(min_percentage_satisfaction_specs>min_satisfaction_percent_2qualify);
MissDistances_qualifying = DesignMetrics(min_percentage_satisfaction_specs>min_satisfaction_percent_2qualify,12);

%Picking the Qualifying Design with minimum miss distance:
minMissDistance=min(MissDistances_qualifying);
winning_design_index = find(DesignMetrics(:,12)==minMissDistance); 

%So, the chosen LTR parameter among the ones tried out is:
rho_design=rho(winning_design_index);

%% Design Charts :  
% The code has already picked out the best penalty matrix for us in terms of minimizing the rise time while meeting the design specifications. 
% Just for the sake of visual confirmation. We plot some design charts
% below:
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(rho,DesignMetrics(:,1),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,1),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[LB1,LB1],'g--','LineWidth',1);
str = xlabel('LTR Parameter $\rho$');set(str,'Interpreter','latex'); 
ylabel('GM at Plant Input [dB]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
leg=legend('GM(dB)','Design Point','Design above line'); set(leg,'Interpreter','latex');

subplot(2,2,2);
plot(rho,DesignMetrics(:,2),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,2),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[LB2,LB2],'g--','LineWidth',1);
str = xlabel('LTR Parameter $\rho$');set(str,'Interpreter','latex'); 
 ylabel('PM at Plant Input [degrees]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([30,55]);
leg=legend('PM(degrees)','Design Point','Design above line'); set(leg,'Interpreter','latex');

subplot(2,2,3);
plot(rho,DesignMetrics(:,3),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,3),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB3,UB3],'k--','LineWidth',1);
str = xlabel('LTR Parameter $\rho$');set(str,'Interpreter','latex'); 
 str=ylabel('Max $\delta_e$(degrees)'); grid on;
set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
leg=legend('$\delta_e$(degrees)','Design Point','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,4);
plot(rho,DesignMetrics(:,4),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,4),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB4,UB4],'k--','LineWidth',1);
str = xlabel('LTR Parameter $\rho$');set(str,'Interpreter','latex'); 
str=ylabel('Max $\dot{\delta_e}$(in degrees/sec)'); grid on;
set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([150,400]);
leg=legend('$\dot{\delta_e}$(degrees\sec)','Design Point','Design below line'); set(leg,'Interpreter','latex');

%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(rho,DesignMetrics(:,5),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,5),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB5,UB5],'k--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex');
str=ylabel('$W_{cp}$ [rad/s] '); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
leg=legend('$W_{cp}$ [rad/s]','Design Point','$\frac{1}{3}\omega_n$ [rad/s] '); set(leg,'Interpreter','latex');

subplot(2,2,2);
plot(rho,DesignMetrics(:,6),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,6),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB6,UB6],'k--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
ylabel('Percentage Undershoot'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([0,12]);
leg=legend('% Undershoot','Design Point','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,3);
plot(rho,DesignMetrics(:,7),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,7),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB7,UB7],'k--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
ylabel('Percentage Overshoot'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([0,12]);
leg=legend('% Overshoot','Design Point','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,4);
plot(rho,DesignMetrics(:,8),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,8),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB8,UB8],'k--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
str=ylabel('Max Sensitivity [dB]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([0,8]);
legend('Max S in A_z and q channels [dB]','Design Point','Design below line');

%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(rho,DesignMetrics(:,9),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,9),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[UB9,UB9],'k--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
str=ylabel('Max Co-Sensitivity [dB]'); grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([0,8]);
legend('Max T in A_z and q channels [dB]','Design Point','Design below line');

subplot(2,2,2);
plot(rho,DesignMetrics(:,10),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,10),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[LB10,LB10],'g--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
str=ylabel('Minimum Re $(\lambda_{cl})$'); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
leg=legend('Minimum Re($\lambda_{cl}$','Design Point','Design above line'); set(leg,'Interpreter','latex');

subplot(2,2,3);
plot(rho,Gain_CrossOverFreqs/(2*pi),'b','LineWidth',1.5);hold on;
plot(rho,Wcp_LQR/(2*pi)*ones(size(rho)),'m--','LineWidth',1.5);hold on;
plot(rho(winning_design_index),Gain_CrossOverFreqs(winning_design_index)/(2*pi),'ro','LineWidth',1.5);
plot([rho(1) rho(end)],[Wcp_LQR/(2*pi)-0.25,Wcp_LQR/(2*pi)-0.25],'g--','LineWidth',1);
plot([rho(1) rho(end)],[Wcp_LQR/(2*pi)+0.25,Wcp_LQR/(2*pi)+0.25],'k--','LineWidth',1);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
str=ylabel('$W_{cp}$[Hz]'); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
ylim([Wcp_LQR/(2*pi)-0.4,Wcp_LQR/(2*pi)+0.4]);
leg=legend('$W_{cp_{Comp}}$','$W_{cp_{LQR}}$','Design Point','Design above line','Design below line'); set(leg,'Interpreter','latex');

subplot(2,2,4);
plot(rho,DesignMetrics(:,11),'b','LineWidth',1.5);hold on;
plot(rho(winning_design_index),DesignMetrics(winning_design_index,11),'ro','LineWidth',1.5);
str1 = xlabel('LTR Parameter $\rho$');set(str1,'Interpreter','latex'); 
str=ylabel('Miss Distance [ft]'); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
xlim([min(rho),max(rho)]);
leg=legend('Miss Distance [ft]','Design Point'); set(leg,'Interpreter','latex');

%% Recomputing time and frequency responses for the design observer penalty matrix

 Qe = Qo + 1/rho_design*(B_tilda_obs)*(B_tilda_obs');
 [L_dash,~,~] = lqr(A_tilda_obs',C_tilda_obs',Qe,Re,0);
 L = L_dash';
 
%%COMPENSATOR DYNAMICS:
% Compensator Inputs: y=[A_z;q] and r=Az_cmd
%Compensator Output: dele_command
%Compensator States: [eI_A_z; eI_A_z_hat; alpha_hat; q_hat]
%This is obtained by putting together 
%1)an integrator block
%2)Observer dynamics and
%3) Controller Dynamics
%Integrator Matrices:
Ai=0; Bi1=[1,0]; Bi2= -1;
Ci=[1;0]; Di1=[0,0;0,1]; Di2=0;
%Observer Matrix:
A_o_cl = (A_tilda_obs - B_tilda_obs*K_lqr - L*C_tilda_obs);
%Matrices for Compensator Dynamics:
A_comp  = [Ai,zeros(1,size(A_o_cl,2)) ; L*Ci,A_o_cl];
B1_comp = [Bi1; L*Di1];
B2_comp = [Bi2; B_command_obs];
C_comp  = [0,-K_lqr];
D1_comp = [0,0];
D2_comp = 0;
sys_comp = ss(A_comp,B1_comp,C_comp,D1_comp);
sys_comp_tocomputeloopgain = ss(A_comp,B1_comp,-C_comp,-D1_comp);

%%Closed Loop System
yp2yclosed = [1,0];%Picking up only A_z as the plant output.
[Aclosed, Bclosed, Cclosed, Dclosed]=closedloop_from_plantplusctrler(Ap,Bp,Cp_OF,Dp_OF,A_comp,B1_comp,B2_comp,C_comp,D1_comp,D2_comp,yp2yclosed);
sys_closedloop = ss(Aclosed, Bclosed, Cclosed, Dclosed);
closed_loop_poles(h,:) = eig(Aclosed)';

%Step Response related design metrics
StepRes= stepinfo(sys_closedloop);

%Loop Gain Model at Plant Input:
Lu_ss=series(sys_plant,sys_comp_tocomputeloopgain); 

%Loop Gain Model at Plant Output:
Ly_ss=series(sys_comp_tocomputeloopgain,sys_plant); 

%Loop Gain Lu : Frequency response
% |w| = frequencies at which frequency response of Lu is evaluated
w=1i*logspace(-1,3,1000);
Lu = freqresp(Lu_ss,w); 
Lu=squeeze(Lu(1,:,:)); %Loop Gain Model Frequency Response at Plant Input
Ly = freqresp(Ly_ss,w);
Ly_Az_Az = squeeze(Ly(1,1,:));%A_z channel Ly frequency response
Ly_q_q = squeeze(Ly(2,2,:)); % q channel Ly frequency response
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

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
ax=gca;
plot(imag(w),20*log10(abs(S_Az_Az)),'b',imag(w),20*log10(abs(S_q_q)),'b--','LineWidth',1.5);hold on;
plot(imag(w),6*ones(size(w)),'k--','LineWidth',1.5);
xlabel('Frequency(in rad/s)'); str1=ylabel('Sensitivity [dB]');
set(str1,'Interpreter','latex');
str2=title('Sensitivity Vs Frequency for the final design'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; 
leg=legend('$A_z-A_z$','$q-q$','Design Upper Cap'); set(leg,'Interpreter','latex');

subplot(1,2,2);
ax=gca;
plot(imag(w),20*log10(abs(T_Az_Az)),'r',imag(w),20*log10(abs(T_q_q)),'r--','LineWidth',1.5);hold on;
plot(imag(w),6*ones(size(w)),'k--','LineWidth',1.5);
xlabel('Frequency(in rad/s)'); str1=ylabel('Co-Sensitivity [dB]');
set(str1,'Interpreter','latex');
str2=title('Co-Sensitivity Vs Frequency for the final design'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; 
leg=legend('$A_z-A_z$','$q-q$','Design Upper Cap'); set(leg,'Interpreter','latex');

%Plot of poles of the system:
ctrllr_poles = eig(A_tilda_ctrl - B_tilda_ctrl*K_lqr);
obsrvr_poles = eig(A_tilda_obs - L*C_tilda_obs );
open_loop_poles = eig(A_tilda_ctrl);
figure('units','normalized','outerposition',[0 0 1 1]);
plot(real(open_loop_poles),imag(open_loop_poles),'ro','MarkerSize',12,'LineWidth',1.5);hold on;
plot(real(ctrllr_poles),imag(ctrllr_poles),'bx','MarkerSize',12,'LineWidth',1.5);
plot(real(obsrvr_poles),imag(obsrvr_poles),'g^','MarkerSize',12,'LineWidth',1.5);
str1 = xlabel('Real Axis');set(str1,'Interpreter','latex'); 
str=ylabel('Imaginary Axis'); grid on;set(str,'Interpreter','latex');
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
ylim([-25,25]);
leg=legend('Open-Loop poles','Controller Poles','Observer Poles'); set(leg,'Interpreter','latex');


%% Nyquist Plot for the Final Design: 
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
title('Observer Design -Nyquist Plot of L_{u_{Comp}}'); grid('on');
ax2.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
xlim([-25 15]); ylim([-30 30]); axis square;
circle(-1,0,MinRD);
leg=legend('$L_{u_{Comp}}$ for positive $\omega$','$L_{u_{Comp}}$ for negative $\omega$','Return-Difference-Vector','Critical Point','Return-Difference-Disk'); set(leg,'Interpreter','latex');

%% Design Time Response:
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
    
% X0 defined before
    
%%Augmented Initial State Vector
X0 =[X_PN0;X_comp0;Xp0];
options = odeset('AbsTol',1e-6,'RelTol',1e-6, 'Events',@obs_event_small_miss_distance);
[time,X_Sol,time_event,X_Sol_event,ie]=ode45(@ode_augmented_pronav_plant_ctrllr_observer,tspan,X0,options);

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

%Extracting Compensator and Plant States from the ode45 solution matrix: 
%Compensator States:
eIAz = X_Sol(:,sel_eIAz);
eIAz_hat = X_Sol(:,sel_eIAz_hat);
alpha_hat = X_Sol(:,sel_alpha_hat);
q_hat = X_Sol(:,sel_q_hat);
%Plant States:
alpha = X_Sol(:,sel_alpha);
q = X_Sol(:,sel_q);
dele = X_Sol(:,sel_dele);
dele_dot = X_Sol(:,sel_dele_dot);
%Calculating A_z_actual, Max elevator displacement and rate:
A_z_actual = Z_alpha*alpha + Z_dele*dele; 
MaxElevatorDisplacement_degrees = (180/pi)*max(abs(dele));
MaxElevatorRate_degrees_per_second = (180/pi)*max(abs(dele_dot));

%% Plots
% Engagement: 
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Final Design (with compensator)');
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
suptitle('Final Design (with compensator)');
ax1=gca;
hold off;
plot(ax1,time,A_z_actual./32.17,'b','LineWidth',1.2); hold on
plot(ax1,time,A_z_cmd./32.17,'r--','LineWidth',1);
xlabel(ax1,'Time (seconds)'); str1=ylabel(ax1,'$A_{z_{cmd}}$,$A_{z_{actual}}$ (in G''s)');
set(str1,'Interpreter','latex');
str2=title(ax1,'Tracking of $A_{z_{cmd}}$ by Missile''s $A_{z_{actual}}$'); grid(ax1,'on');
set(str2,'Interpreter','latex');
ax1.GridLineStyle = ':'; ax1.GridAlpha = 0.3; ax1.FontSize = 16; ax1.LineWidth = 1.4;
xlim([0,11]); ylim([-10 20]); axis square;

%Elevator Deflection
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('Final Design (with compensator)');
subplot(1,2,1)
ax=gca;
plot(time,180/pi*dele,'b','LineWidth',1.5);hold on;
xlabel('Time (seconds)'); str1=ylabel('$\delta_e$ in degrees');
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
xlabel('Time (seconds)'); str1=ylabel('$\dot{\delta_e}$ in degrees/sec');
set(str1,'Interpreter','latex');
str2=title('Elevator Rate Vs Time'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; axis square;
leg=legend('$\dot{\delta_e}$ in degrees'); set(leg,'Interpreter','latex');

%Observer Estimates Tracking actual states:
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(3,1,1)
ax=gca;
plot(time,eIAz,'b',time,eIAz_hat,'r--','LineWidth',1.5);hold on;
xlabel('Time (seconds)'); str1=ylabel('$e_{I_{A_{z}}},\widehat{e_{I_{A_{z}}}}$ in ft/sec');
set(str1,'Interpreter','latex');
str2=title('Integral error in vertical acceleration Vs Time'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; 
leg=legend('$e_{I_{A_{z}}}$','$\widehat{e_{I_{A_{z}}}}$'); set(leg,'Interpreter','latex');

subplot(3,1,2)
ax=gca;
plot(time,alpha,'b',time,alpha_hat,'r--','LineWidth',1.5);hold on;
xlabel('Time (seconds)'); str1=ylabel('$\alpha,\widehat{\alpha}$ in radians');
set(str1,'Interpreter','latex');
str2=title('Angle of Attack  Vs Time'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on;
leg=legend('$\alpha$','$\widehat{\alpha}$'); set(leg,'Interpreter','latex');

subplot(3,1,3)
ax=gca;
plot(time,q,'b',time,q_hat,'r--','LineWidth',1.5);hold on;
xlabel('Time (seconds)'); str1=ylabel('$q,\widehat{q}$ in rad/sec');
set(str1,'Interpreter','latex');
str2=title('Pitch Rate Vs Time'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; 
leg=legend('$q$','$\widehat{q}$'); set(leg,'Interpreter','latex');
