% My ProNav MainScript
% State is y = [beta, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]

% beta - target flight path angle
% RT1  - horizontal position of target wrt inertial cs
% RT2  - vertical position of target wrt inertial cs
% RM1  - horizontal position of missile wrt inertial cs
% RM2  - vertical position of missile wrt inertial cs
% VT1  - horizontal velocity of target wrt inertial cs
% VT2  - vertical velocity of target wrt inertial cs
% VM1  - horizontal velocity of missile wrt inertial cs
% VM2  - vertical velocity of missile wrt inertial cs
clear;
close all;
global HE_rad Np nT
% Pointers to states [beta, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]
sel_beta = 1;
sel_RT1  = 2;
sel_RT2  = 3;
sel_RM1  = 4;
sel_RM2  = 5;
sel_VT1  = 6;
sel_VT2  = 7;
sel_VM1  = 8;
sel_VM2  = 9;
    
%Integration Time Span
tspan = [0,11];

% Parameters
nT = 3*32.2;%Evasive 3g acceleration normal to it's velocity
HE_rad = -20*pi/180;
Np = 3.4; %Typically ranges between 3 to 5

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
    
% Initial condition vector
    y0 = [beta_rad, RT1, RT2, RM1, RM2, VT1, VT2, VM1, VM2]';
    options = odeset('abstol', 1e-6, 'reltol', 1e-6, 'Events',@Pronav_event_small_miss_distance);

% Integrate nonlinear 2-D engagement situation
    [t,y,te,ye,ie] = ode45(@nlinpronav, tspan, y0, options);

    % target and missile velocity magnitudes
    VT = sqrt( y(:,sel_VT1).^2 + y(:,sel_VT2).^2 );

    % relative positions and velocities
    RTM1 = y(:,sel_RT1) - y(:,sel_RM1);
    RTM2 = y(:,sel_RT2) - y(:,sel_RM2);
    VTM1 = y(:,sel_VT1) - y(:,sel_VM1);
    VTM2 = y(:,sel_VT2) - y(:,sel_VM2);

    % relative distance
    RTM = sqrt(RTM1.^2 + RTM2.^2);

    % line of sight angle and time derivative
    lambda     = atan2( RTM2, RTM1 );
    lambda_dot = (RTM1.*VTM2 - RTM2.*VTM1)./RTM.^2;

    % closing velocity
    VC = -(RTM1.*VTM1 + RTM2.*VTM2)./RTM;

    % Compute acc commands
    nc = Np*VC.*lambda_dot;
    
    MissDistance = min(RTM);
 %Plots : 
figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('ProNav Only (No controller and observer involved)');
ax=gca;
plot(y(:,sel_RM1), y(:,sel_RM2),'b','LineWidth',1.5);hold on 
plot(y(:,sel_RT1), y(:,sel_RT2),'r','LineWidth',1.5); 
xlabel('Downrange [ft]'); str1=ylabel('Crossrange [ft]');
set(str1,'Interpreter','latex');
str2=title('Missile and Drone Paths'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; axis square;
leg=legend('Missile Path','Drone Path'); set(leg,'Interpreter','latex');

figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('ProNav Only (No controller and observer involved)');
ax=gca;
plot(t, nc./32.2,'b','LineWidth',1.5);hold on 
xlabel('Time [s]'); str1=ylabel('Acceleration of the Missile(in G)');
set(str1,'Interpreter','latex');
str2=title('Drone pulls a 3 G evasive Maneuver'); grid('on');
set(str2,'Interpreter','latex');
ax.GridLineStyle = ':'; ax.GridAlpha = 0.3; ax.FontSize = 16; ax.LineWidth = 1.4;
hold on; axis square;

