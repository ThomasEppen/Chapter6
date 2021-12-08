%================================================
%     Matlab Script File used to linearize the 
%     non-linear F-16 model. The program will 
%     Extract the longitudal and lateral 
%     direction matrices.  These system matrices 
%     will be used to create pole-zero mapping
%     and the bode plots of each to each control
%     input.
% Author: Richard S. Russell
% 
% Edit: Ewoud Smeur (2021)
%================================================
clear;

global fi_flag_Simulink

newline = sprintf('\n');

%% Trim aircraft to desired altitude and velocity
%%
altitude = input('Enter the altitude for the simulation (ft)  :  ');
velocity = input('Enter the velocity for the simulation (ft/s):  ');

%% Initial guess for trim
%%
thrust = 5000;          % thrust, lbs
elevator = -0.09;       % elevator, degrees
alpha = 8.49;              % AOA, degrees
rudder = -0.01;             % rudder angle, degrees
aileron = 0.01;            % aileron, degrees

%% Find trim for Hifi model at desired altitude and velocity
%%
disp('Trimming High Fidelity Model:');
fi_flag_Simulink = 1;
[trim_state_hi, trim_thrust_hi, trim_control_hi, dLEF, xu_hi] = trim_F16(thrust, elevator, alpha, aileron, rudder, velocity, altitude);

%% Find the state space model for the hifi model at the desired alt and vel.
trim_state_lin = trim_state_hi; trim_thrust_lin = trim_thrust_hi; trim_control_lin = trim_control_hi;
operating_point = operpoint('LIN_F16Block'); % retrieves initial conditions from integrators
operating_point.Inputs(1).u = trim_thrust_lin; operating_point.Inputs(2).u = trim_control_lin(1);
operating_point.Inputs(3).u = trim_control_lin(2); operating_point.Inputs(4).u = trim_control_lin(3);

SS_hi = linearize('LIN_F16Block');

%% Find trim for lofi model at desired altitude and velocity
%%
disp('Trimming Low Fidelity Model:');
fi_flag_Simulink = 0;
[trim_state_lo, trim_thrust_lo, trim_control_lo, dLEF, xu_lo] = trim_F16(thrust, elevator, alpha, aileron, rudder, velocity, altitude);

%% Find the state space model for the lofi model at the desired alt and vel.
%%
trim_state_lin = trim_state_lo; trim_thrust_lin = trim_thrust_lo; trim_control_lin = trim_control_lo;
operating_point = operpoint('LIN_F16Block'); % retrieves initial conditions from integrators
operating_point.Inputs(1).u = trim_thrust_lin; operating_point.Inputs(2).u = trim_control_lin(1);
operating_point.Inputs(3).u = trim_control_lin(2); operating_point.Inputs(4).u = trim_control_lin(3);

SS_lo = linearize('LIN_F16Block');

%%%%%%%%%%%%%%%%%%%%%%%
%% Longitudal Direction
%%%%%%%%%%%%%%%%%%%%%%%

long_states = [3 5 7 8 11 13 14];
long_inputs = [1 2];
long_outputs = [3 5 7 8 11];

SS_long_lo = ss(SS_lo.A(long_states,long_states), SS_lo.B(long_states,long_inputs), SS_lo.C(long_outputs,long_states), SS_lo.D(long_outputs,long_inputs));
SS_long_hi = ss(SS_hi.A(long_states,long_states), SS_hi.B(long_states,long_inputs), SS_hi.C(long_outputs,long_states), SS_hi.D(long_outputs,long_inputs));

SS_long_lo.StateName = SS_lo.StateName(long_states);
SS_long_hi.StateName = SS_hi.StateName(long_states);

SS_long_lo.InputName= SS_lo.InputName(long_inputs);
SS_long_hi.InputName= SS_hi.InputName(long_inputs);

%%%%%%%%%%%%%%%%%%%%
%% Lateral Direction
%%%%%%%%%%%%%%%%%%%%

lat_states = [4 6 7 9 10 12 13 15 16];
lat_inputs = [1 3 4];
lat_outputs = [4 6 7 9 10 12];

SS_lat_lo = ss(SS_lo.A(lat_states,lat_states), SS_lo.B(lat_states,lat_inputs), SS_lo.C(lat_outputs,lat_states), SS_lo.D(lat_outputs,lat_inputs));
SS_lat_hi = ss(SS_hi.A(lat_states,lat_states), SS_hi.B(lat_states,lat_inputs), SS_hi.C(lat_outputs,lat_states), SS_hi.D(lat_outputs,lat_inputs));

SS_lat_lo.StateName = SS_lo.StateName(lat_states);
SS_lat_hi.StateName = SS_hi.StateName(lat_states);

SS_lat_lo.InputName= SS_lo.InputName(lat_inputs);
SS_lat_hi.InputName= SS_hi.InputName(lat_inputs);

newAlongitude = SS_lat_lo.A([3 4 2 5],[3 4 2 5]);
newBlongitude = SS_lat_lo.A([3 4 2 5],7);
newClongitude = SS_lat_lo.C([3 4 2 5],[3 4 2 5]);
newDlongitude = SS_lat_lo.D([3 4 2 5],2);
 
newSSlongitude = ss(newAlongitude, newBlongitude, newClongitude, newDlongitude);
[V_eig_longitude, D_eig_longitude] = eig(newAlongitude);
eig_longitude_shortperiod = D_eig_longitude(1,1);
eig_longitude_phugoid = D_eig_longitude(3,3);
 
%% Reduced state Space Matrixes
%%
SS_longitude_lo = modred(SS_long_lo, [1,6,7],'Truncate');
A_long_lo = SS_longitude_lo.A;
B_long_lo = SS_long_lo.A([2,3,4,5],[6,7]);
C_long_lo = SS_long_lo.C([1 2 3 4 ], [1 2 3 4 ]);
D_long_lo = SS_long_lo.D([1 2 3 4 ], [1 2]);

SS_lateral_lo = modred(SS_lat_lo, [2,3,7,8,9], 'Truncate');
A_lat_lo = SS_lateral_lo.A;
B_lat_lo = SS_lat_lo.A([3,4,5,6],[8,9]);
C_lat_lo = SS_lat_lo.C([1 2 3 4 ], [1 2 3 4 ]);
D_lat_lo = SS_lat_lo.D([1 2 3 4 ], [1 2]);

[Omega_n_long,Z_long,~]=damp(A_long_lo);

T_long1=-log(0.5)/(Z_long(1)*Omega_n_long(1));
T_long2=-log(0.5)/(Z_long(3)*Omega_n_long(3));

Omega_d_long1=Omega_n_long(1)*sqrt(1-Z_long(1)^2);
Omega_d_long2=Omega_n_long(3)*sqrt(1-Z_long(3)^2);

Period_long1=2*pi/Omega_d_long1;
Period_long2=2*pi/Omega_d_long2;

[Omega_n_lat,Z_lat,~]=damp(A_lat_lo);

T_lat1=-log(0.5)/(Z_lat(1)*Omega_n_lat(1));
T_lat2=-log(0.5)/(Z_lat(3)*Omega_n_lat(3));
T_lat3=-log(0.5)/(Z_lat(4)*Omega_n_lat(4));

Omega_d_lat1=Omega_n_lat(1)*sqrt(1-Z_lat(1)^2);
Period_lat1=2*pi/Omega_d_lat1;

tao1=1/(Z_lat(3)*Omega_n_lat(3));
tao2=1/(Z_lat(4)*Omega_n_lat(4));

%% plot eigenmotions

%transferfunctions elevator and rudder
[Numerator_long,Denominator_long]=ss2tf(A_long_lo,B_long_lo,C_long_lo,D_long_lo,2);
[Numerator_lat,Denominator_lat]=ss2tf(A_lat_lo,B_lat_lo,C_lat_lo,D_lat_lo,2);

tfphugoid=tf(Numerator_long(1,:),Denominator_long);

figure(12);
t=0:0.1:1200;
step(-(180/pi)*tfphugoid,t);
axis([0 1200 -30 35]);
xlabel('Time');
ylabel('Theta [deg]');
title('Phugoid');
%# horizontal line
%hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
%changedependvar(hy,'y');
grid on;

tfshortp=tf(Numerator_long(4,:),Denominator_long);

figure(13);
t=0:0.1:10;
step(-tfshortp,t);
axis([0 10 0 5]);
xlabel('Time');
ylabel('q [deg/s]');
title('Short Period');
%horizontal line
%hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
%changedependvar(hy,'y');
grid on;

tfdrollrate=tf(-Numerator_lat(3,:),Denominator_lat);
tfdyawrate=tf(-Numerator_lat(2,:),Denominator_lat);

figure(14);
t=0:0.1:15;
impulse(tfdrollrate,t);
plot((180/pi)*impulse(tfdrollrate,t),(180/pi)*impulse(tfdyawrate,t))

grid on;
title('Dutch Roll');
xlabel('p [deg/s]');
ylabel('r [deg/s]');


tfspiral=tf(Numerator_lat(1,:),Denominator_lat);

figure(15);
t=0:0.1:600;
step(-tfspiral,t);
axis([0 600 -40 360]);
xlabel('Time');
ylabel('phi [deg]');
title('Spiral');
%# horizontal line
%hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
%changedependvar(hy,'y');
grid on;

%transferfunction aileron
[NUMlatail,DENlatail]=ss2tf(A_lat_lo,B_lat_lo,C_lat_lo,D_lat_lo,1);

tfaperiod=tf(NUMlatail(3,:),DENlatail);

figure(16);
t=0:0.01:10;
step(-(180/pi)*tfaperiod,t);
axis([0 10 0 20]);
xlabel('Time');
ylabel('p [deg/s]');
title('Aperiodic Roll');
%# horizontal line
%hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
%changedependvar(hy,'y');
grid on;


sys_long = ss(A_long_lo,B_long_lo,C_long_lo,D_long_lo);
sys_lat = ss(A_lat_lo,B_lat_lo,C_lat_lo,D_lat_lo);


%% Longitudinal
Eig_long=eig(A_long_lo);

I_Eig_long=imag(Eig_long);
RE_Eig_long=real(Eig_long);
c_long=[RE_Eig_long(1,1) RE_Eig_long(3,1)]';
n_long=[I_Eig_long(1,1) I_Eig_long(3,1)]';


T_12_long=-0.693./c_long; % [s]
P_long= 2*pi./n_long; % [s]

Damp_long=-c_long./(sqrt(n_long.^(2)+c_long.^(2)));
wn_long=(sqrt(n_long.^(2)+c_long.^(2))).*sqrt(1-Damp_long.^2); % [rad/s]

T_long = 0:0.01:300;
u_long =[0*T_long' (-1*(pi/180)*ones(1,length(T_long))-1.7797*(pi/180))'];
Y_long = lsim(sys_long,u_long,T_long');

V_long = Y_long(:,1) + 900*0.3048;
Theta_long = Y_long(:,3) + 0;
Alpha_long = Y_long(:,2) + 1.6235*pi/180;
Pitch_long = Y_long(:,4);

figure(1);
subplot(2,1,1);
plot(T_long,(Alpha_long*(180/pi)));
axis([0 5 0 15]);
xlabel('time (s)'); ylabel('angle of attack [deg]');
title('Short period characteristics')
grid on;
subplot(2,1,2);
plot(T_long,Pitch_long);
axis([0 5 0 0.4]);
xlabel('time (s)'); ylabel('pitch rate [rad/s]');
grid on;

figure(2);
subplot(2,1,1);
plot(T_long,V_long);
xlabel('Time [s]'); ylabel('Airspeed [m/s]');
grid on;
title('Phugoid characteristics');
subplot(2,1,2);
plot(T_long,Pitch_long);
xlabel('Time [s]'); ylabel('Pitch rate [rad/s]');
grid on;

%% Lateral
Eig_lat=eig(A_lat_lo);

I_Eig_lat=imag(Eig_lat);
RE_Eig_lat=real(Eig_lat);
c_lat=[RE_Eig_lat(2,1) RE_Eig_lat(3,1) RE_Eig_lat(4,1)]';
n_lat = [I_Eig_lat(1) 0 0]';

T_12_lat=-0.693./c_lat; % [s]
P_lat= 2*pi./n_lat(1); % [s]

Damp_lat=-c_lat./(sqrt(n_lat.^(2)+c_lat.^(2)));
wn_lat=(sqrt(n_lat.^(2)+c_lat.^(2))).*sqrt(1-Damp_lat.^2); % [rad/s]
tau_lat=-1./Eig_lat(3:4);


T_lat1 = 0:0.01:20;
T_lat2 = 0:0.01:80;
u_lat_DR =[0*T_lat1' (-2*(pi/180)*ones(1,length(T_lat1)))'];
u_lat_AR=[(-5*(pi/180)*ones(1,length(T_lat1)))' 0*T_lat1'];
u_lat_SP= zeros(2,length(T_lat2));
x_lat_SP=pi/180*[0 -13 0 0];

Y_lat_DR = lsim(sys_lat,u_lat_DR,T_lat1');
Y_lat_AR = lsim(sys_lat,u_lat_AR,T_lat1');
Y_lat_SP = lsim(sys_lat,u_lat_SP,T_lat2',x_lat_SP);

figure(3);
subplot(2,1,1);
plot(T_lat1,Y_lat_DR(:,3))
ylabel(' Roll rate [rad/s]')
hold on
grid on
title('Dutch Roll characteristics')
subplot(2,1,2);
plot(T_lat1,Y_lat_DR(:,4))
hold on
xlabel('Time [s]')
ylabel(' Yaw rate [rad/s]')
grid on;


figure(4);
subplot(2,1,1);
plot(T_lat1,Y_lat_AR(:,3))
hold on
grid on
ylabel('Roll rate [rad/s]')
title('Aperiodic roll characteristics')
subplot(2,1,2);
plot(T_lat1,Y_lat_AR(:,4))
hold on
xlabel('Time [s]')
ylabel('Yaw rate [rad/s]')
grid on;


figure(5);
subplot(2,1,1);
plot(T_lat2,Y_lat_SP(:,2))
ylabel('Roll angle [deg]')
grid on
hold on
title('Spiral characteristics')
subplot(2,1,2);
plot(T_lat2,Y_lat_SP(:,4))
hold on
xlabel('Time [s]')
ylabel('Yaw rate [rad/s]')
grid on;

%% Display SP

disp('eigenvalue short period')
disp(Eig_long(1:2))

disp('Period short period [s]')
disp(P_long(1))

disp('half amplitude short period [s]')
disp(T_12_long(1))

disp('Damping ratio short period')
disp(Damp_long(1))

disp('natural frequency short period [rad/s]')
disp(wn_long(1))

disp('-----------------------------------------------')
%% Display phugoid
disp('eigenvalue Phugoid')
disp(Eig_long(3:4))

disp('Period short Phugoid [s]')
disp(P_long(2))

disp('half amplitude Phugoid [s]')
disp(T_12_long(2))

disp('Damping ratio Phugoid')
disp(Damp_long(2))

disp('natural frequency Phugoid [rad/s]')
disp(wn_long(2))

disp('-----------------------------------------------')

%% Display Dutch roll
disp('eigenvalue Dutch roll')
disp(Eig_lat(1:2))

disp('Period short Dutch roll [s]')
disp(P_lat(1))

disp('half amplitude Dutch roll [s]')
disp(T_12_lat(1))

disp('Damping ratio Dutch roll')
disp(Damp_lat(1))

disp('natural frequency Dutch roll [rad/s]')
disp(wn_lat(1))

disp('-----------------------------------------------')

%% Display aperiodic roll
disp('eigenvalue aperiodic roll')
disp(Eig_lat(3))

disp('half amplitude aperiodic roll [s]')
disp(T_12_lat(2))

disp('Time constant tau aperiodic roll [s]')
disp(tau_lat(1))


disp('natural frequency aperiodic roll [rad/s]')
disp(wn_lat(2))

disp('-----------------------------------------------')
%% Display spiral
disp('eigenvalue spiral')
disp(Eig_lat(4))

disp('half amplitude spiral [s]')
disp(T_12_lat(3))

disp('Time constant tau spiral [s]')
disp(tau_lat(2))


disp('natural frequency spiral [rad/s]')
disp(wn_lat(2))

disp('-----------------------------------------------')





















%% Longitudinal Poles
%%
figure(2); 
pzmap(SS_long_hi, 'r', SS_long_lo, 'b');
title_string = sprintf('Altitude = %.2f ft Velocity = %.2f ft/s\nLongitudal Directional Poles\n Blue = lofi Red = hifi.', altitude, velocity);
title(title_string);
sgrid;

%% Lateral Poles
%%
figure(3); 
pzmap(SS_lat_hi, 'r', SS_lat_lo, 'b');
title_string = sprintf('Altitude = %.2f ft Velocity = %.2f ft/s\nLateral Directional Poles\n Blue = lofi Red = hifi.', altitude, velocity);
title(title_string);
sgrid;

%% Bode plot longitudinal 

% Choose an input and 
input = 2;
output = 3;

omega = logspace(-2,2,100);

figure
bode(SS_long_hi(output,input),omega)
hold on;
bode(SS_long_lo(output,input),omega)
legend('hifi','lofi')

%% Bode plot lateral 

% Choose an input and 
input = 2;
output = 3;

omega = logspace(-2,2,100);

figure
bode(SS_lat_hi(output,input),omega)
hold on;
bode(SS_lat_lo(output,input),omega)
legend('hifi','lofi')