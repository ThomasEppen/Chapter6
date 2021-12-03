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

[Wn_long,Z_long,P_long]=damp(A_long_lo);


T_long1=-log(0.5)/(Z_long(1)*Wn_long(1));
T_long2=-log(0.5)/(Z_long(3)*Wn_long(3));

Wd_long1=Wn_long(1)*sqrt(1-Z_long(1)^2);
Wd_long2=Wn_long(3)*sqrt(1-Z_long(3)^2);

Period_long1=2*pi/Wd_long1;
Period_long2=2*pi/Wd_long2;

[Wn_lat,Z_lat,P_lat]=damp(A_lat_lo);


T_lat1=-log(0.5)/(Z_lat(1)*Wn_lat(1));
T_lat2=-log(0.5)/(Z_lat(3)*Wn_lat(3));
T_lat3=-log(0.5)/(Z_lat(4)*Wn_lat(4));

Wd_lat1=Wn_lat(1)*sqrt(1-Z_lat(1)^2);
Period_lat1=2*pi/Wd_lat1;

tao1=1/(Z_lat(3)*Wn_lat(3));
tao2=1/(Z_lat(4)*Wn_lat(4));


%% plot eigenmotions

%transferfunctions elevator and rudder
[NUMlong,DENlong]=ss2tf(A_long_lo,B_long_lo,C_long_lo,D_long_lo,2);
[NUMlat,DENlat]=ss2tf(A_lat_lo,B_lat_lo,C_lat_lo,D_lat_lo,2);

tfphugoid=tf(NUMlong(1,:),DENlong);

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

tfshortp=tf(NUMlong(4,:),DENlong);

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

tfdrollrate=tf(-NUMlat(3,:),DENlat);
tfdyawrate=tf(-NUMlat(2,:),DENlat);

figure(14);
t=0:0.1:15;
impulse(tfdrollrate,t);
plot((180/pi)*impulse(tfdrollrate,t),(180/pi)*impulse(tfdyawrate,t))

grid on;
title('Dutch Roll');
xlabel('p [deg/s]');
ylabel('r [deg/s]');


tfspiral=tf(NUMlat(1,:),DENlat);

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