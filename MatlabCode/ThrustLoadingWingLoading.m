% Thrust Loading vs. Wing Loading graph as part of initial aircraft sizing 
% author: Annie Price 
% aoprice@princeton.edu 

%constants

%get rho
h = zeros(1,1000);
for i= 1:length(h)
rho = TitanAtmosphere(h);
end 

s_to = 500; %m 
s_land = 500; %m
g = 1.352; %m/s^2, acc due to gravity on Titan
eff_prop = 0.87; %propeller efficiency, from NACA chart 
P = 1; %MW, constant electric power
P1 = 0.7; %Mw, used to size propellers 
hp1 = 1341.02; %mechanical horsepower conversion for 1 MW
hp2 = 938.715; %power conversion for 0.8 MW
dry_mass = 2000; %kg, dry mass
engine_mass = 1000; %kg 
TO_weight = 2704; %N
V_cruise = 100; %m/s, cruise velocity
M = 0;85; %cruise mach number 
Clmax = 3; %max lift coefficient estimation
a = 194; %m/s, speed of sound Titans surface
T_to = 94.6254; %N, thrust from engine at takeoff
T_cruise = 142.6507; %N, thrust at 1km from the engine
Cd0 = 0.006; %from RHS2DTitan 
Clo = 0.3; %just an estimation
Cl_to = Clmax/1.2; %from Raymer(?)
AR = 1.7; %from RHS2Dtitan
eps = 0.9; %from RHS2DTitan
k = 1/(pi * AR * eps); 
Clalpha = 6.2832; %from RHS2DTitan
Cdto = Cd0 + k*Cl_to^2; %takeoff drag
rhosl = 1.8798; %from wikepedia, g/cm^3
TOP = 170; %from raymer graph, inclue in latex
q =  1/2*V_cruise^2*rho(1);
q_climb = 1/2*50^2*rhosl;
sigma = rho(1)/ rhosl;
n = 1; %estimation
LiftDragCruise = 18;
Vstall = 30.641; %m/s
Beta = 1; %fuel fraction
a_lapse = T_cruise/T_to; 


xload = linspace(1,50,50);
wingload = (Vstall*Clmax/rho(1))^1/2;
wingloadproptakeoff = (TOP)*sigma*Cl_to*(hp2/TO_weight);
Thrustload_cruise = q./(xload).*(Cd0) +  k.*(xload./q-Clo).^2;
Thrustload_cruise1 = 1/(LiftDragCruise); %Raymer
Thrustload_takeoff = Thrustload_cruise *(T_to / T_cruise);
%Thrustload_cruise = T_cruise / TO_weight; 
%Thrustload_cruise = LiftDragInverse_cruise; %thrustmatching
Thrustload_TO = T_to/ TO_weight;
Thrustload_climb = (q_climb*Cd0)./xload + xload.*(n/q_climb*pi*AR*eps);
Thrustload_takeoff2 = 1/ a_lapse * (2*(((Cd0 + k*Cl_to^2)*k)^(1/2)+k));


%powerloading equations for propeller aicraft
%Powerload_takeoff = (V_cruise/550)*(1/LiftDragCruise)*1*1;
%Thrustload_propeller = (550*0.87/ V_cruise)*(hp2/TO_weight);


%convert above thrust-to-weight equtions into horsepower-to-weight
%equations
wingload_maxprop = q *(pi*AR*eps*Cd0/3)^1/2;
Thrustload_propeller = (550*0.87/ V_cruise)*(hp2/TO_weight); %turn this into a conversion factor
ConvCruise = (V_cruise/(550*eff_prop));
Powerload_cruise = Thrustload_cruise * ConvCruise; 
Powerload_takeoff = Thrustload_takeoff * (ConvCruise)*(T_cruise/T_to); 
Powerload_climb = Thrustload_climb * ConvCruise; 


figure (2) 
x = linspace(1,50,50)
y = linspace(1,2,5)
plot(x, Thrustload_cruise)
hold on
plot(x, Thrustload_takeoff)
hold on 
plot (x, Thrustload_climb) 
hold on 
xline (wingload)   
hold on
yline (Thrustload_TO)
hold off
xlim([0 50])
ylim([0 8])
legend ('cruise', 'takeoff', 'climb', 'stall')



figure(3) 
x = linspace(1,50,50)
y = linspace(1,2,5)
plot (x, Powerload_takeoff)
hold on 
plot (x, Powerload_climb)
hold on 
plot(x, Powerload_cruise)
hold on 
xline(wingload)
hold on 
xline (wingload_maxprop)
hold off
xlim([0 50])
ylim([0 10])
legend ('takeoff', 'climb', 'cruise', 'stall')


%%Tail sizing calculations -- very initial 
C_vt = 0.07; %from raymer, for a jet fighter plane
b_w = 44.487; % feet 
S_w = 1164.22455; %ft 
L_vt = 0.55* 54; %from raymer, 55% times fuselage length (rough_
S_Vt = C_vt * b_w *S_w / L_vt;
