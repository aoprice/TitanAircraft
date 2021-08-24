%Sizing of the propeller based off of the Fusion Engine Data
P = 600; % kW
P_real = P *1000;
P1 = P *737.56; %convert to lb*ft/s
P2 = P1 * 0.001818; %horsepower
rho    = TitanAtmosphere(1);
rho2 = atmosisa(3.048); %compare with earth
rho3 = rho2 * 0.062428; %compare with earth, convert to lb/ft^3
rho1 = rho * 0.062428; %convert to lb/ft^3
Vc = 170; %m/s, from LevelFLight
Vc1 = Vc*3.28084; %convert to ft/s
n = 66.6; %rev/s, 4000RPM 
n1 = 251.327; % rad/ s

x =(( rho1 / (P1*n*n) )^( 1 / 5 ));
Cs = Vc1 * x;
Cs_1 = Vc * ((rho / (P_real* n1* n1))^(1/5));
J = 0.7; %from NACA 640 Charts
eff = 0.82; %from NACA 640 charts
D = Vc1/ (n*J); %ft
D_metric = Vc/(n*J); % meters
D1 = 14*(P2^(1/4));

V_tip_static = n *(D/2)^2;
V_tip_helical = (V_tip_static^2 + Vc1^2)^(1/2);

a = 194; %m/s
a1 = a *3.28084; %convert to ft/s

Mach_Limit = V_tip_helical / a1;


%Cp = P_real / (rho * (n^3)*(D^5))
%Ct = Cp * eff/ (J)

%other way of calculating power coefficient 
Cp = P_real/ (rho* n1^3 *D_metric^5);
Cp_alternate = P_real / (.5*rho*(Vc^3)*pi*((3.6465/2)^2))
Ct_alternate = eff * Cp_alternate
Ct = eff * Cp_alternate / J

%thrust needed 
T = 2 * pi*(D/2)^2 *Ct;

