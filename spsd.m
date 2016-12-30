function varagout = spsd(system_velocity,chamber_temperature,chamber_pressure,chamber_geometry,exit_radius,area_ratio,propellant,altitude,initial_mass)
% Title: Spacecraft Propulsion System Design
% Purpose: This script takes a series of mostly physical atributes and
% calculates the performace of the propulsion system.
% Author: James Emerson Parkus
% Date: December 26th, 2016
% Notes: See the document that goes over this script. It will explain how
% the chamber geometry, propellant, and altitude factors are determined.

%% Constants
standard_gravity = 9.8066; %m/s^2
universal_gas_constant = 8.3144598; % [J/(mol*K)]
earth_radius = 6.37*10^6;
gravity = standard_gravity*(earth_radius/(earth_radius+altitude))^2;

%% Nozzle Conditions
exit_area = pi*exit_radius^2;

throat_area = exit_area/area_ratio;

throat_radius = sqrt(throat_area/pi);

conical_half_angle = 15;

nozzle_length = (exit_radius-throat_radius)/tan(deg2rad(conical_half_angle));

%% Units
distance_unit = '[m]';
velocity_unit = '[m/s]';
acceleration_unit = '[m/s^2]';
force_unit = '[N]';
energy_unit = '[J]';
time_unit = '[s]';
mass_unit = '[kg]';
pressure_unit = '[Pa]';
unitless = '[-]';
massflow_unit = '[kg/s]';
density_unit = '[kg/m^3]';
temperature_unit = '[K]';
gas_constant_unit = '[J/kgK]';
molar_mass_unit = '[kg/mol]';
area_unit = '[m^2]';
volume_unit = '[m^3]';
newton_time = '[Ns]';

%% Chamber Conditions
propellant_data = xlsread('RPSD Propellant Info.xlsx');
propellant_name_data = propellant_data(:,2);
molar_mass_data = propellant_data(:,3);
heat_ratio_data = propellant_data(:,4);

switch propellant
    case 1
        PropellantName = 'Hydrogen';
    case 2
        PropellantName = 'Helium';
    case 3
        PropellantName = 'Nitrogen';
    case 4
        PropellantName = 'Freon [R-22]';
    case 5
        PropellantName = 'Neon';
    case 6
        PropellantName = 'Xenon';
    case 7
        PropellantName = 'Argon';
end

molar_mass = molar_mass_data(propellant,1);
specific_heat_ratio = heat_ratio_data(propellant,1);

specific_gas_constant = universal_gas_constant/molar_mass;

%% Throat Conditons
chamber_pressure(1) = chamber_pressure;

throat_temperature = chamber_temperature(1)*(2/(specific_heat_ratio+1));

throat_pressure = chamber_pressure(1)*(2/(specific_heat_ratio+1))^(specific_heat_ratio/(specific_heat_ratio-1));

throat_density = throat_pressure/(specific_gas_constant*throat_temperature);

throat_velocity = sqrt(specific_heat_ratio*specific_gas_constant*throat_temperature);

mass_flowrate_I = throat_density*throat_velocity*throat_area;

mass_flowrate_II = throat_area*chamber_pressure(1)*specific_heat_ratio*sqrt((2/(specific_heat_ratio+1))^((specific_heat_ratio+1)/(specific_heat_ratio-1)))/(sqrt(specific_heat_ratio...
    *specific_gas_constant*chamber_temperature));

%% Exit Conditions
exit_mach = solve_mach(exit_area,throat_area,specific_heat_ratio);

exit_temperature(1) = chamber_temperature(1)/(1+(specific_heat_ratio-1)/2*exit_mach^2);

mass_flowrate_III = specific_heat_ratio*chamber_pressure(1)*exit_area*exit_mach/sqrt(specific_heat_ratio*specific_gas_constant*chamber_temperature)*(1+(specific_gas_constant-1)/2*exit_mach^2)^((2-specific_heat_ratio)/2);

temperature_ratio = chamber_temperature/exit_temperature;

exhaust_velocity = exit_mach*sqrt(specific_heat_ratio*specific_gas_constant*exit_temperature);

exit_pressure = chamber_pressure(1)/(temperature_ratio)^(specific_heat_ratio/(specific_heat_ratio-1));

pressure_ratio = exit_pressure/chamber_pressure(1);

%% System Performance

effective_exhaust_velocity = exhaust_velocity+exit_pressure*exit_area/mass_flowrate_I;

thrust_coefficient(1) = sqrt(2*specific_heat_ratio^2/(specific_heat_ratio-1))*(2/(specific_heat_ratio+1))^((specific_heat_ratio+1)/(specific_heat_ratio-1))...
    *(1-(pressure_ratio)^((specific_heat_ratio-1)/(specific_heat_ratio+1)))+pressure_ratio*area_ratio;

thrust_force_I = thrust_coefficient*chamber_pressure(1)*throat_area;

thrust_force_II = effective_exhaust_velocity*mass_flowrate_I;

thrust_force_III = mass_flowrate_I*exhaust_velocity+exit_pressure*exit_area;

specific_impulse_I = thrust_force_I*standard_gravity/mass_flowrate_I;

specific_impulse_II = effective_exhaust_velocity/standard_gravity;

thrust_force_IV(1) = mass_flowrate_I*specific_impulse_II/standard_gravity;

final_mass = initial_mass/exp(system_velocity/effective_exhaust_velocity);

propellant_mass(1) = abs(final_mass - initial_mass);

burn_time = propellant_mass/mass_flowrate_I;

mass_ratio = initial_mass/final_mass;

%% Tank Conditions
switch chamber_geometry
    case 1
        geometry = 'Spherical';
        radius = 0.05;
        spherical_volume = propellant_mass*specific_gas_constant*chamber_temperature/chamber_pressure(1)
        volume = spherical_volume;
    case 2
        geometry = 'Rectangular';
        length = 0.1;
        height = 0.1;
        width = 0.05;
        rectangular_volume = propellant_mass*specific_gas_constant*chamber_temperature/chamber_pressure(1);
        volume = rectangular_volume;
    case 3
        geometry = 'Cylindrical';
        height = 0.1;
        radius = 0.05;
        cylindrical_volume = propellant_mass*specific_gas_constant*chamber_temperature/chamber_pressure(1)
        volume = cylindrical_volume;
    case 4
        geometry = 'Toroidal';
        minor_radius = 1;
        major_radius = 1;
        torus_volume = propellant_mass*specific_gas_constant*chamber_temperature/chamber_pressure(1)
        volume = torus_volume;
end

%% Iteration Sequence
dt = 1*10^-3;
i = 1;
t(1) = 0;
lim = 1*10^3

while chamber_pressure >= lim
    dm = mass_flowrate_I*dt;
    propellant_mass(i+1) = propellant_mass(i) - dm;
    chamber_pressure(i+1) = propellant_mass(i+1)*specific_gas_constant*chamber_temperature/volume;
    exit_pressure(i+1) = chamber_pressure(i+1)/((temperature_ratio)^(specific_heat_ratio/(specific_heat_ratio-1)));
    thrust_coefficient(i+1) = sqrt(2*specific_heat_ratio^2/(specific_heat_ratio-1))*(2/(specific_heat_ratio+1))^((specific_heat_ratio+1)/(specific_heat_ratio-1))...
        *(1-(exit_pressure(i+1)/chamber_pressure(i+1))^((specific_heat_ratio-1)/specific_heat_ratio))+(exit_pressure(i+1)/chamber_pressure(i+1))*(area_ratio);
    thrust_force_I(i+1) = thrust_coefficient(i+1)*chamber_pressure(i+1)*throat_area;
    acceleration(i+1) = thrust_force_I(i+1)*dt*1/dm;
    velocity(i+1) = acceleration(i+1)*dt;
    distance(i+1) = velocity(i+1)*dt;
    
    impulse_bit(i) = thrust_force_I(i+1)*dt;
    velocity_bit(i) = velocity(i+1)*dt;
    
    i = i + 1;
    t(i+1) = t(1) + dt;
end

total_impulse = sum(impulse_bit); 
specific_impulse = total_impulse/(propellant_mass(1)*standard_gravity);

delta_velocity = sum(velocity_bit);

%% Results
linedivider = '------------';
result = {'Propellant','',PropellantName;
    'Specific gas constant',specific_gas_constant, gas_constant_unit;
    'Molar Mass',molar_mass,molar_mass_unit;
    'Specific Heat Ratio',specific_heat_ratio,unitless;
    linedivider,'','';
    'Nozzle Conditions','','';
    linedivider,'','';
    'Nozzle Length',nozzle_length,distance_unit;
    'Exit Radius',exit_radius,distance_unit;
    'Throat Radius',throat_radius,distance_unit;
    'Exit Area',exit_area,area_unit;
    'Throat Area',throat_area,area_unit;
    linedivider,'','';
    'Chamber Conditions','','';
    linedivider,'','';
    'Tank Geometry',geometry,'';
    'Volume',volume,volume_unit;
    'Chamber Pressure',chamber_pressure(1),pressure_unit;
%     'Tank Radius',tank_radius,distance_unit;
    linedivider,'','';
    'Throat Conditions','','';
    linedivider,'','';
    'Throat Area',throat_area,area_unit;
    'Throat Radius',throat_radius,distance_unit;
    'Throat Temperature',throat_temperature,temperature_unit;
    'Throat Pressure',throat_pressure,pressure_unit;
    'Throat Density',throat_density,density_unit;
    'Throat Velocity',throat_velocity,velocity_unit;
    linedivider,'','';
    'Exit Conditions','','';
    linedivider,'','';
    'Exit Temperature',exit_temperature,temperature_unit;
    'Exit Pressure',exit_pressure,pressure_unit;
    'Exit Mach',exit_mach,unitless;
    'Exhaust Velocity',exhaust_velocity,velocity_unit;
    'Effective Exhaust Velocity',effective_exhaust_velocity,velocity_unit;
    linedivider,'','';
    'System Conditions','','';
    linedivider,'','';
    'Mass Flowrate I',mass_flowrate_I,massflow_unit;
    'Mass Flowrate II',mass_flowrate_II,massflow_unit;
    'Mass Flowrate III',mass_flowrate_III,massflow_unit;
    'Propellant Mass',propellant_mass,mass_unit;
    'Temperature Ratio',temperature_ratio,unitless;
    'Pressure Ratio',pressure_ratio,unitless;
    'Thrust Coefficient',thrust_coefficient,unitless;
    'Thrust Force I',thrust_force_I,force_unit;
    'Thrust Force II',thrust_force_II,force_unit;
    'Thrust Force III',thrust_force_III,force_unit;
    'Thrust Force IV',thrust_force_IV,force_unit;
    'Specific Impulse I',specific_impulse_I,time_unit;
    'Specific Impulse II',specific_impulse_II,time_unit
    'Final Mass',final_mass,mass_unit;
    'Initial Mass',initial_mass,mass_unit;
    'Mass Ratio',mass_ratio,unitless;
    'Propellant Mass',propellant_mass,mass_unit;
    'Burn Time',burn_time,time_unit;
    linedivider,'','';
    'Iterated Values','','';
    linedivider,'','';
    'Total Impulse',total_impulse,newton_time;
    'Specific Impulse',specific_impulse,newton_time;
    'Delta Velocity',delta_velocity,velocity_unit};
    
display(result)

exceldata = {'Variable','Magnitude','Unit';'Propellant',PropellantName,unitless;'Specific Gas Constant',specific_gas_constant,gas_constant_unit;...
    'Molar Mass',molar_mass,molar_mass_unit;'Specific Heat Ratio',specific_heat_ratio,unitless;'Nozzle Conditions',linedivider,linedivider;...
    'Nozzle Length',nozzle_length,distance_unit;'Exit Radius',exit_radius,distance_unit;'Throat Radius',throat_radius,distance_unit;...
    'Exit Area',exit_area,area_unit;'Throat Area',throat_area,area_unit;'Chamber Conditions',linedivider,linedivider;...
    'Tank Geometry',geometry,unitless;'Tank Volume',volume,volume_unit;'Chamber Pressure',chamber_pressure(1),pressure_unit;...
    'Throat Conditions',linedivider,linedivider;...
    'Throat Area',throat_area,area_unit;'Throat Radius',throat_radius,distance_unit;'Throat Temperature',throat_temperature,temperature_unit;...
    'Throat Pressure',throat_pressure,pressure_unit;'Throat Density',throat_density,density_unit;'Throat Velocity',throat_velocity,velocity_unit;...
    'Exit Conditions',linedivider,linedivider;'Exit Temperature',exit_temperature,temperature_unit;'Exit Pressure',exit_pressure,pressure_unit;...
    'Exit Mach',exit_mach,unitless;'Exhaust Velocity',exhaust_velocity,velocity_unit;'Effective Exhaust Velocity',effective_exhaust_velocity,velocity_unit;...
    'System Preformance',linedivider,linedivider;'Mass Flowrate I',mass_flowrate_I,massflow_unit;'Mass Flowrate II',mass_flowrate_II,massflow_unit;...
    'Mass Flowrate III',mass_flowrate_III,massflow_unit;'Propellant Mass',propellant_mass,mass_unit;'Temperature Ratio',temperature_ratio,unitless;...
    'Pressure Ratio',pressure_ratio,unitless;'Thrust Coefficient',thrust_coefficient,unitless;'Thrust Force I',thrust_force_I,force_unit;...
    'Thrust Force II',thrust_force_II,force_unit;'Thrust Force III',thrust_force_III,force_unit;'Thrust Force IV',thrust_force_IV,force_unit;...
    'Specific Impulse I',specific_impulse_I,time_unit;'Specific Impulse II',specific_impulse_II,time_unit;'Final Mass',final_mass,mass_unit;...
    'Initial Mass',initial_mass,mass_unit;'Mass Ratio',mass_ratio,unitless;'Burn Time',burn_time,time_unit};
xlswrite('Spacecraft Propulsion System Data',exceldata);

end

function Mach = solve_mach(A,At,k)
%   solve Mach number from area ratio by Newton-Raphson Method. (assume
%   supersonic)
%   https://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
%   Given by Phil Lindin: Thruster Design Utility
%   Github: https://github.com/runphilrun/TDU
P = 2/(k+1);
Q = 1-P;
R = (A/At).^((2*Q)/P);
a = Q.^(1/P);
r = (R-1)/(2*a);
X = 1/((1+r)+sqrt(r*(r+2)));  % initial guess
diff = 1;  % initalize termination criteria
while abs(diff) > .0001
    F = (P*X+Q).^(1/P)-R*X;
    dF = (P*X+Q).^((1/P)-1)-R;
    Xnew = X - F/dF;
    diff = Xnew - X;
    X = Xnew;
end
Mach = 1/sqrt(X);
end

































