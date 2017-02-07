function varagout = spsd(system_velocity,inlet_temperature,inlet_pressure,chamber_geometry,exit_radius,area_ratio,propellant,altitude,initial_mass)
% Title: Spacecraft Propulsion System Design
% Purpose: This script takes a series of mostly physical atributes and
% calculates the performace of the propulsion system.
% Author: James Emerson Parkus
% Date: December 26th, 2016
% Notes: See the document that goes over this script. It will explain how
% the chamber geometry, propellant, and altitude factors are determined.
% Gas List:
% 1 - Hydrogen
% 2 - Helium
% 3 - Nitrogen
% 4 - Freon [R-22]
% 5 - Neon
% 6 - Xenon
% 7 - Argon
% 
% Function Arguement Example = 122 m/s, 300 K, 1*10^6 Pa, 2 [rectangular],
% 0.005 m {diameter of 1 cm), 100, 3 {GN2}, 4*10^5 m [Standard LEO], 2.66kg

clc;

%% Constants
standard_gravity = 9.8066; %m/s^2
universal_gas_constant = 8.3144598; % [J/(mol*K)]
earth_radius = 6.37*10^6;
adjusted_gravity = standard_gravity*(earth_radius/(earth_radius+altitude))^2;

%% Nozzle Conditions
exit_area = pi*exit_radius^2;

throat_area = exit_area/area_ratio;

throat_radius = sqrt(throat_area/pi);

conical_half_angle = 15;

nozzle_length = (exit_radius-throat_radius)/tan(deg2rad(conical_half_angle)); % Length for diverging section

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
propellant_data = xlsread('RPSD Propellant Info.xlsx'); % Pulls data from excel file with propellant information
propellant_name_data = propellant_data(:,2);
molar_mass_data = propellant_data(:,3);
heat_ratio_data = propellant_data(:,4);
specific_heat_constant_pressure_data = propellant_data(:,6);

switch propellant % Pulls and assigns the name to the propellant name variable
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
    case 8
        PropellantName = 'Carbon Dioxide';
end

molar_mass = molar_mass_data(propellant,1); % Assigns the molar mass of the selected gas
specific_heat_ratio = heat_ratio_data(propellant,1); % Assigns the specific heat ratio of the selected gas
specific_heat_constant_pressure = specific_heat_constant_pressure_data(propellant,1);

specific_gas_constant = universal_gas_constant/molar_mass;

%% Tank Conditions
switch chamber_geometry % finds the volume of the selected tank [Does not work for toroidal tank, 1/4/17]
    case 1
        geometry = 'Spherical';
        radius = 0.05;
        volume = 4/3*pi*radius^3;
        propellant_mass(1) = inlet_pressure*volume/(specific_gas_constant*inlet_temperature);
    case 2
        geometry = 'Rectangular';
        length = 0.1;
        height = 0.1;
        width = 0.05;
        volume = length*width*height;
        propellant_mass(1) = inlet_pressure*volume/(specific_gas_constant*inlet_temperature);
    case 3
        geometry = 'Cylindrical';
        height = 0.1;
        radius = 0.05;
        volume = pi*radius^2*height;
        propellant_mass(1) = inlet_pressure*volume/(specific_gas_constant*inlet_temperature);
    case 4
        geometry = 'Toroidal';
        minor_radius = 1;
        major_radius = 1;
        volume = (pi*minor_radius^2)*(2*pi*major_radius);
        propellant_mass(1) = inlet_pressure*volume/(specific_gas_constant*inlet_temperature);
end

% propellant_mass(1) = 6.01; % Just for testing published data
% volume = 0.01; % Testing published data
    %% System Conditions
    
    exhaust_velocity(1) = sqrt(2*specific_heat_constant_pressure*inlet_temperature); % maximum exit velocity for infinite expansion [expansion into vacuum]
    
    final_mass = initial_mass - propellant_mass(1);
    
    propellant_mass(1) = initial_mass - final_mass;
    
    mass_ratio = initial_mass/final_mass;
    %% Throat Conditions
    
    inlet_pressure(1) = inlet_pressure;
    
    throat_temperature = inlet_temperature*(2/(specific_heat_ratio+1));
    
    throat_pressure(1) = inlet_pressure(1)*(2/(specific_heat_ratio+1))^(specific_heat_ratio/(specific_heat_ratio-1));
    
    throat_density(1) = throat_pressure(1)/(specific_gas_constant*throat_temperature);
    
    throat_velocity = sqrt(specific_heat_ratio*specific_gas_constant*throat_temperature);
    
    mass_flowrate(1) = throat_density(1)*throat_velocity*throat_area;
    
    %% Exit Conditions
    % exit_mach = solve_mach(exit_area,throat_area,specific_heat_ratio); % Executes the newton-raphson approximation function to find exit mach
         
    exit_pressure(1) = inlet_pressure(1)*((specific_heat_ratio-1)/2*(1-((specific_heat_ratio-1)/(specific_heat_ratio+1))*(exhaust_velocity(1)*exit_area/(throat_velocity*throat_area))^2))^(specific_heat_ratio/(specific_heat_ratio-1));
    
    exit_mach = sqrt(2*((inlet_pressure(1)/exit_pressure(1))^((specific_heat_ratio-1)/specific_heat_ratio)-1)/(specific_heat_ratio-1));
    
    exit_temperature = inlet_temperature/(1+(specific_heat_ratio-1)/2*exit_mach^2);
    
    temperature_ratio = exit_temperature/inlet_temperature;
%% Iteration Sequence
dt = 1*10^-3; % Time Step [0.001 seconds]
i = 1; % Step Count
t(1) = 0; % Simulated Time
lim = 1*10^4; % Pressure limit [100000 Pa]
% inlet_pressure(i = inlet_pressure;

while inlet_pressure >= lim
         
%     propellant_mass(i) = propellant_mass(1);
    
    dm(i) = mass_flowrate(i)*dt; % Calculates the ejected mass per unit dt
     
    propellant_mass(i+1) = propellant_mass(i) - dm(i); % Updates the amount of propellant in tank 
    
    inlet_pressure(i+1) = propellant_mass(i+1)*specific_gas_constant*inlet_temperature/volume; % Finds the new inlet pressure
    
    exit_pressure(i+1) = inlet_pressure(i+1)/((temperature_ratio)^(specific_heat_ratio/(specific_heat_ratio-1))); % Finds the new exit pressure
    
    thrust_coefficient(i+1) = sqrt(2*specific_heat_ratio^2/(specific_heat_ratio-1))*(2/(specific_heat_ratio+1))^((specific_heat_ratio+1)/(specific_heat_ratio-1))...
        *(1-(exit_pressure(i+1)/inlet_pressure(i+1))^((specific_heat_ratio-1)/specific_heat_ratio))+(exit_pressure(i+1)/inlet_pressure(i+1))*(area_ratio); % Finds the updated thrust coefficient
    
    thrust_force(i+1) = thrust_coefficient(i+1)*inlet_pressure(i+1)*throat_area; % Calculates the thrust using the updated values
    
    acceleration(i) = thrust_force(i+1)*dt*1/dm(i); % Solving the F = m*a for a
    
    velocity(i) = acceleration(i)*dt; % Attempt to iteratively find the system velocity [Does not seem to be working, 1/4/17]
    
%     distance(i+1) = velocity(i+1)*dt; % Just for fun...finds the linear distance traveled [Should be converted to radial distance for real usage]
    
    throat_pressure(i+1) = inlet_pressure(i+1)*(2/(specific_heat_ratio+1))^(specific_heat_ratio/(specific_heat_ratio-1));
    
    throat_density(i+1) = throat_pressure(i+1)/(specific_gas_constant*throat_temperature);
    
    mass_flowrate(i+1) = throat_density(i+1)*throat_velocity*throat_area;
    
    impulse_bit(i) = thrust_force(i+1)*dt; % Finds the impulse per unit dt
    velocity_bit(i) = velocity(i)*dt; % Finds the velocity per unit dt
    
    impulse_bit(i) = thrust_force(i)*dt; % Finds the impulse per unit dt
    velocity_bit(i) = velocity(i); % Finds the velocity per unit dt
    massflow_bit(i) = mass_flowrate(i)*dt;
    
    t(i+1) = t(i) + dt; % Updates the simulation times
    i = i + 1; % Increases the step count
end
total_mass_flowrate = sum(massflow_bit);
total_impulse = sum(impulse_bit);
specific_impulse = total_impulse/(standard_gravity*total_mass_flowrate);

delta_velocity = sum(velocity_bit);

max_thrust_force = max(thrust_force);
max_acceleration = max(acceleration);
max_exhaust_velocity = max(exhaust_velocity);

grid on
title(PropellantName)
plot(t(1:i),thrust_force(1:i))
xlabel('Time [s]')
ylabel('Thrust [N]')


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
    'Inlet Conditions','','';
    linedivider,'','';
    'Tank Geometry','',geometry;
    'Volume',volume,volume_unit;
    'Inlet Pressure',inlet_pressure(1),pressure_unit;
    'Inlet Temperature',inlet_temperature,temperature_unit;
    linedivider,'','';
    'Throat Conditions','','';
    linedivider,'','';
    'Throat Area',throat_area,area_unit;
    'Throat Radius',throat_radius,distance_unit;
    'Throat Temperature',throat_temperature(1),temperature_unit;
    'Throat Pressure',throat_pressure(1),pressure_unit;
    'Throat Density',throat_density(1),density_unit;
    'Throat Velocity',throat_velocity,velocity_unit;
    linedivider,'','';
    'Exit Conditions','','';
    linedivider,'','';
    'Exit Temperature',exit_temperature,temperature_unit;
    'Exit Pressure',exit_pressure(1),pressure_unit;
    'Exit Mach',exit_mach,unitless;
    'Exhaust Velocity',max_exhaust_velocity,velocity_unit;
    linedivider,'','';
    'System Conditions','','';
    linedivider,'','';
    'Mass Flowrate',mass_flowrate(1),massflow_unit;
    'Propellant Mass',propellant_mass(1),mass_unit;
    'Thrust Coefficient',thrust_coefficient(1),unitless;
    'Thrust Force',thrust_force(1),force_unit;
    'Final Mass',final_mass,mass_unit;
    'Initial Mass',initial_mass(1),mass_unit;
    'Mass Ratio',mass_ratio,unitless;
    'Propellant Mass',propellant_mass(1),mass_unit;
    %     'Burn Time',burn_time,time_unit;
    'Adjusted Gravity',adjusted_gravity,acceleration_unit;
    linedivider,'','';
    'Published Data Testing','','';
    linedivider,'','';
    'Delta Velocity',delta_velocity,velocity_unit;
    linedivider,'','';
    'Iterated Values','','';
    linedivider,'','';
    'Total Impulse',total_impulse,newton_time;
    'Specific Impulse',specific_impulse,time_unit;
    'Delta Velocity',delta_velocity,velocity_unit;
    'Max Thrust',max_thrust_force,force_unit;
    'Max Acceleration',max_acceleration,acceleration_unit;
    'Mass Flowrate',mass_flowrate(1),massflow_unit;
    'Iterations',i,unitless;};

display(result);

xlswrite('Spacecraft Propulsion System Data',result);

end

% function Mach = solve_mach(A,At,k)
% %   solve Mach number from area ratio by Newton-Raphson Method. (assume
% %   supersonic)
% %   https://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
% %   Courtesy of Phil Lindin's Thruster Design Utility
% %   Github: https://github.com/runphilrun/TDU
% P = 2/(k+1);
% Q = 1-P;
% R = (A/At).^((2*Q)/P);
% a = Q.^(1/P);
% r = (R-1)/(2*a);
% X = 1/((1+r)+sqrt(r*(r+2)));  % initial guess
% diff = 1;  % initalize termination criteria
% while abs(diff) > .0001
%     F = (P*X+Q).^(1/P)-R*X;
%     dF = (P*X+Q).^((1/P)-1)-R;
%     Xnew = X - F/dF;
%     diff = Xnew - X;
%     X = Xnew;
% end
% Mach = 1/sqrt(X);
%
% end

function display(result)
% Courtesy of Phil Lindin's TDU
[n,~]=size(result);
for i = 1:n
    fprintf('\n%24s\t%12f\t%s',result{i,:});
end
end
























