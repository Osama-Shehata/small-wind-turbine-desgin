function [N_blades, rotor_radius, avg_blade_solidity, avg_chord_length, root_chord_length, tip_chord_length, tip_speed_ratios, blade_sections, pitch_angles, alpha_data, lift_coefficients, drag_coefficients , Coefficient_Of_Power , Coefficient_Of_Torque] = Initialize_Blade_Variables()
%This function initializes and sets various variables related to the blade geometry and aerodynamics of a wind turbine rotor. 
% The function defines the number of blades, rotor radius, blade solidity, chord lengths, tip speed ratios, blade sections, pitch angles, 
% and loads data from an external file containing alpha values, lift coefficients, and drag coefficients. 
% The function also initializes matrices for storing results related to the coefficient of power and torque.
%%
N_blades = 3;
rotor_radius = 20;
avg_blade_solidity = 0.03;
avg_chord_length = avg_blade_solidity*pi*rotor_radius/N_blades;
root_chord_length = 1.06*avg_chord_length;
tip_chord_length = 0.94*avg_chord_length;
tip_speed_ratios = 1:0.5:16;
blade_sections = 1:0.2:20;
pitch_angles = -15:2:5;

%%
% Initialize variables for storing results
Coefficient_Of_Power = zeros(length(pitch_angles), length(tip_speed_ratios));
Coefficient_Of_Torque = zeros(length(pitch_angles), length(tip_speed_ratios));

%%
alpha_data = table2array(readtable('rotor.xlsx','Range','A1:A54'));
lift_coefficients = table2array(readtable('rotor.xlsx','Range','B1:B54'));
drag_coefficients = table2array(readtable('rotor.xlsx','Range','C1:C54'));


end





