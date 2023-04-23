clear;
clc;

%%
%Osama Mahmoud Shehata - 201900776
%Abdelrahaman Tehamer - 201900735

%Rotor Aerodynamics - Project #1
%Under Supervision of: Dr Mostafa Abdallah Gabr 
%Blade Element Momentum analysis of a wind turbine rotor with three blades.

%%
%1- Initialize Blade Variables
%2- For each pitch angle, then for each tip speed ratios loop over
%different blade sections to to calculate the local tip speed ratio, twist angle, chord length, and solidity for each blade section. 
%3- Use a blade element momentum (blade_element_method) function to solve for the induction factor at each blade section
%4- The local power coefficients are then integrated over the blade to get the overall power coefficient and torque coefficient of the rotor
%5- Plot the Coefficient Of Torque and Coefficient Of Power as a function of the tip speed ratio for each pitch angle.
%6- Calculate the tip speed ratio for a wind speed of 8 m/s and rotor speed of 40 rpm

%%
% Initialize blade variables
[N_blades, rotor_radius, avg_blade_solidity, avg_chord_length, root_chord_length, tip_chord_length, tip_speed_ratios, blade_sections, pitch_angles, alpha_data, lift_coefficients, drag_coefficients , Coefficient_Of_Power , Coefficient_Of_Torque] = Initialize_Blade_Variables();

%%
% Iterate through each pitch angle
for b = 1:length(pitch_angles)
    pitch_angle = pitch_angles(b);
    
    % Initialize variables for storing local results
    A = [];
    Coefficient_Of_Power_Local = zeros(size(blade_sections));
    aa = zeros(length(blade_sections), 2);
    
    % Iterate through each tip speed ratio
    for lambda_n = 1:length(tip_speed_ratios)
        lambda = tip_speed_ratios(lambda_n);
        
        % Iterate through each blade section
        for i = 1:length(blade_sections)
            r = blade_sections(i);
            
            % Calculate local variables
            local_tip_speed_ratio = lambda*(r/rotor_radius);
            blade_twist = 27 - (27/20)*r;
            C_r = ((tip_chord_length - root_chord_length)*(r/rotor_radius) + root_chord_length);
            blade_solidity = N_blades*C_r/(2*pi*r);
            
            % Solve for induction factors using blade element momentum theory
            F = @(a) blade_element_method(a, blade_solidity, blade_twist, local_tip_speed_ratio, pitch_angle, alpha_data, lift_coefficients, drag_coefficients);
            x0 = [0.1 0.1];
            a = fsolve(F, x0);
            aa(i,1) = a(1);
            aa(i,2) = a(2);
            
            % Calculate local power coefficient using momentum theory
            Coefficient_Of_Power_Local(i) = 4*aa(i,2)*(1-aa(i,1))*local_tip_speed_ratio^2;
        end
        
        % Integrate local power coefficient over the blade to get total power coefficient
        Coefficient_Of_Power(b, lambda_n) = 2*trapz(blade_sections, blade_sections.*(Coefficient_Of_Power_Local)/rotor_radius^2);
        
        % Integrate local torque coefficient over the blade to get total torque coefficient
        Coefficient_Of_Torque(b, lambda_n) = 2*trapz(blade_sections, blade_sections.*(Coefficient_Of_Power_Local)/(rotor_radius^2*lambda));
    end
    
    % Plot results for each pitch angle
    figure(1)
    plot(tip_speed_ratios,Coefficient_Of_Torque(b,:))
    xlabel('lambda')
    ylabel('Coefficient Of Torque')
    hold on
    figure(2)
    plot(tip_speed_ratios,Coefficient_Of_Power(b,:))
    xlabel('lambda')
    ylabel('Coefficient Of Power')
    hold on
end

%%
% Plot final results
figure(1)
title('Coefficient Of Torque vs lambda @ different pitch angles')
legend('-15','-13','-11','-9','-7' , '-5' , '-3' , '-1' , '1' , '3' , '5')
figure(2)
title('Coefficient Of Power vs lambda @ different pitch angles')
legend('-15','-13','-11','-9','-7' , '-5' , '-3' , '-1' , '1' , '3' , '5')


%%
% Calculate power at 8 m/s and 40 rpm
lambda_m = 40*2*pi*rotor_radius/(8*60);
maximum = max(Coefficient_Of_Power(:,20));
[bb, ll] = find(Coefficient_Of_Power(:,20)==maximum);
CP_m = max(interp1(tip_speed_ratios,Coefficient_Of_Power(bb,:),lambda_m));
P = 0.5*1.2*pi*(rotor_radius^2)*(8^3)*CP_m; 