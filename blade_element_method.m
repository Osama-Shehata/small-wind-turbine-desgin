function BEM_Equation = blade_element_method(a,blade_solidity,blade_twist,local_tip_speed_ratio,pitch_angle,alpha_data,lift_coefficients,drag_coefficients)
%Input Parameters:
%a: a vector of length two representing the axial and tangential induction factors, respectively.
%blade_solidity: the solidity of the blade element.
%blade_twist: the twist angle of the blade element.
%local_tip_speed_ratio: the local tip speed ratio of the blade element.
%pitch_angle: the pitch angle of the blade element.
%alpha_data: a vector of angle of attacks for the airfoil.
%lift_coefficients: a vector of lift coefficients for the airfoil.
%drag_coefficients: a vector of drag coefficients for the airfoil.

%Output:
%The output of the function is a vector containing the difference between the CP and CT calculated using momentum theory and blade element theory, respectively.
%%
% Calculate CP and CT from momentum theory
    cp_momentum = 4 * a(2) * (1 - a(1)) * local_tip_speed_ratio^2;
    a_star = 0.37;
    if a(1) > a_star
        CT_M = 4 * a_star * (1 - a_star) + 4 * (1 - 2 * a_star) * (a(1) - a_star);
    else
        CT_M = 4 * a(1) * (1 - a(1));
    end
%%
% Calculate phi and alpha
phi = atan((1 - a(1)) / ((1 + a(2)) * local_tip_speed_ratio));
a_of_attack = phi * 180 / pi - blade_twist - pitch_angle;
%%
% Interpolate CL and CD
CL_alpha = interp1(alpha_data,lift_coefficients,a_of_attack); 
CD_alpha = interp1(alpha_data,drag_coefficients,a_of_attack); 
%%
% Calculate CX and CY
cx = CL_alpha * sin(phi) - CD_alpha * cos(phi);
cy = CL_alpha * cos(phi) + CD_alpha * sin(phi);
%%
% Calculate CP and CT from blade element theory
W_U_squared = (1 - a(1))^2 + local_tip_speed_ratio^2 * (1 + a(2))^2;
CP_BE = blade_solidity * W_U_squared * local_tip_speed_ratio * cx;
CT_BE = blade_solidity * W_U_squared * cy;
%%
% Calculate delta_CP and delta_CT
delta_CP = cp_momentum - CP_BE;
delta_CT = CT_M - CT_BE;
%%
BEM_Equation = [delta_CP;delta_CT];
end
