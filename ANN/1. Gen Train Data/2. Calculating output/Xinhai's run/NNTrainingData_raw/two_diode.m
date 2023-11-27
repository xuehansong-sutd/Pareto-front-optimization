%%  Two-diode model function ----------------------------------------------
%   Author: Tan Hu Quee
%   Date Created: 11th Jan 2022
%   Date Modified:
%--------------------------------------------------------------------------
function CalculatedIV = two_diode(Jsc)
%% Initial guess for fsolve
% load IV_test.txt
V = 0:0.01:0.7;  
IV_test = zeros(size(V,2),1); 
% Jsc = Jsc*ones(length(IV_test(:,1),1))
%% Two-diode equation
j = fsolve(@(j)two_diode_eqn(j,Jsc,V),IV_test(:,1));
CalculatedIV = [V',j];
end

function F = two_diode_eqn(j,Jsc,V)
%% User input for fitted parameters for two-diode model
J01 = 1.52E-08;
J02 = 0.00173;
n1 = 1;
n2 = 2;
Rs = 6.30E-05;
Rshunt = 0.179;
%% Constants
kB = 1.38064852e-23;                               % Boltzmann constant in m^2 kg s^-2 K^-1
% kB = 8.617e-5;                                     % Boltzmann constant in eV K^-1
T = 298;                                           % Temperature in K
q = 1.602e-19;                                     % Charge in C
%% Two-diode equation
J1 = J01*(exp((q.*(V'+j*Rs)/(n1*kB*T)))-1);
J2 = J02*(exp((q.*(V'+j*Rs)/(n2*kB*T)))-1);
Jshunt = (V'+j.*Rs)/Rshunt;
F = Jsc-J1-J2-Jshunt-j;
end