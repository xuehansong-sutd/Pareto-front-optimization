clc
clear
% %% Top Cell Post Processing
% % Top Cell
load TC_IV.txt
n = size(TC_IV,1);
Vap_COMSOL_TC = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.1 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.2 1.21 1.22 1.23 1.24 1.25];

% Parameter input

n_Vap_TC = size(Vap_COMSOL_TC,2); % Number of applied voltage points simulated in COMSOL

CalculatedIV_TC = TC_IV';

%% Post processing of COMSOL IV curve
% Vap is applied voltage used to simulate in MCS

hinterp = 0.01; % Step size for interpolation
Vap_interp_TC = Vap_COMSOL_TC; % Vap points to be interpolated

% Extract IV curve to obtain Jsc, Voc, FF, Eff

Jsc_TC = zeros(1,n); % A/m^2
Voc_TC = zeros(1,n); % V
FF_TC = zeros(1,n);
Eff_TC = zeros(1,n);

for i=1:1:n

    j_TC = CalculatedIV_TC(1:n_Vap_TC,i)./10; % Convert A/m^2 to mA/cm^2
    %     jinterp_TC = interp1(Vap_COMSOL_TC,j_TC,Vap_interp_TC);
    jinterp_TC = j_TC';

    Jsc_TC(:,i) = jinterp_TC(end);
%     Voc_TC(:,i) = interp1(jinterp_TC(1,end-15:end),Vap_interp_TC(1,end-15:end), 0, 'linear');
    Pmax_TC = max(Vap_interp_TC.*jinterp_TC); % mW/(cm^2)
    FF_TC(:,i) = Pmax_TC/(Jsc_TC(:,i)*Voc_TC(:,i));
    Eff_TC(:,i) = Pmax_TC/100; % Divided by one sun which is 100 mW/(cm^2)

end

Collated_TC = [Jsc_TC;
    Voc_TC;
    FF_TC;
    Eff_TC]; % Collated data used for post analysis in folder 3b.

save('c2_Collated_TC.mat','Collated_TC');
%% Bottom Cell
load BC_IV.txt
n = size(BC_IV,1);
Vap_COMSOL_BC = 0:0.01:0.7;

% Parameter input

n_Vap_BC = size(Vap_COMSOL_BC,2); % Number of applied voltage points simulated in COMSOL

CalculatedIV_BC = BC_IV';

%% Post processing of COMSOL IV curve
% Vap is applied voltage used to simulate in MCS

hinterp = 0.01; % Step size for interpolation
Vap_interp_BC = Vap_COMSOL_BC(1):hinterp:Vap_COMSOL_BC(end); % Vap points to be interpolated

% Extract IV curve to obtain Jsc, Voc, FF, Eff

Jsc_BC = zeros(1,n); % A/m^2
Voc_BC = zeros(1,n); % V
FF_BC = zeros(1,n);
Eff_BC = zeros(1,n);

for i=1:1:n

    j_BC = CalculatedIV_BC(1:n_Vap_BC,i)./10; % Convert A/m^2 to mA/cm^2
    %     jinterp_TC = interp1(Vap_COMSOL_TC,j_TC,Vap_interp_TC);
    jinterp_BC = j_BC';

    Jsc_BC(:,i) = jinterp_BC(1);
    if sum(jinterp_BC(1,1:10)) < 10
        Voc_BC(:,i) = 0;
    else
        Voc_BC(:,i) = interp1(jinterp_BC(1,end-15:end),Vap_interp_BC(1,end-15:end), 0, 'linear');
    end
    Pmax_BC = max(Vap_interp_BC.*jinterp_BC); % mW/(cm^2)
    FF_BC(:,i) = Pmax_BC/(Jsc_BC(:,i)*Voc_BC(:,i));
    Eff_BC(:,i) = Pmax_BC/100; % Divided by one sun which is 100 mW/(cm^2)

end

Collated_BC = [Jsc_BC;
    Voc_BC;
    FF_BC;
    Eff_BC]; % Collated data used for post analysis in folder 3b.

save('c2_Collated_BC.mat','Collated_BC');

Collated_tot = Eff_TC+Eff_BC;

save('c2_Collated_tot.mat','Collated_tot');