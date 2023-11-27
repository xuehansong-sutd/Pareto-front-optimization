function Gx = TMMTC(x,L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,LglassTC,xi)

% Load in 1sun AM 1.5 solar spectrum in mW/cm2nm
load AM15.mat
L1=L1(1,1);
L2=L2(1,1);
L3=L3(1,1);
L4=L4(1,1);
L5=L5(1,1);
L6=L6(1,1);
L7=L7(1,1);
L8=L8(1,1);
L9=L9(1,1);
L10=L10(1,1);
L11=L11(1,1);
L12=L12(1,1);
L13=L13(1,1);
LglassTC = LglassTC(1,1);
xi = xi(1,1);

% Determine x-coordinates to be solved in COMSOL in [m]
[row, column] = size(x);
x_correction = L2.*ones(row,column);
x = x-x_correction;

% Wavelength to be solved in[nm]
lambda_start = 300;
lambda_end = 1400;
no_interval = 1100;
stepsize = (lambda_end-lambda_start)/no_interval;
lambda = lambda_start:stepsize:lambda_end;

no_lambda = size(lambda,2);

% Loading in of optical properties
layers = {'SodaLime' 'ITO' 'NiO' 'PerovHMv2' 'C60HM' 'SnO2' 'ITO' 'LiF'...
    'MgF2' 'AZO' 'ZnO' 'CdS' 'CIS' 'Mo' 'SodaLime' 'Air'};

%------------END USER INPUT PARAMETERS SPECIFICATION-------------------

% Load in index of refraction for each material
n = zeros(size(layers,2),size(lambda,2));
for index = 1:size(layers,2)
    n(index,:) = LoadRefrIndex(layers{index},lambda);
end

% Constants
h = 6.62606957e-34; 	                    % Js Planck's constant
c = 2.99792458e8;	                        % m/s speed of light
q = 1.60217657e-19;	                        % C electric charge
epsilon0 = 8.854e-12;                       % F/m permittivity
perov_n = real(n(4,:))+real(n(4,:)).*xi;   % Fitting parameter n
perov_k = imag(n(4,:))-imag(n(4,:)).*xi;   % Fitting parameter \kappa
n(4,:) = complex(perov_n,perov_k);

%% Calculating the correction terms

% TOP CELL
M11 = zeros(1,no_lambda);
M12 = zeros(1,no_lambda);
M21 = zeros(1,no_lambda);
M22 = zeros(1,no_lambda);
B11 = zeros(1,no_lambda);
B12 = zeros(1,no_lambda);
B21 = zeros(1,no_lambda);
B22 = zeros(1,no_lambda);
GsumnTC = zeros(row,no_lambda);
lam = 300e-9;
AM15_count = 1;

for i = 1:no_lambda
    
    I_0 = I_mat(n(1,i),n(2,i));
    L_1 = L_mat(n(2,i),L1,lam);
    I_1 = I_mat(n(2,i),n(3,i));
    L_2 = L_mat(n(3,i),L2,lam);
    I_2 = I_mat(n(3,i),n(4,i));
    L_3 = L_mat(n(4,i),L3,lam);
    I_3 = I_mat(n(4,i),n(5,i));
    L_4 = L_mat(n(5,i),L4,lam);
    I_4 = I_mat(n(5,i),n(6,i));
    L_5 = L_mat(n(6,i),L5,lam);
    I_5 = I_mat(n(6,i),n(7,i));
    L_6 = L_mat(n(7,i),L6,lam);
    I_6 = I_mat(n(7,i),n(8,i));
    L_7 = L_mat(n(8,i),L7,lam);
    I_7 = I_mat(n(8,i),n(9,i));

    L_8 = L_mat(n(9,i),L8,lam);
    I_8 = I_mat(n(9,i),n(10,i));
    L_9 = L_mat(n(10,i),L9,lam);
    I_9 = I_mat(n(10,i),n(11,i));
    L_10 = L_mat(n(11,i),L10,lam);
    I_10 = I_mat(n(11,i),n(12,i));
    L_11 = L_mat(n(12,i),L11,lam);
    I_11 = I_mat(n(12,i),n(13,i));
    L_12 = L_mat(n(13,i),L12,lam);
    I_12 = I_mat(n(13,i),n(14,i));
    L_13 = L_mat(n(14,i),L13,lam);
    I_13 = I_mat(n(14,i),n(15,i));
  
% Top Cell matrix
    M = I_0*L_1*I_1*L_2*I_2*L_3*I_3*L_4*I_4*L_5*I_5*L_6*I_6*L_7*I_7*...
        L_8*I_8*L_9*I_9*L_10*I_10*L_11*I_11*L_12*I_12*L_13*I_13;
    M11(1,i)=M(1,1);
    M12(1,i)=M(1,2);
    M21(1,i)=M(2,1);
    M22(1,i)=M(2,2);
    R12=abs(M(2,1)/M(1,1))^2;
    T12=abs(1/M(1,1))^2;
    R21=abs(M(1,2)/M(1,1))^2;
    T21=abs(M(2,2)-(M(2,1)*M(1,2)/M(1,1)))^2;

% Incoherent correction term calculation
    R01 = abs((n(end,i)-n(1,i))/(n(end,i)+n(1,i)))^2;
    T01 = abs(2*n(end,i)/(n(end,i)+n(1,i)))^2;
    R10 = abs((n(1,i)-n(end,i))/(n(end,i)+n(1,i)))^2;
    T10 = abs(2*n(1,i)/(n(end,i)+n(1,i)))^2;
%     E1f_sq = abs(sqrt((2*AM15(1,i))/(c*epsilon0)))^2;
    P1 = exp((4*pi*imag(n(1,i))*LglassTC)/lam);
    M11_til = -(- P1^2 + R10*R12)/(P1*T01*T12);
    M12_til = -(R21*P1^2 - R10*R12*R21 + R10*T12*T21)/(P1*T01*T12);
    R(1,i) = -(R01*P1^2 - R01*R10*R12 + R12*T01*T10)/(- P1^2 + R10*R12);
    T(1,i) = -(P1*T01*T12)/(- P1^2 + R10*R12);
% Correction term
    Cf = -(P1*T01)/(- P1^2 + R10*R12);    
% Top cell electric field calculation
    Ef = sqrt((2*Cf*AM15(1,i))/(c*epsilon0));
    Mf = I_0*L_1*I_1*L_2*I_2;
    Mb = L_3*I_3*L_4*I_4*L_5*I_5*L_6*I_6*L_7*I_7*L_8*I_8*L_9*I_9*L_10*I_10*L_11*I_11*L_12*I_12;
    
    Sf11 = Mf(1,1);
    Sf12 = Mf(1,2);
    Sf21 = Mf(2,1);
    Sf22 = Mf(2,2);
    
    Sb11 = Mb(1,1);
    Sb12 = Mb(1,2);
    Sb21 = Mb(2,1);
    Sb22 = Mb(2,2);    
    
    v3f = ((Sb11/(Sb11*Sf11+Sf12*Sb21)))*Ef;
    w3f = ((Sb21/(Sb11*Sf11+Sf12*Sb21)))*Ef;
    
    k_perovTC = (2*pi/lam)*n(4,i);
    nk_perovTC = real(n(4,i))*imag(n(4,i));
    
    EzF(:,1) = v3f.*exp(1i.*k_perovTC.*x(:,1))+w3f.*exp(-1i.*k_perovTC.*x(:,1));
    
    Ez_mag(:,1) = abs(EzF(:,1)).^2;

    GsumnTC(:,i) = ((2*pi*nk_perovTC)/(h))*epsilon0*Ez_mag(:,1);
    
%     GsumnTC(:,i) = 10*AM15(1,AM15_count)*(4*pi*nk_perovTC/(h*c))*Ez_mag(:,1);
    
    lam = lam+stepsize*1e-9;
    AM15_count = AM15_count+stepsize;
    
end

Gx = zeros(row,1);
  
f_0 = GsumnTC(:,1);
f_n = GsumnTC(:,end);

f_even = zeros(row,1);
f_odd = zeros(row,1);

    for j = 2:1:no_lambda-1
        if mod(j,2)==0
            f_odd(:,1) = f_odd(:,1)+GsumnTC(:,j);
        else
            f_even(:,1) = f_even(:,1)+GsumnTC(:,j);
        end
    end
    Gx(:,1) = (stepsize/3)*(f_0+4.*f_odd+2.*f_even+f_n);

    xstepsize = 1e-9;
    jsc=sum(Gx(:,1))*(xstepsize)*q;
end



%------------------- Helper Functions ------------------------------------
% Function I_mat
% This function calculates the transfer matrix, I, for reflection and
% transmission at an interface between materials with complex dielectric 
% constant n1 and n2.
function I = I_mat(n1,n2)
r=(n1-n2)/(n1+n2);
t=2*n1/(n1+n2);
I=[1 r; r 1]/t;
end

% Function L_mat
% This function calculates the propagation matrix, L, through a material of
% complex dielectric constant n and thickness d for the wavelength lambda.
function L = L_mat(n,d,lambda)
xi=2*pi*n/lambda;
L=[exp(-1i*xi*d) 0; 0 exp(1i*xi*d)];
end

function I_til = I_til(n1,n2)
r12=(n1-n2)/(n1+n2);
t12=2*n1/(n1+n2);
R12 = abs(r12)^2;
T12 = abs(t12)^2;
r21=(n2-n1)/(n1+n2);
t21=2*n2/(n1+n2);
R21=abs(r21)^2;
T21=abs(t21)^2;
I_til = [1 -R21; R12 T21*T12-R21*R12]/T12;
end

function L_til = L_til(n,d,lambda)
xi=4*pi*imag(n)/lambda;
L_til=[exp(xi*d) 0; 0 exp(-xi*d)];
end
% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function ntotal = LoadRefrIndex(name,wavelengths)

%Data in IndRefr, Column names in IndRefr_names
load TC_and_BC_nk_data.mat

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr(:,strmatch('Wavelength',IndRefr_names));
n=IndRefr(:,strmatch(strcat(name,'_n'),IndRefr_names));
k=IndRefr(:,strmatch(strcat(name,'_k'),IndRefr_names));

% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
ntotal = n_interp+1i*k_interp; 
end