%________________________________________________
%Simulation Code for Staggered Double Vane       |
%Backward Wave Amplifier (No SCF)                |
%________________________________________________|
%Author: Aakash Bansal                           |
%Dated : 13th July 2017                          |
%Place : MWT Division, CSIR-CEERI, Pilani        |
%Guide : Dr. Vishnu Srivastava, Emeritus Sci.    |
%________________________________________________|

close all
clear all
clc

%________________________________________________
%Constant Values                                 |
%________________________________________________|
c = 299792458;      %Speed of Light              |
e = 1.6*10^(-19);   %Electron Charge             |
ms = 9.11*10^(-31);  %Mass of Electron           |
ep0=8.85*10^(-12);  %Permittivity of Free Space  |
%________________________________________________|

%________________________________________________ 
%Variable Values                                 |
%________________________________________________|
F = 220000000000;   %Operating Frequency         |
I = 5;              %Beam Current in mA          |
V = 14;             %Beam Voltage in kV          |
BW= 50660000000;    %BandWidth of the RF Signal  |
Beam_T = 100;       %Beam Thickness in um        |
Beam_W = 600;       %Beam Width in um            |
ND = 12;            %Number of Disc(or electrons)|
NZ = 8;             %Number of Steps per Lambda_e|
Ln = 0.12;          %Length of Tube in m         |
Pin= 1000;  %Input Power of RF Signal    |
%________________________________________________|
BeamT = Beam_T/1000000;
BeamW = Beam_W/1000000;
I0 = I/1000;
% Calculating Beam Velocity
% Using Equation 1 specified in the documentation
u0 = c*sqrt(1-1/(1+(V/511))^2);
% Calculating Pitch
Pitch = 0.00014;
% Calculating Axial Propagation Constant 
% Using equation 2 specified in the documentation
Beta = 1.4*pi/Pitch;

vp = 2*pi*F/Beta;

% Calculating Gap
Gap = 0.3*Pitch;
% Calculating Vane Thickness
V_Thick = Pitch - Gap;
% Calculating Beam Tunnel Height
A = 1/Beta;
B_T_Height = 2*A;
% Calculating Vane Height
V_Height = 3 * B_T_Height;
% Calculating Total Height
H = (2*V_Height) + B_T_Height;
% Calculating upper and lower cutoff frequency
Fu = F + BW/2;
Fl = F - BW/2;
% Calculating Width of the Rectangular Waveguide
W = c/(2*Fl);
% Calculating Current Density
J = (I0)/(BeamT*BeamW);
% Calculating Wavelength
Lambda = 2*pi/Beta;
% Calculating Charge
q = (J*BeamT*BeamW)/(ND*6*F);
% Calculating Impedance - This is for Helix and needs to be updated for
% staggered later
k = 2*pi*F/c;
Gamma = sqrt((Beta^2) - (k^2));
% N is SlowDown Factor
N = Beta/k;
a = 2*A;
num = (((besseli(0,(Gamma*a)))^2)-((besseli(1,(Gamma*a)))^2))*N;
den = ((2*(besseli(0,(Gamma*a)))*(besseli(1,(Gamma*a))))-((besseli(0,(Gamma*a)))^2)+((besseli(1,(Gamma*a)))^2));
%Z = (120/(Gamma*a))*((Gamma^4)/(Beta^4))*(num/den);
Z = 0.9675;
% Initializing Conditions at plane, k=1

dz = u0/(F*NZ);
NS = round(Ln/dz);

for x=1:ND
    pos(NS,x) = x*dz;
    t(NS,x) = (x-1)*Lambda/(ND*u0);
    ui(NS,x) = u0 - u0*0.001*sin(2*pi*F*t(NS,x));
    % Modulated Velocity
    uine(NS,x) = 1/ui(NS,x);
    % Unmodulated Velocity
    %uine(1,x) = 1/u0;
end

z(1) = 0;
Acc_E(1,:) = 0;

Vin_N(NS) = sqrt(2*Z*Pin);
Z_Ln = 0:Ln/NS:Ln;
irf(1) = 0;
current_constant = 2*F*q;
ang_f = 2*pi*F;
dVin_N = 0;
curr_const2(1:NS) = 0;
Ec(NS) = Beta*sqrt(2*Z*Pin);

for k = 1:NS-1  
    Ec(NS-k) = Ec(NS-k+1)*exp(1i*Beta*Ln/NS);
end

for k = 1:NS-1
    for x = 1:ND
        Ec(NS-k) = Ec(NS-k)*exp(1i*ang_f*t(NS-k+1,x));
         Acc_E(NS-k,x) = (e*Ec(NS-k)/ms)*((1-(u0/c)^2)^(3/2))*((uine(NS-k+1,x))^3)/10^12;
         uine(NS-k,x) = uine(NS-k+1,x) + dz*Acc_E(NS-k,x);
         t(NS-k,x) = (t(NS-k+1,x) + dz*uine(NS-k+1,x));
         pos(NS-k,x) = t(NS-k,x)/uine(NS-k,x);
         curr_const2(NS-k+1) = curr_const2(NS-k+1) + (sin(3*ang_f*BeamT*0.5/u0)/(3*ang_f*BeamT*0.5/u0))*exp(-1i*3*ang_f*t(NS-k+1,x));
    end
    Vin_N(NS-k) = ((Vin_N(NS-k+1)) + (dVin_N));
    Pout_N(NS-k) = abs(Vin_N(NS-k))*abs(Vin_N(NS-k))/(2*Z);
    Gain(NS-k) = 10*log10(Pout_N(NS-k)/Pin);
    z(k+1) = z(k) + dz;
    if(Pin == 0)
        irf(NS-k+1) = 0;
    else
        irf(NS-k+1) = 10*current_constant*curr_const2(NS-k+1);
    end
    Tou_e = (-1/dz)*log((irf(NS-k+1)/irf(NS-k)));
    dVin_N = (25*Z/2)*abs(irf(NS-k+1))*((1-exp(Beta*dz)/exp(-1*Tou_e*dz))/(1+Beta/Tou_e));
end

figure;
plot(Z_Ln(1:NS-1),real(Ec(1:NS-1)));
title('Circuit Electric Field due to RF Voltage');
ylabel('Electric Field, E_c -->');
xlabel('Distance --> ');

figure;
for x = 1:ND
    plot(t(:,x),z);
    hold on
end
title('Electron Tracking with Space Charge Field and Circuit Field Applied Applied');
xlabel('Time -->');
ylabel('Distance -->');

figure;
plot(z(1:NS-1),5*abs(irf(1:NS-1)/I0));
title('Current Representation at Different Planes due to Electron Bunching with Space Charge Field Applied');
xlabel('Distance -->');
ylabel('I_r_f/I_o -->');

figure;
plot(Z_Ln(1:NS-1),(Vin_N(1:NS-1)));
title('Output Backward Voltage of RF Signal in the tube');
xlabel('Distance --> ');
ylabel('Voltage, V_o_u_t -->');

figure;
plot(Z_Ln(1:NS-1),(Pout_N(1:NS-1)));
title('Output Backward Power of RF Signal in the tube');
xlabel('Distance --> ');
ylabel('Output Power, P_o_u_t (dB)-->');

figure;
plot(Z_Ln(1:NS-1),10*log10(Pout_N(1:NS-1)));
title('Output Backward Power of RF Signal in the tube');
xlabel('Distance --> ');
ylabel('Output Power, P_o_u_t -->');

figure;
plot(Z_Ln(1:NS-1),Gain);