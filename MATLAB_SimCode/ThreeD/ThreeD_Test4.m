%________________________________________________
%Simulation Code for Staggered Double Vane       |
%Traveling Wave Tube (No Cylindrical Beam)       |
%Calculating Induced Voltage                     |
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
I = 100;            %Beam Current in mA          |
V = 20;             %Beam Voltage in kV          |
BW= 50660000000;    %BandWidth of the RF Signal  |
Beam_T = 100;       %Beam Thickness in um        |
Beam_W = 600;       %Beam Width in um            |
ND = 12;            %Number of Disc(or electrons)|
NZ = 8;             %Number of Steps per Lambda_e|
Ln = 0.12;          %Length of Tube in m         |
Pin= 0.01;             %Input Power of RF Signal |
%________________________________________________|
BeamT = Beam_T/1000000;
BeamW = Beam_W/1000000;
I0 = I/1000;
% Calculating Beam Velocity
% Using Equation 1 specified in the documentation
vp = c*sqrt(1-1/(1+(V/511))^2);
u0 = vp;
% Calculating Axial Propagation Constant 
% Using equation 2 specified in the documentation
Beta = 2*pi*F/vp;
% Calculating Pitch
Pitch = 2.5*pi/Beta;
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
% Calculating Electric Field
E0 = Beta*sqrt(2*Z*Pin);
% Initializing Conditions at plane, k=1
dz = u0/(F*NZ);
NS = round(Ln/dz);
B_T_Height = B_T_Height;


% Beam Tracking with Space Charge Field

SCF(1:NS,1:ND) = 0;
E_rough = 0;
E_const = 2*I0*c*c/(pi*W*B_T_Height*BeamT*BeamW*10^5);

for x=1:ND
    pos(1,x) = x*dz;
    t(1,x) = (x-1)*Lambda/(ND*u0);
    ui(1,x) = u0 - u0*0.001*sin(2*pi*F*t(1,x));
    % Modulated Velocity
    uine(1,x) = 1/ui(1,x);
    % Unmodulated Velocity
    %uine(1,x) = 1/u0;
end

z(1) = 0;
Acc_E(1,:) = 0;

% Electric Field Variations through the length of the tube

% Calculating Input RF Voltage

Vin(1) = sqrt(2*Z*Pin);

Z_Ln = 0:Ln/NS:Ln;
irf(1) = 0;
current_constant = 2*F*q;
ang_f = 2*pi*F;
dVin = 0;
curr_const2(1:NS) = 0;

ut = sqrt(2*e*V/ms);
ux(1,1:ND) = 0;
uy(1,1:ND) = 0.1*ut;
uz(1,1:ND) = 0.9*ut;
Bz = 1;
posx(1,1:ND) = 0;
posy(1,1:ND) = 0;
posz(1,1:ND) = 0;

for k = 1:NS-1
    for x = 1:ND
        for j = 1:ND
            % Calculating Relative Position of each electron in different
            % planes
            relpos(k,x,j) = posz(k,x) - posz(k,j);
            % Calculating Space Charge Electric Field applied by each
            % electron on the electron in study for different planes
            for l = 1:2
                for m = 1:2
                    E_rough = E_rough + (2*B_T_Height/(l*pi))*(2*W/(m*pi))*sin(m*pi*BeamW/(2*W))*((sin(m*pi/2))^2)*exp(-1*pi*m*sign(relpos(k,x,j))*relpos(k,x,j)/W)*sin(l*pi*BeamT/(2*B_T_Height)*((sin(l*pi/2))^2))*exp(-1*pi*l*sign(relpos(k,x,j))*relpos(k,x,j)/B_T_Height)*relpos(k,x,j);
                end
            end
         SCField(k,x,j) = E_const*E_rough;
         E_rough = 0;
         end
        % Calculating Total Space Charge Field applied on the electron in
        % study by all neighbouring electrons
        SCF(k,x) = sum(SCField(k,x,:));
        Vin(k+1) = (Vin(k)) + sign(dVin)*(dVin);
        Pout(k) = abs(Vin(k))*abs(Vin(k))/(2*Z);
        Gain(k) = 10*log10(Pout(k)/Pin);
        Ec(k) = (1i*Beta*Vin(k)*exp(1i*Beta*Z_Ln(k)));
        Acx(k,x) = (e/u0)*(real(Ec(k)) + (uy(k,x)*Bz))/(ms*uz(k,x));
        Acy(k,x) = (e/u0)*((ux(k,x)*Bz))/(ms*uz(k,x));
        Acz(k,x) = (e/u0)*(real(Ec(k)))/(ms*uz(k,x));
        ux(k,x+1) = ux(k,x) + dz*Acx(k,x);
        uy(k,x+1) = uy(k,x) + dz*Acy(k,x);
        uz(k,x+1) = uz(k,x) + dz*Acz(k,x);
        posx(k,x+1) = posx(k,x) + dz*ux(k,x)/uz(k,x);
        posy(k,x+1) = posy(k,x) + dz*uy(k,x)/uz(k,x);
        posz(k,x+1) = posz(k,x) + dz;
        t(k,x+1) = t(k,x) + posz(k,x)/uz(k,x);
    end
    z(k+1) = z(k) + dz;
    val(k,:) = exp(-1i*ang_f*t(k,:));
    for x = 1:ND
        curr_const2(k) = curr_const2(k) + (sin(3*ang_f*BeamT*0.5/u0)/(3*ang_f*BeamT*0.5/u0))*exp(-1i*3*ang_f*t(k,x));
    end
    %if(Pin == 0)
    %    irf(k) = 0;
    %else
        irf(k+1) = 10*current_constant*curr_const2(k);
    %end
    if(irf(k) == 0)
        dVin = 0;
    else
        Tou_e = (-1/dz)*log((irf(k+1)/irf(k)));
        dVin = (-Z/2)*abs(irf(k))*((1-exp(-1*Beta*dz)/exp(-1*Tou_e*dz))/(1-Beta/Tou_e));
    end
end

scatter3(posx,posy,posz)