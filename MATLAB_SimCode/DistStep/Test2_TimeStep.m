%_________________________________________________________________________
%Simulation Code for Staggered Double Vane                                |
%Traveling Wave Tube                                                      |
%_________________________________________________________________________|
%Author: Aakash Bansal                                                    |  
%Dated : 13th July 2017                                                   |
%Place : MWT Division, CSIR-CEERI, Pilani                                 |
%Guide : Dr. Vishnu Srivastava, Emeritus Sci.                             |
%_________________________________________________________________________|

close all
clear all
clc

%_________________________________________________________________________
%Constant Values                                                          |
%_________________________________________________________________________|
c = 299792458;      %Speed of Light                                       |
e = 1.6*10^(-19);   %Electron Charge                                      |
m = 9.11*10^(-31);  %Mass of Electron                                     |
%_________________________________________________________________________|

%_________________________________________________________________________
%Variable Values                                                          |
%_________________________________________________________________________|
F = 220000000000;               %Operating Frequency                      |
I = 100;                        %Beam Current in mA                       |
V = 20;                         %Beam Voltage in kV                       |
BW= 50660000000;                %BandWidth of the RF Signal               |
Beam_T = 100;                   %Beam Thickness                           |
Beam_W = 600;                   %Beam Width                               |
ND = 24;                        %Number of Disc(or electrons)             |
NS = 16;                        %Number of Steps per Lambda_rf            |
NR = 1;                         %Number of cycles                         |
Pin= 1;                         %Input Power of RF Signal                 |
%_________________________________________________________________________|

%Calculating Beam Velocity
%Using Equation 1 specified in the documentation
u0 = c*sqrt(1-1/(1+(V/511))^2);
%Calculating Axial Propagation Constant 
%Using equation 2 specified in the documentation
Beta = 2*pi*F/u0;
%Calculating Pitch
Pitch = 2.5*pi/Beta;
%Calculating Gap
Gap = 0.3*Pitch;
%Calculating Vane Thickness
V_Thick = Pitch - Gap;
%Calculating Beam Tunnel Height
A = 1/Beta;
B_T_Height = 2*A;
%Calculating Vane Height
V_Height = 3 * B_T_Height;
%Calculating Total Height
H = (2*V_Height) + B_T_Height;
%Calculating upper and lower cutoff frequency
Fu = F + BW/2;
Fl = F - BW/2;
%Calculating Width of the Rectangular Waveguide
W = c/(2*Fl);
%Calculating Current Density
J = (I*10^5)/(Beam_T*Beam_W);
%Calculating Wavelength
Lambda = 2*pi/Beta;
%Calculating Charge
q = (J*Beam_T*Beam_W)/(1000000*ND*6*F);
%Calculating Impedance - This is for Helix and needs to be updated for
%staggered later
k = 2*pi*F/c;
Gamma = sqrt((Beta^2) - (k^2));
%N is SlowDown Factor
N = Beta/k;
a = 2*A;
num = (((besseli(0,(Gamma*a)))^2)-((besseli(1,(Gamma*a)))^2))*N;
den = ((2*(besseli(0,(Gamma*a)))*(besseli(1,(Gamma*a))))-((besseli(0,(Gamma*a)))^2)+((besseli(1,(Gamma*a)))^2));
Z = (120/(Gamma*a))*((Gamma^4)/(Beta^4))*(num/den);
%Calculating Electric Field
E0 = Beta*sqrt(2*Z*Pin);
dt = NR/(F*ND);
t(1,1) = 0;
z(1:ND,1:NS) = 0;

%Calculating Velocity for all electrons

v(1) = u0;

for k = 1:ND-1
    for j = 1:NS-1
        z(k,j+1) = z(k,j) + dt*v(k);
        t(k,j+1) = t(k,j) + dt;
    end
    t(k+1,1) = t(k,2);
    v(k+1) = v(k) + dt*e*E0*1000*sin(2*pi*F*t(k,2))*sqrt(1-(u0/c)^2)/m;
end

for k = 1:ND
    plot(t(k,:),z(k,:));
    hold on
end
xlabel('Time -->');
ylabel('Distance -->');
title('Electron Bunching Representation at each wavelength cycle');

figure;
plot(v);
xlabel('Electron Count -->');
ylabel('Velocity -->');
title('Variation in Velocity with electron number (or Disc Number)');

%RF Beam Current Calculations
current_constant = 2*F*e*Beam_T/1000000;
ang_f = 2*pi*F;
for j=1:ND-1
    val(:,j) = (sin(j*ang_f*Beam_T/(2*v(j)*1000000))/(j*ang_f*Beam_T/(2*v(j)*1000000)))*exp(-i*j*ang_f*t(j,:));
end
i(1)=0;
for j=1:ND-1
    i(j+1) = current_constant*sum(val(:,j));
end
figure;
abs_i = abs(i);
plot(z(:,NS),abs_i);
title('Current Representation at Different Planes due to Electron Bunching');
xlabel('Distance -->');
ylabel('Current -->');