%________________________________________________
%Simulation Code for Staggered Double Vane       |
%Traveling Wave Tube                             |
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
m = 9.11*10^(-31);  %Mass of Electron            |
%________________________________________________|

%________________________________________________
%Variable Values                                 |
%________________________________________________|
F = 220000000000;   %Operating Frequency         |
I = 100;            %Beam Current in mA          |
V = 20;             %Beam Voltage in kV          |
BW= 50660000000;    %BandWidth of the RF Signal  |
Beam_T = 100;       %Beam Thickness              |
Beam_W = 600;       %Beam Width                  |
ND = 24;            %Number of Disc(or electrons)|
NZ = 16;            %Number of Steps per Lambda_e|
Ln = 0.01;          %Length of Tube in m         |
Pin= 1;             %Input Power of RF Signal    |
%________________________________________________|

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
dz = 2*c/(F*NZ);
NS = round(Ln/dz);
t(1,1) = 0;
z(1:NS,1:ND) = 0;

%Calculating Velocity for all electrons
Acc(1) = 0;
uin(1) = 1/u0;
for k = 1:NS
    for j = 1:ND
        z(k,j+1) = z(k,j) + dz;
        uin(k+1) = uin(k) + dz*Acc(k);  
        t(k,j+1) = t(k,j) + dz*uin(k);
        t(k+1,1) = t(k,2) + 0.5*Acc(k)*dz*dz;
    end
    Acc(k+1) = (e*E0*sin(2*pi*F*t(k,2))*210/m)*((1-(u0/c)^2)^(3/2))*((uin(k))^3);
end

for k = 1:NS
    plot(t(k,:),z(k,:));
    hold on
end
xlabel('Time -->');
ylabel('Distance -->');
title('Electron Bunching Representation at each wavelength cycle');
figure;
plot(uin);
title('Inverted Velocity for each Electron');
xlabel('Electron Number');
ylabel('Inverted Velocity');
current_constant = 2*F*e*Beam_T/1000000;
ang_f = 2*pi*F;
for j=1:NS-1
    val(:,j) = (sin(j*ang_f*Beam_T*uin(j)/(2*1000000))/(j*ang_f*Beam_T*uin(j)/(2*1000000)))*exp(-i*j*ang_f*t(:,j));
end
for j=1:NS-1
    i(j) = abs(current_constant*sum(val(:,j)));
end
figure;
plot(i);