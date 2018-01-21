%________________________________________________
%Simulation Code for Look Up Table for           |
%Charge Force in a Traveling Wave Tube           |
%________________________________________________|
%Author: Aakash Bansal                           |
%Dated : 26th July 2017                          |
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
ND = 40;            %Number of Disc(or electrons)|
Ln = 0.005;         %Length of Tube in m         |
Pin= 1;             %Input Power of RF Signal    |
%________________________________________________|
BeamT = Beam_T/1000000;
BeamW = Beam_W/1000000;
I0 = I/1000;
% Calculating Beam Velocity
% Using Equation 1 specified in the documentation
u0 = c*sqrt(1-1/(1+(V/511))^2);
% Calculating Axial Propagation Constant 
% Using equation 2 specified in the documentation
Beta = 2*pi*F/u0;
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
% Initializing Conditions at plane, k=1
dz = u0/(F*ND);
B_T_Height = 1.1*B_T_Height;
% Calculating Rectangular and Cylindrical Space Charge Force
X = 0:Lambda/2000:Lambda;
E_rough = 0;
E_const = 2*I0*c*c/(pi*W*B_T_Height*BeamT*BeamW);
Cyl_rough = 0;
Cyl_b = sqrt(BeamT*BeamW/pi);
Cyl_a = sqrt(W*B_T_Height/pi);
mu = [3.831706 7.015586 10.173468 13.323692 16.470630]./Cyl_a;
Cyl_const = (4*q*q/(pi*ep0*Cyl_a*Cyl_a*Cyl_b*Cyl_b*dz*dz));

for b = 1:2001
    for l = 1:10
        for m = 1:10
            E_rough = E_rough + (2*B_T_Height/(l*pi))*(2*W/(m*pi))*sin(m*pi*BeamW/(2*W))*((sin(m*pi/2))^2)*exp(-1*pi*m*X(b)/W)*sin(l*pi*BeamT/(2*B_T_Height)*((sin(l*pi/2))^2))*exp(-1*pi*l*X(b)/B_T_Height)*X(b);
        end
    end
    
     for m = 1:5
        if(X(b)<dz)
            Cyl_rough = Cyl_rough + ((besselj(1,(mu(m)*Cyl_b))^2)/((mu(m)^4)*(besselj(1,(mu(m)*Cyl_a))^2)))*(1-exp(-1*mu(m)*X(b))-(exp(-1*mu(m)*dz)*sinh(mu(m)*X(b))));
        else
            Cyl_rough = Cyl_rough + ((besselj(1,(mu(m)*Cyl_b))^2)/((mu(m)^4)*(besselj(1,(mu(m)*Cyl_a))^2)))*(2*((sinh(mu(m)*dz/2))^2)*exp(-1*mu(m)*X(b)));
        end
     end
    
    E_sc(b) = E_const*E_rough;
    E_rough = 0;
    Cyl_Force(b) = Cyl_const*Cyl_rough;
    Cyl_rough = 0;
end
Rec_Force = q*E_sc;

plot(X,Rec_Force);
hold on
plot(X,Cyl_Force);
title('Force Representation for Cylindrical and Rectangular Beam');
xlabel('Distance -->');
ylabel('Force -->');
legend('Rectangular SCF','Cylindrical SCF');