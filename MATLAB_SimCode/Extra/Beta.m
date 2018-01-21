%________________________________________________
%Simulation Code for Staggered Double Vane       |
%Backward Wave Amplifier (No Cylindrical Beam)   |
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
C = 299792458;      %Speed of Light              |
e = 1.6*10^(-19);   %Electron Charge             |
ms = 9.11*10^(-31);  %Mass of Electron           |
ep0=8.85*10^(-12);  %Permittivity of Free Space  |
%________________________________________________|

%________________________________________________ 
%Variable Values                                 |
%________________________________________________|
F = 170000000000;   %Operating Frequency         |
I = 100;            %Beam Current in mA          |
V = 20;             %Beam Voltage in kV          |
BW= 50660000000;    %BandWidth of the RF Signal  |
BEAM_T = 100;       %Beam Thickness in um        |
BEAM_W = 600;       %Beam Width in um            |
ND = 12;            %Number of Disc(or electrons)|
NZ = 8;             %Number of Steps per Lambda_e|
Ln = 0.10;          %Length of Tube in m         |
Pin= 0.001;         %Input Power of RF Signal    |
%________________________________________________|
BEAM_T = BEAM_T/1000000;
BEAM_W = BEAM_W/1000000;
I0 = I/1000;
% Calculating Beam Velocity
% Using Equation 1 specified in the documentation
VP = C*sqrt(1-1/(1+(V/511))^2);
% Calculating Axial Propagation Constant 
% Using equation 2 specified in the documentation
BETA = 2*pi*F/VP;
% Calculating Pitch
PITCH = 2*pi/BETA;
% Calculating Gap
GAP = 0.3*PITCH;
% Calculating Vane Thickness
V_THICK = PITCH - GAP;
% Calculating Beam Tunnel Height
A = 1/BETA;
B_T_HEIGHT = 2*A;
% Calculating Vane Height
V_HEIGHT = 3 * B_T_HEIGHT;
% Calculating Total Height
H = (2*V_HEIGHT) + B_T_HEIGHT;
% Calculating upper and lower cutoff frequency
Fu = F + BW/2;
Fl = F - BW/2;
% Calculating Width of the Rectangular Waveguide
W = C/(2*Fl);
% Calculating Current Density
J = (I0)/(BEAM_T*BEAM_W);
% Calculating Wavelength
Lambda = 2*pi/BETA;
% Calculating Charge
q = (J*BEAM_T*BEAM_W)/(ND*6*F);

KX = pi/BEAM_W;
K0 = 2*F*pi/C;

GAMMA_N = sqrt((KX^2)+(BETA^2)-(K0^2));

BETA_N = 2.5*pi/PITCH;

LM = sqrt((KX^2)+((pi/GAP)^2)-(K0^2));

M = ((PITCH*GAP*GAMMA_N*cosh(GAMMA_N*A))-sinh(GAMMA_N*A)*(LM*tanh(LM*V_HEIGHT))*GAP*GAP);
N = ((PITCH*GAP*GAMMA_N*sinh(GAMMA_N*A))-cosh(GAMMA_N*A)*(LM*tanh(LM*V_HEIGHT))*GAP*GAP);
P = ((PITCH*GAP*GAMMA_N*cosh(GAMMA_N*A))-sinh(GAMMA_N*A)*(LM*tanh(LM*V_HEIGHT))*GAP*GAP)*exp(-1j*(BETA_N-BETA)*PITCH/2);
Q = ((PITCH*GAP*GAMMA_N*sinh(GAMMA_N*A))-cosh(GAMMA_N*A)*(LM*tanh(LM*V_HEIGHT))*GAP*GAP)*exp(-1j*(BETA_N-BETA)*PITCH/2);

O = (M*Q)-(N*P);