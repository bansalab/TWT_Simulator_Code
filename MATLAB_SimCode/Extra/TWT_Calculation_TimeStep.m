close all
clear all
clc

c = 299792458;
e = 1.6*10^(-19);
m = 9.11*10^(-31);

F = input('Center Frequency (in GHz): ');
I = input('Beam Current (in mA)     : ');
V = input('Beam Voltage (in kV)     : ');
BW= input('Bandwidth (in GHz)       : ');

F = F*1000000000;
BW= BW*1000000000;

%Calculating Beam Velocity
v0 = c*sqrt(1-1/(1+(V/511))^2);
%Calculating Axial Propagation Constant
Beta = 2*pi*F/v0;
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

%Presenting All the Values
disp('All Values are in meter');
disp(' ');
disp(' ');
disp('Pitch');
disp(Pitch);
disp('Gap');
disp(Gap);
disp('Vane Thickness');
disp(V_Thick);
disp('Beam Tunnel Height');
disp(B_T_Height);
disp('Vane Height');
disp(V_Height);
disp('Total Height');
disp(H);
disp('Width');
disp(W);
disp(' ');
disp('All Values are in Hz');
disp(' ');
disp(' ');
disp('Lower Cut off Frequency');
disp(Fl);
disp('Upper Cut off Frequency');
disp(Fu);
disp('Central Frequency');
disp(F);

%Sheet Beam Calculations
disp(' ');
disp(' ');
disp('Sheet Beam Calculations');
disp(' ');

%Data from User
Beam_T = input('Thickness of the beam in um     : ');
Beam_W = input('Width of the beam in um         : ');
ND = input('Number of Discs per Wavelength  : ');
NR = input('Number of steps per wavelength  : ');
Ln = input('Length of the tube in mm        : ');
Ln = Ln/1000;
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
%Considering Power of input RF Signal being 1W
%Calculating Electric Field
Pin = 1;
E0 = Beta*sqrt(2*Z*Pin);
%Calculating Interchanged Energy of the Macroparticle
h = Lambda/NR;
NZ = Ln/h;
%zn = Lambda/NZ;
%dE = q*E0*zn;
%Calculating Velocity of each macroparticle in each section
tn = 1/(F*ND);
ang_f = 2*pi*F;
Lambdae = v0/F;
%Time Step Calculations
Acc = (e*E0*sin(ang_f*tn)/(m*10))*((1-(v0/c)^2)^(3/2))*((1/v0)^3);
% A will be dt/dz
A(1:NZ,1:ND) = 0;
t(1:NZ,1:ND) = 0;
z(1:NZ,1:ND) = 0;
A(1) = 1/v0;
z(1) = 0;
t(1) = 0;
for i=1:ND
    t(1,i) = (i-1)*Lambdae/(ND*v0);
end
A(1,:) = 1/v0;
for k=1:NZ
    for i=1:ND
        A(k+1,i) = A(k,i) + Acc*h;
        t(k+1,i) = t(k,i) + (A(k,i)+(0.5*Acc*h))*h;
        z(k+1,i) = z(k,i) + h;
    end
end
for k=1:ND
    plot(t(:,k),z(:,k));
    hold on
end
xlabel('Time');
ylabel('z');