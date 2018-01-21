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
ND = input('Number of Time Distributions    : ');
NR = input('Number of Traverse Positions    : ');
NZ = input('Number of Distributions in Z    : ');

%Calculating Current Density
J = (I*10^5)/(Beam_T*Beam_W);
%Calculating Wavelength
Lambda = v0/F;
%Calculating Charge
q = (J*Beam_T*Beam_W)/(1000000*ND*NR*F);
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
zn = Lambda/NZ;
dE = q*E0*zn;
%Calculating Velocity of each macroparticle in each section
tn = 1/(F*ND);
ang_f = 2*pi*F;
%Euler Method Implementation to Calculate Velocity and distance covered by 
%each electrons after each section
a=0;
b=1/F;
h = (b-a)*1000/ND;
vn(1:ND) = 0;
t(1:ND,1:NZ) = 0;
z(1:ND,1:NZ) = 0;
t(1,1) = 0;
vn(1) = v0;
z(1,1) = 0;
for k=1:1:(ND)
for i=1:1:(NZ)
    t(k,i+1) = t(k,i) + h;
    if k == 1
    vn(k+1) = vn(k) + h*e*E0*sin(ang_f*t(k,i))/m;
    else
    vn(k+1) = vn(k) + h*e*E0/m;
    end
    z(k,i+1) = z(k,i) + h*vn(k);
    t(k+1,i) = t(k,i) + h;
    z(k+1,2) = z(k,1) + h*vn(k+1);
end
end
for k=1:1:ND
plot(t(k,:),z(k,:));
hold on
end

%Time Step Calculations
Acc = (-e*E0/m)*(1-(v0/c)^2)^(3/2)*(1/v0)^3;

%Display Data
disp(' ');
disp(' ');
disp('Current Density (J) in A/cm2');
disp(J);
disp('Charge of the macroparticle (q)');
disp(q);
disp('Impedance (Z) in ohms');
disp(Z);
disp('Electric Field in Tube (E) in V/cm');
disp(E0/100);
disp('Interchanged Energy of the Macroparticle (dE)');
disp(dE);
disp('Angular Frequency');
disp(ang_f);
disp('Each Time Iteration');
for i=1:ND
    disp(t(i));
end
disp(t(1));
disp('Velocity of each Macroparticle');
% for i=1:ND
%     disp(vel(i));
% end
