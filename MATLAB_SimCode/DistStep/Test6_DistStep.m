%________________________________________________
%Simulation Code for Staggered Double Vane       |
%Traveling Wave Tube (No Cylindrical Beam)       |
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
ND = 06;            %Number of Disc(or electrons)|
NZ = 16;            %Number of Steps per Lambda_e|
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
% Calculating Impedance - This is for Helix and needs to be updated for
% staggered later
k = 2*pi*F/c;
Gamma = sqrt((Beta^2) - (k^2));
% N is SlowDown Factor
N = Beta/k;
a = 2*A;
num = (((besseli(0,(Gamma*a)))^2)-((besseli(1,(Gamma*a)))^2))*N;
den = ((2*(besseli(0,(Gamma*a)))*(besseli(1,(Gamma*a))))-((besseli(0,(Gamma*a)))^2)+((besseli(1,(Gamma*a)))^2));
Z = (120/(Gamma*a))*((Gamma^4)/(Beta^4))*(num/den);
% Calculating Electric Field
E0 = Beta*sqrt(2*Z*Pin);
% Initializing Conditions at plane, k=1
dz = u0/(F*NZ);
NS = round(Ln/dz);
B_T_Height = 1.1*B_T_Height;

% Electric Field Variations through the length of the tube
Z = 0:Ln/NS:Ln;
for k = 1:NS+1
    Ec(k) = (E0*sin(Beta*Z(k)));
end

figure;
plot(Z,Ec);
title('Circuit Electric Field Representation');
ylabel('Electric Field, E_c -->');
xlabel('Distance --> ');
figure;

for j =1:ND
    t(1,j) = (j-1)*Lambda/(ND*u0);
    ui(1,j) = u0 - u0*0.1*sin(2*pi*F*t(1,j));
    uin(1,j) = 1/ui(1,j);
end
z(1) = 0;
Acc(1:NS) = 0;
%Calculating Conditions for other planes
for k=1:NS-1
    for j=1:ND
        Acc(k+1) = (e*E0*sin(2*pi*F*t(k,2))/ms)*((1-(u0/c)^2)^(3/2))*((uin(k,j))^3);
        uin(k+1,j) = uin(k,j) + dz*Acc(k);
        t(k+1,j) = t(k,j) + dz*uin(k,j);
    end
    z(k+1) = z(k) + dz;
end
%Plotting Electron Bunching
for j = 1:ND
    plot(t(:,j),z);
    hold on
end
xlabel('Time -->');
ylabel('Distance -->');
title('Electron Bunching Representation at each wavelength cycle');
figure;
for k=1:NS
    plot(uin(k,:));
    hold on
end
title('Inverted Velocity for each Electron');
xlabel('Electron Number');
ylabel('Inverted Velocity');
current_constant = 2*F*q*BeamT;
ang_f = 2*pi*F;
for j=1:NS
    %val(j,:) = (sin(j*ang_f*dz/(2*u0))/(j*ang_f*dz/(2*u0)))*exp(-1i*ang_f*t(j,:));
    val(j,:) = exp(-1i*ang_f*t(j,:));
end
for j=1:NS
    irf(j) = 100000*current_constant*abs(sum(val(j,:)));
end
figure;
plot(z,irf/I0);
title('Current Representation at Different Planes due to Electron Bunching');
xlabel('Distance -->');
ylabel('I_r_f/I_o -->');

%Calculating Space Charge Electric Field (E_sc) for rectangular sheet beam

X = 0:Lambda/2000:Lambda;
E_rough = 0;
E_const = 2*I0*c*c/(pi*W*B_T_Height*BeamT*BeamW);
for b = 1:2001
    for l = 1:10
        for m = 1:10
            E_rough = E_rough + (2*B_T_Height/(l*pi))*(2*W/(m*pi))*sin(m*pi*BeamW/(2*W))*((sin(m*pi/2))^2)*exp(-1*pi*m*X(b)/W)*sin(l*pi*BeamT/(2*B_T_Height)*((sin(l*pi/2))^2))*exp(-1*pi*l*X(b)/B_T_Height)*X(b);
        end
    end
    E_sc(b) = E_const*E_rough;
    E_rough = 0;
end
E_sc_total = sum(E_sc);
figure;
plot(X,E_sc);
title('Space Charge Electric Field for Staggered TWT');
xlabel('Wavelength -->');
ylabel('Space  Charge Electric Field -->');

% Calculating Force due to Space Charge Electric Field
Force = q*E_sc;
figure;
plot(X,Force);
title('Force Representation for Rectangular Sheet Beam');
xlabel('Wavelength -->');
ylabel('Force -->');

% Beam Tracking with Space Charge Field

SCF(1:NS,1:ND) = 0;
E_rough = 0;
E_const = 2*I0*c*c/(pi*W*B_T_Height*BeamT*BeamW);

for i=1:ND
    pos(1,i) = i*dz;
    t(1,i) = (i-1)*Lambda/(ND*u0);
    ui(1,i) = u0 - u0*0.1*sin(2*pi*F*t(1,i));
    uine(1,i) = 1/u0;
end

z(1) = 0;
Acc_E(1,:) = 0;

for k = 1:NS-1
    for i = 1:ND
        for j = 1:ND
            % Calculating Relative Position of each electron in different
            % planes
            relpos(k,i,j) = pos(k,i) - pos(k,j);
            % Calculating Space Charge Electric Field applied by each
            % electron on the electron in study for different planes
            for l = 1:10
                for m = 1:10
                    E_rough = E_rough + (2*B_T_Height/(l*pi))*(2*W/(m*pi))*sin(m*pi*BeamW/(2*W))*((sin(m*pi/2))^2)*exp(-1*pi*m*sign(relpos(k,i,j))*relpos(k,i,j)/W)*sin(l*pi*BeamT/(2*B_T_Height)*((sin(l*pi/2))^2))*exp(-1*pi*l*sign(relpos(k,i,j))*relpos(k,i,j)/B_T_Height)*sign(relpos(k,i,j))*relpos(k,i,j);
                end
            end
         SCField(k,i,j) = E_const*E_rough;
         E_rough = 0;
        end
        % Calculating Total Space Charge Field applied on the electron in
        % study by all neighbouring electrons
        SCF(k,i) = sum(SCField(k,i,:));
        Acc_E(k+1,i) = (e*(Ec(k)+SCF(k,i))/ms)*((1-(u0/c)^2)^(3/2))*((uine(k,i))^3)/10^12;
        %Acc_E(k+1,i) = 0;
        uine(k+1,i) = uine(k,i) + dz*Acc_E(k,i);
        t(k+1,i) = (t(k,i) + dz*uine(k,i));
        pos(k+1,i) = t(k+1,i)/uine(k+1,i);
    end
    z(k+1) = z(k) + dz;
end

figure;
for i = 1:ND
    plot(t(:,i),z);
    hold on
end
title('Electron Tracking with Space Charge Field Applied');
xlabel('Time -->');
ylabel('Distance -->');
current_constant = 2*F*q*BeamT;
ang_f = 2*pi*F;
for j=1:NS
    val(j,:) = exp(-1i*ang_f*t(j,:));
end
for j=1:NS
    irf(j) = 100000*current_constant*abs(sum(val(j,:)));
end
figure;
plot(z,irf/I0);
title('Current Representation at Different Planes due to Electron Bunching with Space Charge Field Applied');
xlabel('Distance -->');
ylabel('I_r_f/I_o -->');