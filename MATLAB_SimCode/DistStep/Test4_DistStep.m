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
NZ = 08;            %Number of Steps per Lambda_e|
Ln = 0.12;         %Length of Tube in m         |
Pin= 0;             %Input Power of RF Signal    |
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
        t(k+1,j) = t(k,j) + dz*uin(k,j);
        Acc(k+1) = (e*E0*sin(2*pi*F*t(k+1,2))/m)*((1-(u0/c)^2)^(3/2))*((uin(k,j))^3);
        uin(k+1,j) = uin(k,j) + dz*Acc(k);
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
% figure;
% plot(X,E_sc);
% title('Space Charge Electric Field for Staggered TWT');
% xlabel('Wavelength -->');
% ylabel('Space  Charge Electric Field -->');

%Calculating Force due to Space Charge Electric Field
Force = q*E_sc;
figure;
plot(X,Force);
title('Force Representation for Rectangular Sheet Beam');
xlabel('Wavelength -->');
ylabel('Force -->');

%Force Calculations using VS Formula
vs_rough = 0;
%vs_a = sqrt(W*B_T_Height/pi);
vs_a = B_T_Height;
vs_b = BeamT;
mu = [3.831706 7.015586 10.173468 13.323692 16.470630]./vs_a;
vs_const = (4*q*q/(pi*ep0*vs_a*vs_a*vs_b*vs_b*dz*dz));
for b = 1:2001
    for m = 1:1
        if(X(b)<=dz)
            vs_rough = vs_rough + ((besselj(1,(mu(m)*vs_b))^2)/((mu(m)^4)*(besselj(1,(mu(m)*vs_a))^2)))*(1-exp(-1*mu(m)*sign(X(b))*X(b))-(exp(-1*mu(m)*dz)*sinh(mu(m)*sign(X(b))*X(b))));
        else
            %vs_rough = vs_rough + exp(-mu(m)*X(b));
            vs_rough = vs_rough + ((besselj(1,(mu(m)*vs_b))^2)/((mu(m)^4)*(besselj(1,(mu(m)*vs_a))^2)))*(2*((sinh(mu(m)*dz/2))^2)*exp(-1*mu(m)*X(b)));
        end
    end
    vs_Force(b) = vs_const*vs_rough;
    vs_rough = 0;
end
figure;
plot(X,vs_Force);
hold on
plot(X,Force,'--b');
%title('Force Representation for Cylindrical and Rectangular Beam');
xlabel('Distance (m) -->');
ylabel('Force (N) -->');

CD_Rec = I0/(pi*vs_b*vs_b*u0);
CD_Cyl = I0/(BeamT*BeamW*u0);

%DC Electric Field according to P. C. Panda
d = B_T_Height;

for j=1:10
    dc_k(j) = (2*(j-1)+1)*pi/d;
end
X = 0:Lambda/2000:Lambda/2;
Ex_const = 2*CD_Rec/(d*ep0);
Ex_rough = 0;

for b = 1:1001
    for j = 1:10
        Ex_rough = Ex_rough + sin(dc_k(j)*BeamT/2)*exp(-dc_k(j)*W/2)*cos(dc_k(j)*0)*sinh(dc_k(j)*X(b));
    end
    Ex(b) = Ex_const*Ex_rough;
    Ex_rough = 0;
end
figure;
plot(X,Ex);
title('X-Component of Space Charge Field as explained by Panda');
ylabel('Ex (in V/m^2) -->');
xlabel('Distance -->');
hold on
Ey_rough = 0;
Ey_const = 4*CD_Rec/(d*ep0);

Y = 0:B_T_Height/2000:B_T_Height/2;

for b = 1:1001
    for j=1:10
        Ey_rough = Ey_rough + sin(dc_k(j)*BeamT/2)*(sin(dc_k(j)*Y(b))/(dc_k(j)^2));
    end
    Ey(b) = Ey_const*Ey_rough;
    Ey_rough = 0;
end
figure;
plot(Y,Ey);
title('Y-Component of Space Charge Field as explained by Panda');
ylabel('E_y (in V/m^2) -->');
xlabel('Distance -->');
hold on
vs_r = 0:vs_a/2000:vs_a/2;
for b = 1:1001
    Ey_Cyl(b) = J*vs_r(b)/(2*u0*ep0);
end
% figure;
plot(vs_r,Ey_Cyl);
% title('R-Component of Space Charge Field as explained by Panda');
% ylabel('E_r (in V/m^2) -->');
% xlabel('Distance -->');