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
ND = 24;            %Number of Disc(or electrons)|
NZ = 16;            %Number of Steps per Lambda_e|
Ln = 0.004;         %Length of Tube in m         |
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
%Initializing Conditions at plane, k=1
dz = c/(F*NZ);
NS = round(Ln/dz);

for j =1:ND
    t(1,j) = (j-1)*Lambda/(ND*u0);
    ui(j) = u0 - u0*0.1*sin(2*pi*F*t(1,j));
    uin(j) = 1/ui(j);
end
z(1) = 0;
Acc(1:NS) = 0;
%Calculating Conditions for other planes
for k=1:NS-1
    for j=1:ND
        %Acc(k+1) = (e*E0*sin(2*pi*F*t(k+1,2))/m)*((1-(u0/c)^2)^(3/2))*((uin(k))^3);
        %uin(k+1,j) = uin(k,j) + 0.1*sin();
        t(k+1,j) = t(k,j) + dz*uin(j);
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
plot(uin);
title('Inverted Velocity for each Electron');
xlabel('Electron Number');
ylabel('Inverted Velocity');
current_constant = 2*F*q*Beam_T;
ang_f = 2*pi*F;
for j=1:NS
    val(j,:) = exp(-1i*ang_f*t(j,:));
end
for j=1:NS
    irf(j) = current_constant*abs(sum(val(j,:)))/I;
end
figure;
plot(z,irf);
title('Current Representation at Different Planes due to Electron Bunching');
xlabel('Distance -->');
ylabel('Current -->');

%Calculating Space Charge Electric Field (E_sc) for rectangular sheet beam

X = 0:Lambda/2000:Lambda/2;
E_rough = 0;
E_const = 2*I*c*c*(10^(-12))/(pi*W*B_T_Height*Beam_T*Beam_W);
for b = 1:1001
    for l = 1:20
        for m = 1:20
            E_rough = E_rough + (2*B_T_Height/(l*pi))*(2*W/(m*pi))*sin(m*pi*Beam_W/(2*W*1000000))*((sin(m*pi/2))^2)*exp(-1*pi*m*X(b)/W)*sin(l*pi*Beam_T/(2*B_T_Height*1000000)*((sin(l*pi/2))^2))*exp(-1*pi*l*X(b)/B_T_Height)*X(b);
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

%Calculating Force due to Space Charge Electric Field
Force = q*E_sc;
figure;
plot(X,Force);
title('Force Representation for Rectangular Sheet Beam');
xlabel('Wavelength -->');
ylabel('Force -->');

%Force Calculations using VS Formula
vs_rough = 0;
vs_a = sqrt(W*B_T_Height/pi);
vs_b = sqrt(Beam_T*Beam_W/pi);
mu = [3.831706 7.015586 10.173468 13.323692 16.470630]./vs_a;
vs_const = (4*q*q/(pi*ep0*vs_a*vs_a*vs_b*vs_b*dz*dz));
for b = 1:1001
    for m = 1:5
        if(X(b)<dz)
            vs_rough = vs_rough + ((besselj(1,(mu(m)*vs_b))^2)/((mu(m)^4)*(besselj(1,(mu(m)*vs_a))^2)))*(1-exp(-1*mu(m)*sign(X(b))*X(b))-(exp(-mu(m)*dz)*sinh(mu(m)*sign(X(b))*X(b))));
        else
            %vs_rough = vs_rough + exp(-mu(m)*X(b));
            vs_rough = vs_rough + ((besselj(1,(mu(m)*vs_b))^2)/((mu(m)^4)*(besselj(1,(mu(m)*vs_a))^2)))*(2*((sinh(mu(m)*dz/2))^2)*exp(-mu(m)*X(b)));
        end
    end
    vs_Force(b) = vs_const*vs_rough;
    vs_rough = 0;
end
figure;
plot(X,vs_Force);
title('Force Representation for Cylindrical Beam');
xlabel('Distance -->');
ylabel('Force -->');

CD_Rec = I/(pi*vs_b*vs_b*u0);
CD_Cyl = I/(Beam_T*Beam_W*u0);