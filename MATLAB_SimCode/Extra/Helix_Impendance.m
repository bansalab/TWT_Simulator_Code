clc
close all;

c = 299792458;

f = input('Enter Frequnency in THz            : ');
p = input('Enter Pitch of the Waveguide in um : ');
a = input('Enter Beam Tunnel Height in um     : ');

f = f*1000000000000;
p = p/1000000;
a = a/1000000;

%Calculating Beta
Beta = 2.5*pi/p;

%Calculating K
k = 2*pi*f/c;

%Calculating Gamma
Gamma = sqrt((Beta^2) - (k^2));

%Calculating SlowDown Factor N
N = Beta/k;

%Calculating Impedance
num = (((besseli(0,(Gamma*a)))^2)-((besseli(1,(Gamma*a)))^2))*N;
den = ((2*(besseli(0,(Gamma*a)))*(besseli(1,(Gamma*a))))-((besseli(0,(Gamma*a)))^2)+((besseli(1,(Gamma*a)))^2));
Imp = (120/(Gamma*a))*((Gamma^4)/(Beta^4))*(num/den);

disp(' ');
disp(' ');
disp('Impedance');
disp(Imp);