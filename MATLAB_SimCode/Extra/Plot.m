G = [21.56 22 23.1 22.9 22.5];
F = [332 338 342 346 348];
V = [12 13 14 15 16];

Gd = [22 21 23 22.9 22.5];

plot(V,F,'-b');
yyaxis right
plot(V,G,'-o');
hold on
plot(V,Gd,'-*');