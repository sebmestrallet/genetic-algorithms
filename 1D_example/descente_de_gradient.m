clear all
close all

%pour obtenir la dérivée
%syms f(x)
%f(x) = sin(2*pi*0.01*x).*(-x)*0.02 -4;
%Df = diff(f);

lambda = 40;%10

%fonctionnelle avec minimums locaux
f = @(x) sin(2*pi*0.01*x).*(-x)*0.02 -4;
Df = @(x) - sin((pi*x)/50)/50 - (x*pi*cos((pi*x)/50))/2500;

%fonctionnelle parabolique
% f = @(x) 0.001.*x.*x - 0.255.*x;
% Df = @(x) x/500 - 51/200;

x_vecteur = 0:1:2^8;
y_vecteur = f(x_vecteur);

step = 0;
max_steps = 200;
score_vector = NaN(1,max_steps+1);%+1 case pour l'interation 0

x = randi(2^8);%individu au hasard
y = f(x);

figure
plot(x_vecteur,y_vecteur);
hold on
plot(x,y,'.k');
hold off

y_precedent = inf;

while abs(y-y_precedent) > 1e-4
    step = step + 1;
    y_precedent = y;
    
    x = x - lambda*Df(x);
    y = f(x);
    score_vector(1,step+1) = -y;
    
    fprintf("\nPas %i --------------\n",step);
    fprintf("x = %f\n",x);
    fprintf("y = %f\n",y);
    
    plot(x_vecteur,y_vecteur);
    hold on
    plot(x,y,'.k',"MarkerSize",10);
    hold off
    xlim([0 2^8]);
    xticks([0 64 128 192 256]);
    drawnow limitrate
    
    pause(0.05);%attente de 50ms
end

fprintf("\nCritère de convergence atteint\n");

%tracé de l'évolution du score
figure
score_vector = score_vector(1,1:step+1);
plot(0:step,score_vector);
xlabel("Pas");
ylabel("Score");
title("Evolution du score");
xlim([0 step])



