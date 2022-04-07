L = 100;%largeur du rectangle
H = 50;%hauteur du rectangle
nb_points_par_snake = 40;
delta_theta = (2*pi)/nb_points_par_snake;%pas en radians
thetas = 0:delta_theta:2*pi-delta_theta;

figure
polarplot([],[]);
hold on
for i = 1:nb_points_par_snake
    theta = thetas(i);
    polarplot(theta,rho_max(L,H,theta),'.k',"MarkerSize",10);
end
hold off
title("Calcul du rayon maximal rho_{max} selon l'angle theta (image rectangulaire)")
thetaticks(thetas*360/(2*pi));

function valeur = rho_max(L,H,theta)

    a = L/2;
    b = H/2;
    if abs(tan(theta)) <= b/a
        valeur = a/abs(cos(theta));
    else
        valeur = b/abs(sin(theta));
    end
end