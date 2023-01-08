% figure
% polarplot([],[]);
% nb_theta_divisions = 30;
% thetaticks(0:360/nb_theta_divisions:360);

image_filename = "tile.tif";
image = imread(image_filename);
[H,L] = size(image);


%codage en coordonnées cartésiennes

figure
imagesc(image);
colormap("gray")

hold on
x1 = 40;
y1 = 60;
plot(x1,y1,'.r','MarkerSize',20);
line([0,x1],[y1 y1],'Color','red','LineStyle',':','LineWidth',1);
line([x1,x1],[0 y1],'Color','red','LineStyle',':','LineWidth',1);
x2 = 150;
y2 = 125;
plot(x2,y2,'.r','MarkerSize',20);
line([0,x2],[y2 y2],'Color','red','LineStyle',':','LineWidth',1);
line([x2,x2],[0 y2],'Color','red','LineStyle',':','LineWidth',1);
hold off
title("Codage en coordonnées cartésiennes")
xlabel("x");
ylabel("y");

%problème de la représentation de coordonnées non valides

figure
image_filename = "tile_extended.png";
image = imread(image_filename);
imagesc(image);
hold on
x1 = 40;
y1 = 60;
plot(x1,y1,'.r','MarkerSize',20);
line([0,x1],[y1 y1],'Color','red','LineStyle',':','LineWidth',1);
line([x1,x1],[0 y1],'Color','red','LineStyle',':','LineWidth',1);
x2 = 150;
y2 = 125;
plot(x2,y2,'.r','MarkerSize',20);
line([0,x2],[y2 y2],'Color','red','LineStyle',':','LineWidth',1);
line([x2,x2],[0 y2],'Color','red','LineStyle',':','LineWidth',1);
hold off
title("toutes les représentations ne sont pas des solutions valides")
xlabel("x (puissance de 2 supérieure)");
ylabel("y (puissance de 2 supérieure)");
xticks([0 64 128 192 255])
yticks([0 64 128 192 255])


%codage en coordonnées polaires

nb_points_par_snake = 40;
delta_theta = (2*pi)/nb_points_par_snake;%pas en radians
thetas = 0:delta_theta:2*pi-delta_theta;

figure
image_filename = "tile.tif";
image = imread(image_filename);
imagesc(image);
colormap("gray")
hold on
for i = 1:nb_points_par_snake
    theta = thetas(i);
    
    [x,y] = pol2cart(theta,rho_max(L,H,theta));
    x = x+L/2;
    y = y+H/2;
    line([L/2 x],[H/2 y],'Color','red','LineStyle',':','LineWidth',1);
    
    [point_x,point_y] = pol2cart(theta,50);
    point_x = point_x+L/2;
    point_y = point_y+H/2;
    plot(point_x,point_y,'.r','MarkerSize',20);
end
hold off
title_str = sprintf("Codage en coordonnées polaires, %d points par snake",nb_points_par_snake);
title(title_str);
xlabel("x");
ylabel("y");


function valeur = rho_max(L,H,theta)

    a = L/2;
    b = H/2;
    if abs(tan(theta)) <= b/a
        valeur = a/abs(cos(theta));
    else
        valeur = b/abs(sin(theta));
    end
end