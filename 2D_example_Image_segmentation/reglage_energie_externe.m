%la réussite de la segmentation dépend des coefficients de l'énergie
%externe, dont l'influence peut être observer avec ce script

clear variables
close all

emplacement_image = "tile.tif";

sigma = 5;
weline = 0.05;%contribution de l'intensite 
weedge = 1;%contribution des bords 
weterm = 0;%contribution des terminaisons 

image = imread(emplacement_image);
[H,L] = size(image);%dimensions de l'image : Hauteur et Largeur

%affichage de l'image initiale
figure(1)
colormap gray;
imagesc(image);

image_lissee = double(imgaussfilt(image,sigma));
figure(2)
colormap gray;
imagesc(image_lissee);
title("Image filtrée")

%Computing external forces

eline = image_lissee; %eline is simply the image intensities

[grady,gradx] = gradient(image_lissee);

eedge = -1 * sqrt ((gradx .* gradx + grady .* grady)); %eedge is measured by gradient in the image

%masks for taking various derivatives
m1 = [-1 1];
m2 = [-1;1];
m3 = [1 -2 1];
m4 = [1;-2;1];
m5 = [1 -1;-1 1];

cx = conv2(image_lissee,m1,'same');
cy = conv2(image_lissee,m2,'same');
cxx = conv2(image_lissee,m3,'same');
cyy = conv2(image_lissee,m4,'same');
cxy = conv2(image_lissee,m5,'same');

eterm = NaN(H,L);%preallocation
for i = 1:H
    for j= 1:L
        % eterm as deined in Kass et al Snakes paper
        eterm(i,j) = (cyy(i,j)*cx(i,j)*cx(i,j) -2 *cxy(i,j)*cx(i,j)*cy(i,j) + cxx(i,j)*cy(i,j)*cy(i,j))/((1+cx(i,j)*cx(i,j) + cy(i,j)*cy(i,j))^1.5);
    end
end

figure(3)
imagesc(eterm);
title("eterm");

figure(4)
imagesc(abs(eedge));
title("abs(eedge)");

eext = (weline*eline + weedge*eedge -weterm * eterm); %eext as a weighted sum of eline, eedge and eterm

figure(5)
colormap turbo;
imagesc(eext);
colorbar
title("Energie externe à minimiser")

figure(6)
surf(eext,'EdgeColor','interp','FaceColor','interp');
colormap turbo;