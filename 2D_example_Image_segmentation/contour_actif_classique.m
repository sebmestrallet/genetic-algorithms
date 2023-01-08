%basé sur activeContoursSnakesDemo, mais sans interface graphique

clear variables
close all

emplacement_image = "tile.tif";

%parametre du filtre
sigma = 2;
%parametres du snake
alpha = 0.08;
beta = 0.01;
gamma = 1;
kappa = 0.2;
%parametres GVF
weline = 0;%-0.2
weedge = 1.2;
weterm = 0;%0.2
iterations = 400;%nb iterations

%ouverture de l'image --------------------
image = imread(emplacement_image);
[row,col] = size(image);

figure
colormap gray;
imagesc(image);

%Filtrage ---------------------
%image_filtered = filter_function(image,sigma);
image_filtered = double(imgaussfilt(image,sigma));

%Snake initial
%ici le snake initial est prédéfini, mais on peut appeler la fonction
%getsnake.m de activeContoursSnakesDemo pour le choisir (graphiquement) à chaque execution

xs=[76.0817 86.9347 96.9514 106.1726 114.6387 122.3904 129.4684 135.9132 141.7655 147.0658 151.8548 156.1731...
	160.0613 163.5600 166.7098 169.5514 172.1253 174.4722 176.6327 178.6473 180.5567 182.3866 184.1027 185.6557...
    186.9965 188.0758 188.8445 189.2532 189.2529 188.7942 187.8279 186.3257 184.3426 181.9546 179.2374 176.2671...
    173.1195 169.8705 166.5960 163.3719 160.2740 157.3410 154.4614 151.4865 148.2675 144.6557 140.5025 135.6589...
    129.9763 123.3059 115.4990 106.4822 96.4834 85.8058 74.7526 63.6272 52.7328 42.3726 32.8499 24.4678...
    17.5298 12.2592 8.5601 6.2570 5.1742 5.1361 5.9670 7.4913 9.5334 11.9176 14.4683 17.0406...
    19.6127 22.1936 24.7921 27.4174 30.0782 32.7836 35.5426 38.3639 41.2567 44.2299 47.2924 50.4531...
    53.7211 57.1052 60.6145 64.2578 68.0442 71.9825 76.0817];

ys=[14.1769 23.8426 31.8539 38.3477 43.4611 47.3311 50.0950 51.8896 52.8521 53.1195 52.8288 52.1173...
	51.1218 49.9794 48.8273 47.8025 47.0420 46.6829 46.8623 47.7172 49.3846 51.9573 55.3499 59.4329...
    64.0767 69.1516 74.5280 80.0762 85.6667 91.1697 96.4558 101.4218 106.0712 110.4344 114.5412 118.4221...
	122.1070 125.6261 129.0096 132.2877 135.4904 138.6310 141.6551 144.4914 147.0685 149.3151 151.1599 152.5316...
	153.3587 153.5700 153.0942 151.8830 149.9801 147.4526 144.3672 140.7909 136.7906 132.4331 127.7853 122.9142...
	117.8865 112.7626 107.5755 102.3518 97.1181 91.9008 86.7264 81.6215 76.6125 71.7260 66.9885 62.4226...
	58.0361 53.8329 49.8168 45.9917 42.3615 38.9301 35.7013 32.6791 29.8673 27.2698 24.8905 22.7333...
	20.8020 19.1006 17.6329 16.4028 15.4142 14.6709 14.1769];

%Computing external forces

eline = image_filtered; %eline is simply the image intensities

[grady,gradx] = gradient(image_filtered);

eedge = -1 * sqrt ((gradx .* gradx + grady .* grady)); %eedge is measured by gradient in the image
%masks for taking various derivatives
m1 = [-1 1];
m2 = [-1;1];
m3 = [1 -2 1];
m4 = [1;-2;1];
m5 = [1 -1;-1 1];

cx = conv2(image_filtered,m1,'same');
cy = conv2(image_filtered,m2,'same');
cxx = conv2(image_filtered,m3,'same');
cyy = conv2(image_filtered,m4,'same');
cxy = conv2(image_filtered,m5,'same');

eterm = NaN(row,col);%preallocation
for i = 1:row
    for j= 1:col
        % eterm as deined in Kass et al Snakes paper
        eterm(i,j) = (cyy(i,j)*cx(i,j)*cx(i,j) -2 *cxy(i,j)*cx(i,j)*cy(i,j) + cxx(i,j)*cy(i,j)*cy(i,j))/((1+cx(i,j)*cx(i,j) + cy(i,j)*cy(i,j))^1.5);
    end
end

% imview(eterm);
% imview(abs(eedge));

eext = (weline*eline + weedge*eedge -weterm * eterm); %eext as a weighted sum of eline, eedge and eterm

[fx, fy] = gradient(eext); %computing the gradient

%initializing the snake
xs=xs';
ys=ys';
[m n] = size(xs);
[mm nn] = size(fx);

%populating the penta diagonal matrix
A = zeros(m,m);
b = [(2*alpha + 6 *beta) -(alpha + 4*beta) beta];
brow = zeros(1,m);
brow(1,1:3) = brow(1,1:3) + b;
brow(1,m-1:m) = brow(1,m-1:m) + [beta -(alpha + 4*beta)]; % populating a template row
for i=1:m
    A(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end

[L U] = lu(A + gamma .* eye(m,m));
Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization


for i=1:iterations%pour chaque iteration
    
    ssx = gamma*xs - kappa*interp2(fx,xs,ys);
    ssy = gamma*ys - kappa*interp2(fy,xs,ys);
    
    %calculating the new position of snake
    xs = Ainv * ssx;
    ys = Ainv * ssy;
    
    
    %Displaying the snake in its new position
    imagesc(image);
    hold on;
    plot([xs; xs(1)], [ys; ys(1)], 'r-','LineWidth',1);
    hold off;
    
    pause(0.001)
end