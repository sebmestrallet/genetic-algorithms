%fonctionnalites codees en dur
%- codage en coordonnees polaires des points du snake, a thetas implicites (seul
%  rho code)
%- selection : elistime (-> implementer roue biaisee et K-tournois ?)
%- croisement : par coupure ou croisement uniforme a proba variable
%               (-> implementer croisement à k site ?)

clear variables
close all

%PARAMETRES -------------------------------------------

%IMAGE A SEGMENTER
emplacement_image = "tile.tif";

%PARAMETRES DE LA FONCTION D'ENERGIE
%parametre du filtre gaussien
sigma = 10;
coarse_to_fine.sigmas = [5 3 1];
coarse_to_fine.generations = [20 40 60];
%parametres de l'energie externe
energie_externe.coeffs.intensite = 0;%contribution de l'intensite 
energie_externe.coeffs.bord = 1;%contribution des bords
energie_externe.coeffs.terminaison = 0;%contribution des terminaisons
%penalites sur rho pour simuler une energie interne de regularite
penalites_rho.mean = 0;
penalites_rho.max = 0;
penalites_rho.std = 0;
penalites_rho.somme_diff = -0.03;

%PARAMETRE DE L'AG
taille_population = 100;%nb de snake par generation
nb_points_par_snake = 40;
calcul_nb_bits_auto = true;%calcul de nb_bits d'après H et L
nb_bits = 6;%valeur si calcul_nb_bits_auto==false
croisement.proba = 0.8;%appliquee a deux parents
croisement.methode = "uniforme";% "1point" ou "uniforme"
croisement.uniforme.proba = 0.5;%si croisement uniforme
mutation_proba = 0.01;%appliquee a chaque gene de chaque enfant
nb_individus_selectionnes = floor(taille_population/2);
nb_individus_remplaces = floor(taille_population/2);%nb enfants qui remplacent des parents
population_initiale.rho_moy = 3/4;%entre 0 et 1, relatif au rayon minimal
population_initiale.dispertion = 5;
max_generations = 200;
afficher_chaque_generation = true;%si false, n'affichera que la derniere generation 
%hybridation avec voisinage
voisinage.theta.utiliser = false;
voisinage.theta.periode = 10;%en generations
voisinage.rho.utiliser = true;
voisinage.rho.periode = 10;%en generations
voisinage.rho.ecart_au_voisin = 5;%pixels
afficher_moyenne_genes = false;%affichage a la fin de la moyenne des genes
%critere de convergence
convergence_atteinte = @(scores,num_generation) (std(scores) <= 1) || (num_generation > max_generations);

% INITIALISATION : OUVERTURE DE L'IMAGE ET FILTRAGE ------------------
%ouverture de l'image
image = imread(emplacement_image);
[H,L] = size(image);%dimensions de l'image : Hauteur et Largeur

%affichage de l'image initiale
figure(1)
colormap gray;
imagesc(image);
title("Image à segmenter")

%Filtrage ---------------------
image_floutee = double(imgaussfilt(image,sigma));
figure(2)
colormap gray;
imagesc(image_floutee);
title(sprintf("Image après le filtre gaussien, sigma=%i",sigma));

% INITIALISATION : CALCUL DES FORCES EXTERNES ------------------

energie_externe.matrice = calcul_energie_externe(image_floutee,energie_externe.coeffs);

figure(3)
colormap turbo;
imagesc(energie_externe.matrice);
colorbar
title(sprintf("Energie externe à minimiser, sigma=%i",sigma));

% INITIALISATION : POPULATION INITIALE ------------------

if calcul_nb_bits_auto
    %calcul du nombre minimal de bits pour coder les pixels sur une
    %demi-diagonale
    pixels_sur_demi_diagonale = sqrt( (H/2)^2 + (L/2)^2 );
    nb_bits = ceil(log2(pixels_sur_demi_diagonale));
end

%repartition des thetas sur [0:2pi]
delta_theta = (2*pi)/nb_points_par_snake;%pas en radians
thetas = 0:delta_theta:2*pi-delta_theta;

num_generation = 0;%compteur de generation

%generation d'une population de snakes aléatoires avec rhos equiprobables
%population = randi(2,taille_population,nb_bits*nb_points_par_snake)-1;
%randi(2)   -> 1 ou 2
%randi(2)-1 -> 0 ou 1

%generation selon une loi normale sur le rayon
population = NaN(taille_population,nb_bits*nb_points_par_snake);%preallocation
for num_individu = 1:taille_population %pour chaque individu
    for num_point = 1:nb_points_par_snake %pour chaque point
        theta = thetas(num_point);
        val_max_rho = rho_max(L,H,theta);
        
        %forme rectangulaire -> relatif à val_max_rho (depend de theta)
        %rho = population_initiale.rho_moy*val_max_rho + population_initiale.dispertion*randn();
        
        %forme circulaire -> relatif au rayon min
        rho = population_initiale.rho_moy*min(L/2,H/2) + population_initiale.dispertion*randn();
        
        if rho > val_max_rho
            rho = val_max_rho;%rho ne doit pas depasser val_max_rho pour le codage
        elseif rho < 0
            rho = 0;%rho ne doit pas etre negatif pour le codage
        end
        code_gray = point_polaire_vers_binaire(rho,val_max_rho,nb_bits);
        population(num_individu,index_point(num_point,nb_bits)) = code_gray;%enregistrement dans la matrice population
    end
end

%calcul des scores de la population initiale
scores = calcul_scores_population(population,energie_externe.matrice,thetas,nb_bits,L,H,penalites_rho);

affichage_stats_scores(scores,num_generation);

while ~convergence_atteinte(scores,num_generation)
    
    num_generation = num_generation+1;
    
    %coarse to fine de l'energie externe
    if length(coarse_to_fine.generations)>=1 && num_generation==coarse_to_fine.generations(1)
        %s'il reste des valeurs dans coarse_to_fine.generations et qu'on se
        %trouve a la prochaine generation ou il faut adapter le sigma
        sigma = coarse_to_fine.sigmas(1);%nouveau sigma
        fprintf("\nNouveau sigma : %i, recalcul de l'énergie externe\n",sigma);
        image_floutee = double(imgaussfilt(image,sigma));
        energie_externe.matrice = calcul_energie_externe(image_floutee,energie_externe.coeffs);
        
        %mise a jour de la fenetre de l'image floutee
        figure(2)
        colormap gray;
        imagesc(image_floutee);
        title(sprintf("Image après le filtre gaussien, sigma=%i",sigma));
        
        %mise a jour de le fenetre avec l'energie externe
        figure(3)
        colormap turbo;
        imagesc(energie_externe.matrice);
        colorbar
        title(sprintf("Energie externe à minimiser, sigma=%i",sigma));
        
        %on supprime la colonne 1
        if length(coarse_to_fine.generations)==1
            coarse_to_fine.generations = [];
            coarse_to_fine.sigmas = [];
        else
            coarse_to_fine.generations = coarse_to_fine.generations(2:end);
            coarse_to_fine.sigmas = coarse_to_fine.sigmas(2:end);
        end
    end
    
    %---------- SELECTION --------
    %sélection de nb_individus_selectionnes individus selon leur score
    parents = selection(population,scores,nb_individus_selectionnes);
    
    %---------- CROISEMENTS --------
    %pour chaque paire de parents
    
    parents_melanges = melange_parents(parents);
    enfants = zeros(size(parents_melanges));
    for num_individu = 1:2:size(parents_melanges,1)-mod(size(parents_melanges,1),2)%de 2 en 2
        
        parent1 = parents_melanges(num_individu,:);
        parent2 = parents_melanges(num_individu+1,:);
        
        if rand()<croisement.proba
            %alors croisement
            
            if croisement.methode == "1point"
                [enfant1,enfant2] = croisement_1point(parent1,parent2);
            elseif croisement.methode == "uniforme"
                [enfant1,enfant2] = croisement_uniforme(parent1,parent2,croisement.uniforme.proba);
            else
                error("Methode de croisement inconnue (variable croisement.methode)");
            end
            
            %ecriture dans le tableau enfants
            enfants(num_individu,:) = enfant1;
            enfants(num_individu+1,:) = enfant2;
        else
            %pas de croisement
            %ecriture dans le tableau enfants
            enfants(num_individu,:) = parent1;
            enfants(num_individu+1,:) = parent2;
        end
    end
    if mod(size(parents_melanges,1),2) ~= 0
       %le nombre de parents n'est pas pair, le parent orphelin est clone
       enfants(end,:) = parents_melanges(end,:);
    end
    
    %---------- MUTATIONS --------
    for num_individu = 1:size(enfants,1)%pour chaque enfant
        for num_gene = 1:nb_bits%pour chaque gene
            if rand()<mutation_proba
                enfants(num_individu,num_gene) = double(~enfants(num_individu,num_gene));%complement
            end
        end
    end
    
    %---------- REMPLACEMENT --------
    %nouvelle population
    %elitisme : les meilleurs enfants sont gardes
    %et remplace les moins bons individus de la génération d'avant
    if nb_individus_selectionnes < nb_individus_remplaces
       error("Erreur : il y a plus de parents à remplacer que d'enfant\n");
    end
    index_des_meilleurs_enfants = index_des_meilleurs_individus(calcul_scores_population(enfants,energie_externe.matrice,thetas,nb_bits,L,H,penalites_rho),nb_individus_remplaces);
    index_des_pires_parents = index_des_pires_individus(scores,nb_individus_remplaces);
    for i = 1:nb_individus_remplaces
        %remplacement dans la population
        %d'un des moins bon par un des meilleurs enfants
        population(index_des_pires_parents(i),:) = enfants(index_des_meilleurs_enfants(i),:);
    end
    
    %---------- VOISINAGE -----------------
    %/!\ l'ordre dans lequel sont appliques les voisinages a une influence
    if voisinage.theta.utiliser && mod(num_generation,voisinage.theta.periode)==0
        fprintf("\nApplication d'une étape de voisinage selon theta\n");
        population = voisinage_selon_thetas(population,energie_externe.matrice,thetas,nb_bits,nb_points_par_snake,L,H);
    end
    if voisinage.rho.utiliser && mod(num_generation,voisinage.rho.periode)==0
        fprintf("\nApplication d'une étape de voisinage selon rho\n");
        population = voisinage_selon_rhos(population,voisinage.rho.ecart_au_voisin,energie_externe.matrice,thetas,nb_bits,nb_points_par_snake,L,H);
    end
    
    %calcul des nouveaux scores
    scores = calcul_scores_population(population,energie_externe.matrice,thetas,nb_bits,L,H,penalites_rho);
    
    affichage_stats_scores(scores,num_generation);
    
    if afficher_chaque_generation
        [~,index_du_max] = max(scores);
        meilleur_individu = population(index_du_max,:);
        [snake_x,snake_y,~] = snake_binaire_vers_xy(meilleur_individu,thetas,nb_bits,L,H);
        figure(4)
        affichage_image_et_snake(image,snake_x,snake_y,sprintf("Meilleur snake de la génération %i",num_generation));
    end
    
    %waitforbuttonpress;
    %pause(0.1);
end

fprintf("\nCritère de convergence atteint\n");

if ~afficher_chaque_generation
	[~,index_du_max] = max(scores);
	meilleur_individu = population(index_du_max,:);
	[snake_x,snake_y] = snake_binaire_vers_xy(meilleur_individu,thetas,nb_bits,L,H);
    figure(4)
	affichage_image_et_snake(image,snake_x,snake_y,"Meilleur snake de la dernière génération");
end

if afficher_moyenne_genes
    figure(5)
    bar(mean(population,1),"k")%,"LineWidth",1
    xlim([0 size(population,2)])
    ylim([0 1])
    title("Moyenne des gènes")
    xlabel("Numéro de gène")
    ylabel("Moyenne")
end

function valeur = rho_max(L,H,theta)
    a = L/2;
    b = H/2;
    if abs(tan(theta)) <= b/a
        valeur = a/abs(cos(theta));
    else
        valeur = b/abs(sin(theta));
    end
end

function affichage_image_et_snake(I,x,y,titre)
    colormap gray;
    imagesc(I);
    hold on
    %segments entre points
    for i = 1:length(x)-1
        line([x(i) x(i+1)],[y(i) y(i+1)],"Color","red","LineWidth",1);
    end
    %segment entre le premier et le dernier
    line([x(1) x(end)],[y(1) y(end)],"Color","red","LineWidth",1);
    title(titre);
    hold off
    drawnow% limitrate
end

function [x,y,rho] = point_binaire_vers_xy(code_gray,theta,L,H)
    nb_bits = length(code_gray);

	valeur_decimale = gc2dec(code_gray);
	%(rappel : la val max en code gray 7 bits est [1 0 0 0 0 0 0])

	%normalisation selon rho max
	rho = ( valeur_decimale/(2^nb_bits-1) );%entre 0 et 1
    rho = rho * rho_max(L,H,theta);%entre 0 et rho_max

	%conversion en coordonnees cartesiennes
	[x,y] = pol2cart(theta,rho);
    %le point (0,0) est au milieu en coordonnees polaires, mais en haut a
    %gauche en coordonnees cartesiennes
    x = x+L/2;%ajout d'une demi largeur
	y = y+H/2;%ajout d'une demi hauteur
	%maintenant (0,0) est bien en haut a gauche de l'image

end

function gray_code = point_polaire_vers_binaire(rho,rho_max,nb_bits)
    rho_normalise = rho/rho_max;%entre 0 et 1;
    valeur_decimale = round(rho_normalise*(2^nb_bits-1));%entre 0 et 2^nb_bits
    gray_code = dec2gc(valeur_decimale,nb_bits);
end

function [x_vecteur,y_vecteur,rho_vecteur] = snake_binaire_vers_xy(snake,thetas,nb_bits,L,H)

    %deduction du nb de points depuis la longueur de theta
    nb_points_par_snake = length(thetas);

    %preallocation
    x_vecteur = NaN(1,nb_points_par_snake);%coordonnees x des points du snake
    y_vecteur = NaN(1,nb_points_par_snake);%coordonnees y des points du snake
    rho_vecteur = NaN(1,nb_points_par_snake);%magnitudes des points du snake

    %pour chaque point du snake
    for num_point = 1:nb_points_par_snake
        [x_vecteur(num_point),y_vecteur(num_point),rho_vecteur(num_point)] = point_binaire_vers_xy(snake(index_point(num_point,nb_bits)),...
                                             thetas(num_point),...
                                             L,...
                                             H);
    end
end

function score = calcul_score_point(x,y,matrice_energie_ext)
    try
        x = round(x);
        %limitation entre 1 et le nombre de colonnes dans la matrice
        if x < 1
            x = 1;
        elseif x > size(matrice_energie_ext,2)
            x = size(matrice_energie_ext,2);
        end
        y = round(y);
        %limitation entre 1 et le nombre de lignes dans la matrice
        if y < 1
            y = 1;
        elseif y > size(matrice_energie_ext,1)
            y = size(matrice_energie_ext,1);
        end
        score = -matrice_energie_ext(y,x);%un bon score est une énergie faible
    catch
        fprintf("x = %d, y = %d\n",x,y);
        error("problème d'index, fonction compute_score_point\n");
    end
end

function score = calcul_score_snake(x_vecteur,y_vecteur,rho_vecteur,penalites_rho,matrice_energie_ext)
    score = 0;
    for i = 1:length(x_vecteur)%pour chaque point du snake
        score = score + calcul_score_point(x_vecteur(i),y_vecteur(i),matrice_energie_ext);
    end
    
    %ajout d'une composante sur la std des magnitudes
    score = score + penalites_rho.std*std(rho_vecteur);
    
    %ajout d'une composante sur la somme des differences successives de
    %magnitudes
    rho_vector_extended = [rho_vecteur rho_vecteur(1)];
    differences = diff(rho_vector_extended);
    score = score + penalites_rho.somme_diff*sum(differences.^2);
    %-> privilegie les snakes reguliers, sans zigzag
    
    %ajout d'une composante sur la moyenne des magnitudes
    score = score + penalites_rho.mean*mean(rho_vecteur);
    
    %ajout d'une composante sur la moyenne des magnitudes
    score = score + penalites_rho.max*max(rho_vecteur);
    %-> privilegie les snakes proche du centre de l'image
end

function scores = calcul_scores_population(population,matrice_energie_ext,thetas,nb_bits,L,H,penalites_rho)

    scores = NaN(size(population,1),1);%vecteur colonne avec taille_population éléments
    
    for num_individu = 1:size(population,1)%pour chaque individu
        
        %calcul des coordonnes cartesiennes
        [x_vector,y_vector,rho_vector] = snake_binaire_vers_xy(population(num_individu,:),thetas,nb_bits,L,H);
        
        %calcul du score
        scores(num_individu) = calcul_score_snake(x_vector,y_vector,rho_vector,penalites_rho,matrice_energie_ext);
    end
end

function affichage_stats_scores(scores,num_generation)
    fprintf("\nGénération %i ----------------\n",num_generation);
    fprintf("Score moyen :    %f\n",mean(scores));
    fprintf("Meilleur score : %f\n",max(scores));
    fprintf("Ecart-type :     %f\n",std(scores));
end

function individus_selectionnes = selection(population,scores,nb_individus_a_choisir)    
    scores_et_population = [scores population];
    scores_et_population_trie = sortrows(scores_et_population,1,'descend');%tri des individus par score
    individus_selectionnes = scores_et_population_trie(1:nb_individus_a_choisir,:);%on ne garde que les nb_child premiers
    individus_selectionnes = individus_selectionnes(:,2:end);%on supprime la colonne des scores
end

function parents_melanges = melange_parents(parents)
    new_index = randperm(size(parents,1)); % permutation des index
    parents_melanges = parents(new_index,:);
end

function [enfant1,enfant2] = croisement_1point(parent1,parent2)
    enfant1 = parent1;
    enfant2 = parent2;
    
    %tirage au sort d'un site de croisement
    point_croisement = randi(length(enfant1)-1);%entre 1 compris et length(child1)-1 compris
    
    temp = enfant1(1:point_croisement);
    enfant1(1:point_croisement) = enfant2(1:point_croisement);
    enfant2(1:point_croisement) = temp;
end

%function croisement_multipoints
%random_index = randperm(1:length(parent1));
%crossing_points = random_index(1:nb_points_coupure);

function [enfant1,enfant2] = croisement_uniforme(parent1,parent2,proba_point_croisement)

    %preallocation
    enfant1 = NaN(size(parent1));
    enfant2 = NaN(size(parent2));
    
    masque = double(rand(1,length(parent1))>proba_point_croisement);
    %ex avec proba_point_croisement=0.1
    % si rand > proba_point_croisement -> masque à 1 = pas permutation
    % si rand <= proba_point_croisement -> masque à 0 = permutation
    %on aura bien 10% de permutation en moyenne
    
    for i = 1:length(parent1)%pour chaque num de gene
        if masque(i) == 1
            enfant1(i) = parent1(i);
            enfant2(i) = parent2(i);
        else
            enfant1(i) = parent2(i);
            enfant2(i) = parent1(i);
        end
    end
    
%     %affichage detaille
%     fprintf("%i",mask);
%     fprintf("\n")
%     fprintf("%i",parent1);
%     fprintf("\n")
%     fprintf("%i",parent2);
%     fprintf("\n")
%     fprintf("%i",child1);
%     fprintf("\n")
%     fprintf("%i",child2);
%     fprintf("\n")
%     error("fin")
end

function index = index_des_meilleurs_individus(scores,nb_individus_a_choisir)
    index_originaux = (1:size(scores,1))';
    index_et_scores = [index_originaux scores];
    sorted = sortrows(index_et_scores,2,'descend');%tri des individus par score
    individus_selectionnes = sorted(1:nb_individus_a_choisir,:);%on ne garde que les number premiers
    index = individus_selectionnes(:,1);%1ere colonne (index originaux)
end

function index = index_des_pires_individus(scores,nb_individus_a_choisir)
    index_originaux = (1:size(scores,1))';
    index_et_scores = [index_originaux scores];
    sorted = sortrows(index_et_scores,2,'ascend');%tri des individus par score
    individus_selectionnes = sorted(1:nb_individus_a_choisir,:);%on ne garde que les number premiers
    index = individus_selectionnes(:,1);%1ere colonne (index originaux)
end

function index = index_point(num_point,nb_bits)
	index = (num_point-1)*nb_bits+1:num_point*nb_bits;
end

function nouvelle_popupation = voisinage_selon_thetas(population,matrice_energie_ext,thetas,nb_bits,nb_points_par_snake,L,H)
    nouvelle_popupation = NaN(size(population));
    for num_individu = 1:size(population,1) %pour chaque individu
        snake = population(num_individu,:);
        [x_vector,y_vector,~] = snake_binaire_vers_xy(snake,thetas,nb_bits,L,H);

        %on commence par calculer et stocker le score de chaque point (pas
        %d'energie interne)
        scores_vecteur = NaN(size(x_vector));
        for index = 1:nb_points_par_snake %pour chaque point
            %on calcule son score
            scores_vecteur(index) = calcul_score_point(x_vector(index),y_vector(index),matrice_energie_ext);
        end

        for index = 1:nb_points_par_snake%pour chaque point
            %point precedent
            if index==1
                index_du_precedent = nb_points_par_snake;
            else
                index_du_precedent = index-1;
            end
            %point suivant
            if index==nb_points_par_snake
                index_du_suivant = 1;
            else
                index_du_suivant = index+1;
            end

            %comparaison
            if scores_vecteur(index_du_precedent) > scores_vecteur(index) && scores_vecteur(index_du_precedent) > scores_vecteur(index_du_suivant)
                %le point precedent est le meilleur des trois
                code_gray_du_precedent = population(num_individu,index_point(index_du_precedent,nb_bits));
                [~,~,rho] = point_binaire_vers_xy(code_gray_du_precedent,thetas(index_du_precedent),L,H);
                rho = min(rho,rho_max(L,H,thetas(index)));
                nouvelle_popupation(num_individu,index_point(index,nb_bits)) = point_polaire_vers_binaire(rho,rho_max(L,H,thetas(index)),nb_bits);
            elseif scores_vecteur(index_du_suivant) > scores_vecteur(index) && scores_vecteur(index_du_suivant) > scores_vecteur(index_du_precedent)
                %le point suivant est le meilleur des trois
                code_gray_du_suivant = population(num_individu,index_point(index_du_precedent,nb_bits));
                [~,~,rho] = point_binaire_vers_xy(code_gray_du_suivant,thetas(index_du_suivant),L,H);
                rho = min(rho,rho_max(L,H,thetas(index)));
                nouvelle_popupation(num_individu,index_point(index,nb_bits)) = point_polaire_vers_binaire(rho,rho_max(L,H,thetas(index)),nb_bits);
            else
                %le point actuel est le meilleur des trois
                nouvelle_popupation(num_individu,index_point(index,nb_bits)) = population(num_individu,index_point(index,nb_bits));
            end
        end
    end
end

function nouvelle_popupation = voisinage_selon_rhos(population,ecart_pixels,matrice_energie_ext,thetas,nb_bits,nb_points_par_snake,L,H)
    nouvelle_popupation = NaN(size(population));
    for num_individu = 1:size(population,1) %pour chaque individu
        snake = population(num_individu,:);
        [x_vector,y_vector,~] = snake_binaire_vers_xy(snake,thetas,nb_bits,L,H);

        %on commence par calculer et stocker le score de chaque point (pas
        %d'energie interne)
        scores_vecteur = NaN(size(x_vector));
        for index = 1:nb_points_par_snake %pour chaque point
            %on calcule son score
            scores_vecteur(index) = calcul_score_point(x_vector(index),y_vector(index),matrice_energie_ext);
        end

        for index = 1:nb_points_par_snake%pour chaque point
            valeur_decimale = gc2dec(snake(index_point(index,nb_bits)));

            %point interieur
            valeur_decimale_inferieure = max(valeur_decimale-ecart_pixels,1);
            %point exterieur
            valeur_decimale_superieure = min(valeur_decimale+ecart_pixels,2^nb_bits-1);

            valeur_binaire_inferieure = dec2gc(valeur_decimale_inferieure,nb_bits);
            valeur_binaire_superieure = dec2gc(valeur_decimale_superieure,nb_bits);
            [x_inf,y_inf,~] = point_binaire_vers_xy(valeur_binaire_inferieure,thetas(index),L,H);
            [x_sup,y_sup,~] = point_binaire_vers_xy(valeur_binaire_superieure,thetas(index),L,H);
            score_point_interieur = calcul_score_point(x_inf,y_inf,matrice_energie_ext);
            score_point_exterieur = calcul_score_point(x_sup,y_sup,matrice_energie_ext);

            %comparaison
            if score_point_interieur > scores_vecteur(index) && score_point_interieur > score_point_exterieur
                %le point interieur est le meilleur des trois
                nouvelle_popupation(num_individu,index_point(index,nb_bits)) = valeur_binaire_inferieure;
            elseif score_point_exterieur > scores_vecteur(index) && score_point_exterieur > score_point_interieur
                %le point exterieur est le meilleur des trois
                nouvelle_popupation(num_individu,index_point(index,nb_bits)) = valeur_binaire_superieure;
            else
                %le point actuel est le meilleur des trois
                nouvelle_popupation(num_individu,index_point(index,nb_bits)) = population(num_individu,index_point(index,nb_bits));
            end
        end
    end
end

function energie_externe = calcul_energie_externe(image_floutee,coeffs)
    [H,L] = size(image_floutee);

    %Computing external forces

    eline = image_floutee; %eline is simply the image intensities

    [grady,gradx] = gradient(image_floutee);

    eedge = -1 * sqrt ((gradx .* gradx + grady .* grady)); %eedge is measured by gradient in the image

    %masks for taking various derivatives
    m1 = [-1 1];
    m2 = [-1;1];
    m3 = [1 -2 1];
    m4 = [1;-2;1];
    m5 = [1 -1;-1 1];

    cx = conv2(image_floutee,m1,'same');
    cy = conv2(image_floutee,m2,'same');
    cxx = conv2(image_floutee,m3,'same');
    cyy = conv2(image_floutee,m4,'same');
    cxy = conv2(image_floutee,m5,'same');

    %preallocation
    eterm = NaN(H,L);
    for i = 1:H
        for j= 1:L
            % eterm as deined in Kass et al Snakes paper
            eterm(i,j) = (cyy(i,j)*cx(i,j)*cx(i,j) -2 *cxy(i,j)*cx(i,j)*cy(i,j) + cxx(i,j)*cy(i,j)*cy(i,j))/((1+cx(i,j)*cx(i,j) + cy(i,j)*cy(i,j))^1.5);
        end
    end

    % imview(eterm);
    % imview(abs(eedge));

    energie_externe = (coeffs.intensite*eline + coeffs.bord*eedge -coeffs.terminaison * eterm); %eext as a weighted sum of eline, eedge and eterm
end