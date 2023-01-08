clear all
close all

taille_population = 20;%doit être multiple de 4 (/2 -> enfants, /2 -> paires d'enfants)
proba_croisement = 0.8;
proba_mutation = 0.08;%initialement 0.01
nb_individus_remplaces = 5;
afficher_etapes = false;
afficher_evolution_du_score = false;

%compteur de générations
num_generation = 0;
max_generations = 100;
%preallocation
meilleurs_scores = NaN(1,max_generations+1);%+1 case pour la génération 0
scores_moyens = NaN(1,max_generations+1);

%création de la population
%une ligne pour chaque individu, une colonne pour chaque gène (8 bits ici)
%population = randi(2,taille_population,8)-1;
population = [1 1 0 1 0 0 0 1;
              1 0 0 0 1 1 1 1;
              0 1 1 1 1 0 0 0;
              1 1 1 1 1 1 1 1;
              1 1 0 1 1 0 1 0;
              0 1 0 1 0 1 1 0;
              0 1 0 1 1 0 0 0;
              1 0 1 0 1 1 0 0;
              1 1 1 0 0 1 0 0;
              1 0 1 0 1 1 1 0;
              0 1 0 1 0 0 0 0;
              1 0 1 0 0 0 1 0;
              1 0 1 1 1 0 0 1;
              0 0 0 0 1 1 0 1;
              1 0 0 1 1 0 0 0;
              0 1 0 0 0 1 0 0;
              0 1 1 0 1 1 1 0;
              1 0 0 0 0 1 1 1;
              1 1 1 1 0 0 1 0;
              1 0 0 0 0 0 0 0];

%calcul des scores
scores = calcul_scores_population(population);

figure
if afficher_etapes
    fprintf("Population initiale\n");
    afficher_population_et_scores(population,scores);
    waitforbuttonpress;
end

afficher_stats_scores(scores,num_generation);
scores_moyens(1,num_generation+1) = mean(scores);
meilleurs_scores(1,num_generation+1) = max(scores);

%affichage des individus
tracer_population(population,sprintf("Population à la génération %i",num_generation));

%critère de convergence : population homogène
%écart-type faible
while (std(scores) > 0.01) && (num_generation <= max_generations)

    num_generation = num_generation+1;
    
    %---------- SELECTION --------
    %sélection de taille_population/2 individus selon leur score
    parents = selection(population,scores);
    if afficher_etapes
        fprintf("\nParents sélectionnés (%i meilleurs)\n",taille_population/2);
        afficher_population_et_scores(parents,calcul_scores_population(parents));
        waitforbuttonpress;
    end
    
    %---------- CROISEMENTS --------
    %pour chaque paire de parents
    if mod(size(parents,1),2) ~= 0
       error("Le nombre de parents n'est pas pair");
       %TODO : à la place, juste cloner l'individu orphelin
    end
    parents_melanges = melange_parents(parents);
    if afficher_etapes
        fprintf("\nParents mélangés\n");
        afficher_population_et_scores(parents_melanges,calcul_scores_population(parents_melanges));
        waitforbuttonpress;
    end
    enfants = zeros(size(parents_melanges));
    for i = 1:2:size(parents_melanges,1)%de 2 en 2
        parent1 = parents_melanges(i,:);
        parent2 = parents_melanges(i+1,:);
        if afficher_etapes
            fprintf("\nPaire des parents %i&%i\n",i,i+1);
            fprintf("Parent 1 : ");
            afficher_genes_individu(parent1);
            fprintf("\n");
            fprintf("Parent 2 : ");
            afficher_genes_individu(parent2);
            fprintf("\n");
        end
        if rand()<proba_croisement
            %alors croisement.
            [enfant1,enfant2,point_croisement] = croisement(parent1,parent2);
            %ecriture dans le tableau enfants
            enfants(i,:) = enfant1;
            enfants(i+1,:) = enfant2;
            if afficher_etapes
                fprintf("Croisement au point %i\n",point_croisement);
                fprintf("Enfant 1 : ");
                afficher_genes_individu(enfant1);
                fprintf("\n");
                fprintf("Enfant 2 : ");
                afficher_genes_individu(enfant2);
                fprintf("\n");
                waitforbuttonpress;
            end
        else
            %ecriture dans le tableau enfants
            enfants(i,:) = parent1;
            enfants(i+1,:) = parent2;
            if afficher_etapes
                fprintf("Pas de croisement\n");
                waitforbuttonpress;
            end
        end
    end
    
    %---------- MUTATIONS --------
    %pour chaque gène
    for i = 1:size(enfants,1)%pour chaque enfant
        if afficher_etapes
            fprintf("\nEnfant %i\n",i);
            afficher_genes_individu(individu_depuis_index(enfants,i));
            fprintf("\n");
        end
        emplacement_mutations = zeros(1,8);
        for j = 1:8%pour chaque gène
            if rand()<proba_mutation
                emplacement_mutations(j) = 1;
                enfants(i,j) = double(~enfants(i,j));%complément
            end
        end
        if afficher_etapes
            afficher_genes_individu(individu_depuis_index(enfants,i));
            fprintf("\n");
            if sum(emplacement_mutations) ~= 0 %s'il y a eu une mutation
                for j = 1:8%pour chaque gène
                    if emplacement_mutations(j) == 1
                        fprintf("^");%marqueur de l'emplacement de la mutation
                    else
                        fprintf(" ");
                    end
                end
                fprintf("\n");
            else
                fprintf("Pas de mutation\n");
            end
            waitforbuttonpress;
        end
    end
    
    %---------- REMPLACEMENT --------
    %nouvelle population
    %elitisme : les meilleurs enfants sont gardés
    %et remplacent les moins bons individus de la génération d'avant
    index_des_meilleurs_enfants = index_des_meilleurs(calcul_scores_population(enfants),nb_individus_remplaces);
    index_des_pires_parents = index_des_pires(scores,nb_individus_remplaces);
    for i = 1:nb_individus_remplaces
        %remplacement dans la population
        %d'un des moins bons parents par un des meilleurs enfants
        population(index_des_pires_parents(i),:) = enfants(index_des_meilleurs_enfants(i),:);
    end
    if afficher_etapes
        fprintf("\nPopulation avant sélection, croisements et mutations\n");
        afficher_population_et_scores(population,scores)
        fprintf("index des pires :");
        fprintf("% d",index_des_pires_parents);%espacements entre les index
        fprintf("\nEnfants\n");
        afficher_population_et_scores(enfants,calcul_scores_population(enfants))
        fprintf("index des meilleurs :");
        fprintf("% d",index_des_meilleurs_enfants);%espacements entre les index
        fprintf("\n");
        waitforbuttonpress;
    end
    
    %calcul des nouveaux scores
    scores = calcul_scores_population(population);
    
    afficher_stats_scores(scores,num_generation);
    scores_moyens(1,num_generation+1) = mean(scores);
    meilleurs_scores(1,num_generation+1) = max(scores);
    
    if afficher_etapes
        fprintf("\nNouvelle population\n");
        afficher_population_et_scores(population,scores)
    end
    
    tracer_population(population,sprintf("Population à la génération %i",num_generation));
    %waitforbuttonpress;
    pause(0.5);
end

fprintf("\nCritère de convergence atteint\n");
[~,index_du_max] = max(scores);
fprintf("Meilleur individu : ");
meilleur_individu = population(index_du_max,:);
afficher_genes_individu(meilleur_individu);
fprintf("\n= %i\n",valeur_decimale(meilleur_individu));

if afficher_evolution_du_score
    %évolution des scores max et moyen au cours des générations
    figure
    scores_moyens = scores_moyens(1,1:num_generation+1);
    meilleurs_scores = meilleurs_scores(1,1:num_generation+1);
    plot(0:num_generation,scores_moyens);
    hold on
    plot(0:num_generation,meilleurs_scores);
    hold off
    xlabel("Générations");
    ylabel("Score");
    title("Evolution des scores moyen et max");
    legend("Score moyen","Score max","Location","best");
end

function y = fonctionnelle_a_minimiser(x)
    y = sin(2*pi*0.01*x).*(-x)*0.02 - 4;
end

function gray_code = code_binaire(decimal)
    gray_code = dec2gc(decimal,8);%8 bits
end

function decimal = valeur_decimale(binaire)
    decimal = gc2dec(binaire);
end

function score = calcul_score(x)
    score = -fonctionnelle_a_minimiser(x);
    %fonctionnelle minimale = score maximal
end

function individual = individu_depuis_index(population,index)
    individual = population(index,:);
end

function scores = calcul_scores_population(population)
    scores = zeros(size(population,1),1);%vecteur colonne avec taille_population éléments
    for individual_index = 1:size(population,1)%pour chaque individu
        scores(individual_index) = calcul_score(...
                                    valeur_decimale(...
                                     individu_depuis_index(population,individual_index)));
    end
end

function tracer_population(population,titre)
    x = 0:1:2^8;
    y = fonctionnelle_a_minimiser(x);
    plot(x,y,"LineWidth",1);
    hold on
    for index_individu = 1:size(population,1)%pour chaque individu
        val_decimal_individu = valeur_decimale(individu_depuis_index(population,index_individu));
        plot(val_decimal_individu,...
             fonctionnelle_a_minimiser(val_decimal_individu),...
             '.k',...
             "MarkerSize",20);
    end
    title(titre);
    hold off
    xlim([0 2^8]);
    xticks([0 64 128 192 256]);
end

function afficher_stats_scores(scores,num_generation)
    fprintf("\nGénération %i ----------------\n",num_generation);
    fprintf("Score moyen :    %f\n",mean(scores));
    fprintf("Meilleur score : %f\n",max(scores));
    fprintf("Ecart-type :     %f\n",std(scores));
end

function individus_selectionnes = selection(population,scores)
    nb_individus = size(population,1);
    if mod(nb_individus,2) ~= 0
       error("Le nombre d'individus n'est pas pair");
    end
    nb_enfants = nb_individus/2;
    
    scores_et_population = [scores population];
    scores_et_population_trie = sortrows(scores_et_population,1,'descend');%tri des individus par score
    individus_selectionnes = scores_et_population_trie(1:nb_enfants,:);%on ne garde que les nb_child premiers
    individus_selectionnes = individus_selectionnes(:,2:end);%on supprime la colonne des scores
end

function parents_melanges = melange_parents(parents)
    new_index = randperm(size(parents,1)); % permutation des index
    parents_melanges = parents(new_index,:);
end

function [enfant1,enfant2,point_croisement] = croisement(parent1,parent2)
    enfant1 = parent1;
    enfant2 = parent2;
    
    %tirage au sort d'un site de croisement
    point_croisement = randi(7);%entre 1 et 7 compris
    
    temp = enfant1(1:point_croisement);
    enfant1(1:point_croisement) = enfant2(1:point_croisement);
    enfant2(1:point_croisement) = temp;
end

function index = index_des_meilleurs(scores,quantite)
    index_originaux = (1:size(scores,1))';
    index_et_scores = [index_originaux scores];
    index_et_scores_trie = sortrows(index_et_scores,2,'descend');%tri des individus par score
    individus_selectionnes = index_et_scores_trie(1:quantite,:);%on ne garde que les quantite premiers
    index = individus_selectionnes(:,1);%1ere colonne (index originaux)
end

function index = index_des_pires(scores,quantite)
    index_originaux = (1:size(scores,1))';
    index_et_scores = [index_originaux scores];
    index_et_scores_trie = sortrows(index_et_scores,2,'ascend');%tri des individus par score
    individus_selectionnes = index_et_scores_trie(1:quantite,:);%on ne garde que les quantite premiers
    index = individus_selectionnes(:,1);%1ere colonne (index originaux)
end

function afficher_genes_individu(individual)
    fprintf("%i",individual);
end

function afficher_population_et_scores(population,scores)
    for index_individu = 1:size(population,1)%pour chaque individu
        fprintf("%2i | ",index_individu);
        afficher_genes_individu(individu_depuis_index(population,index_individu));
        fprintf(" | %f\n",scores(index_individu));
    end
end