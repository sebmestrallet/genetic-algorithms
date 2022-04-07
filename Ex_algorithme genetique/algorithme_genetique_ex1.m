%https://towardsdatascience.com/introduction-to-genetic-algorithms-including-example-code-e396e98d8bf3
%traduction du code Java vers MATLAB

%création de la population
%une ligne pour chaque individu (10), une colonne pour chaque gène (5)
%population = randi(2,10,5)-1;
population = [0 0 1 0 0;
              0 1 0 0 1;
              0 0 0 1 0;
              1 1 0 1 0;
              0 1 1 0 1;
              1 0 0 1 1;
              0 0 1 1 0;
              0 0 1 0 1;
              0 0 0 0 0;
              0 1 1 0 1];

%compteur de générations
num_generation = 0;

%calcul du score d'aptitude de chaque individu
scores = score_aptitude(population);

%meilleur individu
[meilleur_score,meilleur_id] = le_plus_adapte(scores);
meilleur_individu = individu_depuis_population(population,meilleur_id);

%affichage du numéro de génération et du meilleur score
fprintf("Génération: " + num_generation + " Meilleur score: " + meilleur_score+"\n");

%tant que la population n'a pas un individu avec le score maximal
%d'aptitude
while (meilleur_score < 5)
    num_generation = num_generation + 1;
    
    %étape de sélection
    [meilleur_individu,second_meilleur_individu] = selection(population);
    
    %étape de croisement
    [enfant1,enfant2] = croisement(meilleur_individu,second_meilleur_individu);
    
    %étape de mutation (probabiliste)
    if rand()<0.01
        fprintf("Mutation !\n");
        [enfant1,enfant2] = mutation(enfant1,enfant2);
    end
    
    %remplacement du pire individu par un clone du meilleur
    [~,pire_id] = le_moins_adapte(scores);
    population = ajout_meilleur_enfant(population,pire_id,choix_meilleur_enfant(enfant1,enfant2));
    
    %calcul des nouveaux scores d'aptitude
    scores = score_aptitude(population);
    [meilleur_score,meilleur_id] = le_plus_adapte(scores);
    meilleur_individu = individu_depuis_population(population,meilleur_id);
    
    %affichage du numéro de génération et du meilleur score
    fprintf("Génération: " + num_generation + " Meilleur score: " + meilleur_score+"\n");
    
end

fprintf("Solution trouvée à la génération : "+num_generation+"\n");
fprintf("Score d'aptitude : "+meilleur_score+"\n");
fprintf("Gènes : ");
fprintf("%d",meilleur_individu);%affichage de chaque gène 
fprintf("\n");


function [meilleur_individu,second_meilleur_individu] = selection(population)
    scores = score_aptitude(population);
    
    %sélection de l'individu le plus adapté
    [~,meilleur_id] = le_plus_adapte(scores);
    meilleur_individu = individu_depuis_population(population,meilleur_id);

    %sélection du second individu le plus adapté
    [~,second_meilleur_id] = second_le_plus_adapte(scores);
    second_meilleur_individu = individu_depuis_population(population,second_meilleur_id);
end


function [meilleur_individu,second_meilleur_individu] = croisement(meilleur_individu,second_meilleur_individu)
    %choix d'un site de croisement au hasard
    point_croisement = randi(4);%entre 1 et 4 compris

    temp =  meilleur_individu(1:point_croisement);
    meilleur_individu(1:point_croisement) = second_meilleur_individu(1:point_croisement);
    second_meilleur_individu(1:point_croisement) = temp;
end


function [meilleur_individu,second_meilleur_individu] = mutation(meilleur_individu,second_meilleur_individu)
    %choix d'un site de mutation au hasard
    point_mutation_1 = randi(5);%entre 1 compris et 5 compris
    
    meilleur_individu(point_mutation_1) = double(~meilleur_individu(point_mutation_1));
    
    %choix d'un site de mutation au hasard
    point_mutation_2 = randi(5);%entre 1 compris et 5 compris
    
    second_meilleur_individu(point_mutation_2) = double(~second_meilleur_individu(point_mutation_2));
end

function meilleur_enfant = choix_meilleur_enfant(enfant1,enfant2)
    if score_aptitude_individu(enfant1) > score_aptitude_individu(enfant2)
        meilleur_enfant = enfant1;
    else
        meilleur_enfant = enfant2;
    end
end

function population = ajout_meilleur_enfant(population,pire_id,meilleur_enfant)
    population(pire_id,:) = meilleur_enfant;%remplacement du pire par le meilleur
end

function score = score_aptitude(population)
    score = sum(population,2);%somme par ligne
end

function score = score_aptitude_individu(individu)
    score = sum(individu);
end

function [meilleur_score,meilleur_id] = le_plus_adapte(scores)
    [meilleur_score,meilleur_id] = max(scores);
end

function [second_meilleur_score,deuxieme_meilleur_id] = second_le_plus_adapte(scores)
    [~,meilleur_id] = max(scores);%on cherche le meilleur
    scores(meilleur_id) = min(scores)-1;%on le change en min
    [second_meilleur_score,deuxieme_meilleur_id] = max(scores);%on cherche le nouveau min
end

function [pire_score,pire_id] = le_moins_adapte(scores)
    [pire_score,pire_id] = min(scores);
end

function individu = individu_depuis_population(population,index)
    individu = population(index,:);
end