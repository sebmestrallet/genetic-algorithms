from random import random
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from time import sleep
from IPython.display import clear_output

#decimal to grey code
def dec2gc(dec,N):
    binary = dec
    binary ^= (binary >> 1)#conversion happens here
    #convert to string with bin(), remove '0b', pad with '0's for fixed width
    binary = '{:0>{width}}'.format(bin(binary)[2:], 'b', width=N)
    #separate chars with ' ' and use this separator to get an int array
    return np.fromstring(" ".join(binary),dtype=int,sep=" ")

def gc2dec(gc):
    #create a string from the array
    gc = np.array2string(gc, separator='')[1:-1]
    gc = int(gc,base=2)
    inv = 0
    while(gc):
        inv = inv ^ gc
        gc = gc >> 1
    return inv

fonctionnelle_a_minimiser = lambda x: np.sin(2*np.pi*0.01*x)*(-x)*0.02 - 4

code_binaire = lambda decimal: dec2gc(decimal,8);#8 bits

valeur_decimale = lambda binaire: gc2dec(binaire)

def calcul_score(x):
    return -fonctionnelle_a_minimiser(x)

def individu_depuis_index(population,index):
    return population[index,:]

def calcul_scores_population(population):
    scores = np.zeros((population.shape[0],1))
    for index_individu in range(1,population.shape[0]):
        scores[index_individu] = calcul_score(valeur_decimale(individu_depuis_index(population,index_individu)))
    return scores

def tracer_population(population,titre) -> go.Figure:
    x=np.arange(0,2**8)
    y=fonctionnelle_a_minimiser(x)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x,y=y,mode='lines'))
    x = np.apply_along_axis((lambda x: float(valeur_decimale(x))), axis=1, arr=population)
    y = fonctionnelle_a_minimiser(x)
    fig.add_trace(go.Scatter(x=x,y=y,mode='markers',marker_color='black'))
    fig.update_layout(title_text=titre)
    return fig

def afficher_stats_scores(scores,num_generation):
    print("\nGénération {} ----------------".format(num_generation))
    print("Score moyen :    {}".format(scores.mean()))
    print("Meilleur score : {}".format(scores.max()))
    print("Ecart-type :     {}".format(scores.std()))

def selection(population,scores):
    nb_individus = population.shape[0]
    if(nb_individus%2 != 0):
        print("Le nombre d'individus n'est pas pair")
    nb_enfants = int(nb_individus/2)
    index = np.argsort(scores,0)#index qui trierais les scores
    index = np.flip(index)
    population_triee = population[-index,:]#l'appliquer sur population
    population_triee = population_triee.squeeze()#pourquoi cela est necessaire ???????
    individus_selectionnes = population_triee[:nb_enfants]#on ne garde que les nb_enfants premiers
    return individus_selectionnes

def melange_parents(parents):
    return parents[np.random.permutation(parents.shape[0])]

def croisement(parent1,parent2):
    enfant1 = np.copy(parent1)
    enfant2 = np.copy(parent2)
    #tirage au sort d'un point de croisement
    #chromosomes de 8 genes [0:7]
    # |0|1|2|3|4|5|6|7|
    #   ^ ^ ^ ^ ^ ^ ^
    #   0 1 2 3 4 5 6
    #-> 7 points de croisement possible [0:6]
    point_croisement = np.random.randint(0,7)#7 is excluded
    temp = np.copy(enfant1[0:point_croisement+1])#copy the first chunck of enfant1 
    enfant1[0:point_croisement+1] = np.copy(enfant2[0:point_croisement+1])#replace the first chunck of enfant1 by the first chunck of enfant2
    enfant2[0:point_croisement+1] = np.copy(temp)#replace the first chunck of enfant2 by the saved first chunck of enfant1
    return enfant1,enfant2,point_croisement

def index_des_meilleurs(scores,quantite):
    index_originaux = np.arange(scores.shape[0])
    index = np.argsort(scores,0)#index qui trierais les scores
    population_triee = index_originaux[index[::-1]]#l'appliquer sur population
    return population_triee[:quantite]#on ne garde que les quantite premiers

def index_des_pires(scores,quantite):
    index_originaux = np.arange(scores.shape[0])
    index = np.argsort(scores,0)#index qui trierais les scores
    population_triee = index_originaux[index[::1]]#l'appliquer sur population
    return population_triee[:quantite]#on ne garde que les quantite premiers

def afficher_genes_individu(individu):
    print(np.array2string(individu, separator='')[1:-1],end='')

def afficher_population_et_scores(population,scores):
    for index_individu in range(population.shape[0]):
        print('{0: >2} | '.format(index_individu),end='')#fixed width of 2 chars
        afficher_genes_individu(individu_depuis_index(population,index_individu))
        print(' | {}'.format(scores[index_individu]))

def boucle_optimisation(population,proba_croisement,proba_mutation,max_generations,nb_individus_remplaces,afficher_etapes,afficher_evolution_du_score,scores_moyens,meilleurs_scores):
    num_generation = 0
    #calcul des scores
    scores = calcul_scores_population(population)
    if(afficher_etapes):
        print("Population initiale")
        afficher_population_et_scores(population,scores)
    scores_moyens[0,num_generation] = np.mean(scores)
    meilleurs_scores[0,num_generation] = np.max(scores)

    #critère de convergence : population homogène
    #écart-type faible
    while (np.std(scores) > 0.01) and (num_generation < max_generations):
        clear_output(wait=True)#lazy way to "overwrite" the output each generation but makes the output blink

        num_generation = num_generation+1
        
        #---------- SELECTION --------
        #selection de taille_population/2 individus selon leur score
        parents = selection(population,scores)
        if afficher_etapes:
            print("\nParents sélectionnés ({} meilleurs)".format(population.shape[0]/2))
            afficher_population_et_scores(parents,calcul_scores_population(parents))
        
        #---------- CROISEMENTS --------
        #pour chaque paire de parents
        if parents.shape[0]%2 != 0:
            print("Erreur : Le nombre de parents n'est pas pair")
            #TODO : à la place, juste cloner l'individu orphelin
        parents_melanges = melange_parents(parents)
        if afficher_etapes:
            print("\nParents mélangés")
            afficher_population_et_scores(parents_melanges,calcul_scores_population(parents_melanges))
        enfants = np.zeros((parents_melanges.shape[0],parents_melanges.shape[1]))#meme dimensions que parents_melanges
        for i in np.arange(0,parents_melanges.shape[0],2):#de 2 en 2
            parent1 = parents_melanges[i,:]
            parent2 = parents_melanges[i+1,:]
            if afficher_etapes:
                print("\nPaire des parents {}&{}".format(i,i+1))
                print("Parent 1 : ",end='')
                afficher_genes_individu(parent1)
                print("\nParent 2 : ",end='')
                afficher_genes_individu(parent2)
                print()
            if random()<proba_croisement:
                #alors croisement.
                enfant1,enfant2,point_croisement = croisement(parent1,parent2)
                #ecriture dans le tableau enfants
                enfants[i,:] = enfant1
                enfants[i+1,:] = enfant2
                if afficher_etapes:
                    print("Croisement au point {}".format(point_croisement))
                    print("Enfant 1 : ")
                    afficher_genes_individu(enfant1)
                    print("\nEnfant 2 : ")
                    afficher_genes_individu(enfant2)
                    print("")
            else:
                #ecriture dans le tableau enfants
                enfants[i,:] = parent1
                enfants[i+1,:] = parent2
                if afficher_etapes:
                    print("Pas de croisement")
        
        #---------- MUTATIONS --------
        #pour chaque gène
        for i in range(enfants.shape[0]):#pour chaque enfant
            if afficher_etapes:
                print("\nEnfant {}".format(i))
                afficher_genes_individu(individu_depuis_index(enfants,i))
                print("")
            emplacement_mutations = np.zeros((1,8),dtype=np.uint8)
            for j in range(8):#pour chaque gène
                if random()<proba_mutation:
                    emplacement_mutations[0,j] = 1
                    enfants[i,j] = np.uint8(not bool(enfants[i,j]))#complement
            if afficher_etapes:
                afficher_genes_individu(individu_depuis_index(enfants,i))
                print("")
                if np.sum(emplacement_mutations) != 0:#s'il y a eu une mutation
                    for j in range(8):#pour chaque gene
                        if emplacement_mutations[0,j] == 1:
                            print("^",end='')#marqueur de l'emplacement de la mutation
                        else:
                            print(" ",end='')
                    print("")
                else:
                    print("Pas de mutation")
        
        #---------- REMPLACEMENT --------
        #nouvelle population
        #elitisme : les meilleurs enfants sont gardes
        #et remplacent les moins bons individus de la generation d'avant
        enfants = enfants.astype(np.uint8)#pourquoi cela est necessaire ???????
        scores_enfants = calcul_scores_population(enfants)
        afficher_population_et_scores(enfants,scores_enfants)
        index_des_meilleurs_enfants = index_des_meilleurs(scores_enfants,nb_individus_remplaces)
        index_des_pires_parents = index_des_pires(scores,nb_individus_remplaces)
        for i in range(nb_individus_remplaces):
            #remplacement dans la population
            #d'un des moins bons parents par un des meilleurs enfants
            population[index_des_pires_parents[i],:] = enfants[index_des_meilleurs_enfants[i],:]
        if afficher_etapes:
            print("\nPopulation avant selection, croisements et mutations")
            afficher_population_et_scores(population,scores)
            print("index des pires :",end='')
            print(" {}".format(index_des_pires_parents))#espacements entre les index
            print("\nEnfants")
            afficher_population_et_scores(enfants,scores_enfants)
            print("index des meilleurs :",end='')
            print(" {}".format(index_des_meilleurs_enfants))#espacements entre les index
            print("")
        
        #calcul des nouveaux scores
        scores = calcul_scores_population(population)
        
        afficher_stats_scores(scores,num_generation)
        scores_moyens[0,num_generation] = np.mean(scores)
        meilleurs_scores[0,num_generation] = np.max(scores)
        
        if afficher_etapes:
            print("\nNouvelle population")
            afficher_population_et_scores(population,scores)
        
        tracer_population(population,"Population a la generation {}".format(num_generation))
        sleep(0.5)

    print("\nCritere de convergence atteint")
    index_du_max = np.argmax(scores)
    print("Meilleur individu : ",end='')
    meilleur_individu = population[index_du_max,:]
    afficher_genes_individu(meilleur_individu)
    print("\n= {}".format(valeur_decimale(meilleur_individu)))

    if afficher_evolution_du_score:

        #évolution des scores max et moyen au cours des générations
        fig2 = go.Figure()
        x=np.arange(0,num_generation+1)
        fig2.add_trace(go.Scatter(x=x,y=meilleurs_scores[0:num_generation+1][0],mode='lines+markers'))
        fig2.update_layout(title_text="Evolution du score")
        fig2.layout.yaxis.range = [0, 9]
        fig2.show()