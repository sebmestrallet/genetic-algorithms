import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

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

def tracer_population(population,titre):
    x=np.arange(0,2**8)
    y=fonctionnelle_a_minimiser(x)
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x,y=y,mode='lines'))
    x = np.apply_along_axis((lambda x: float(valeur_decimale(x))), axis=1, arr=population)
    y = fonctionnelle_a_minimiser(x)
    fig.add_trace(go.Scatter(x=x,y=y,mode='markers',marker_color='black'))
    fig.update_layout(title_text=titre)
    fig.show()

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
    population_triee = population[index[::-1]]#l'appliquer sur population
    individus_selectionnes = population_triee[:nb_enfants]#on ne garde que les nb_enfants premiers
    return individus_selectionnes

def afficher_genes_individu(individu):
    print(np.array2string(individu, separator='')[1:-1],end='')

def afficher_population_et_scores(population,scores):
    for index_individu in range(population.shape[0]):
        print('{0: >2} | '.format(index_individu),end='')#fixed width of 2 chars
        afficher_genes_individu(individu_depuis_index(population,index_individu))
        print(' | {}'.format(scores[index_individu]))
    