import marimo

__generated_with = "0.18.4"
app = marimo.App()

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import numpy as np
    import plotly.graph_objects as go
    from IPython.display import clear_output
    from rich.table import Table
    from rich.theme import Theme
    from rich.console import Console
    from rich.panel import Panel
    from shutil import copyfile


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Genetic algorithms: 1D example
    """)
    return


@app.cell
def _():
    # create a Rich Console object
    console = Console( 
        color_system="truecolor", # force coloring while https://github.com/Textualize/rich/pull/3651 isn't merged
        theme=Theme({"repr.number": ""}) # no special style for numbers, affecting chromosome printing
    )
    return (console,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Objective function

    Define the function to minimize : $f(x) = −0.02x \times sin(0.01 x \times 2 \pi) − 4$

    With $x$ an integer between $0$ and $2^8 = 255$ included.

    Within the framework of genetic algorithms, we seek to maximize the fitness score of individuals, regarding their environment &rarr; the score is therefore $-f(x)$.

    $f(x)$ is defined so that the score is always positive.

    *Put functions in separate cells for them to be reusable in the other notebook*
    """)
    return


@app.function
def function_to_minimize(x):
    return -0.02 * x * np.sin(0.01 * x * 2 * np.pi) - 4


@app.function
def fitness_score(x):
    return -function_to_minimize(x)


@app.function
# define how to plot the function
def plot_objective_function(x: np.ndarray, y: np.ndarray) -> go.Figure:  
    fig = go.Figure(
        go.Scatter(x=x, y=y, mode='lines', name='objective function'),
        layout_xaxis_range=[0, 255],
        layout_yaxis_range=[-9, 0]
    )
    fig.update_layout(title_text='Function to minimize')
    return fig


@app.cell
def _():
    # compute its value over the domain
    x = np.arange(0, 2 ** 8) # range [0, 2^8=256[ -> [0, 255]
    y = function_to_minimize(x) # evaluate all values in x

    # plot the function
    _fig = plot_objective_function(x, y)
    _fig.show()

    print(f'The minimum is {np.min(y):0.3f} at x={np.argmin(y)}')
    return x, y


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Binary representation

    In order to use a genetic algorithm, we must be able to represent solutions (= individuals) as chromosomes -> alleles vector.

    The simplest way to do so with integers is to use their binary representation.

    But there is a catch: to enforce adjacent inputs to have adjacent binary codes, we must use the [Gray code](https://en.wikipedia.org/wiki/Gray_code).

    This way, the binary representation of two successive integers will differ in only one bit.

    We need to functions, `dec2gc()` to compute the Gray code of a decimal number, and `gc2dec()` the inverse.

    https://www.geeksforgeeks.org/decimal-equivalent-gray-code-inverse/
    """)
    return


@app.function
def dec2bin(dec, N):
    """
    Decimal to Binary conversion
    """
    # convert to string with bin(), remove '0b', pad with '0's for fixed width
    binary = '{:0>{width}}'.format(bin(dec)[2:], 'b', width=N)
    # separate chars with ' ' and use this separator to get an int array
    return np.fromstring(' '.join(binary), dtype=int, sep=' ')


@app.function
def bin2dec(bin):
    """
    Binary to decimal conversion
    """
    # create a string from the array
    bin = np.array2string(bin, separator='')[1:-1] # remove leading and trailing square brackets
    return int(bin, base=2)


@app.function
def dec2gc(dec, N):
    """
    Decimal to Gray code conversion
    """
    binary = dec
    binary = binary ^ binary >> 1 # conversion happens here
    # convert to string with bin(), remove '0b', pad with '0's for fixed width
    binary = '{:0>{width}}'.format(bin(binary)[2:], 'b', width=N)
    # separate chars with ' ' and use this separator to get an int array
    return np.fromstring(' '.join(binary), dtype=int, sep=' ')


@app.function
def gc2dec(gc):
    """
    Gray code to decimal conversion
    """
    #create a string from the array
    gc = np.array2string(gc, separator='')[1:-1] # remove leading and trailing square brackets
    gc = int(gc, base=2)
    inv = 0
    while gc:
        inv = inv ^ gc
        gc = gc >> 1
    return inv


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Compact string representation of chromosomes
    """)
    return


@app.function
# `chromosome` being a numpy array of int (row of `population` defined below)
def chromosome2str(chromosome: np.ndarray) -> str:
    return np.array2string(chromosome, separator='')[1:-1]


@app.cell
def _():
    chromosome2str(np.array([1,1,0,1,0,0,0,1],dtype=int))
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Test the conversion functions
    """)
    return


@app.cell
def _(console):
    decimal_value_1 = 95
    binary_value_1 = dec2bin(decimal_value_1, 8)
    gray_code_1 = dec2gc(decimal_value_1, 8)
    back_to_decimal = gc2dec(gray_code_1)
    console.print(f'{decimal_value_1} has {chromosome2str(binary_value_1)} as binary value\n   and {chromosome2str(gray_code_1)} as Gray code')

    # test inverse conversion
    assert back_to_decimal == decimal_value_1
    back_to_decimal = bin2dec(binary_value_1)
    assert back_to_decimal == decimal_value_1

    decimal_value_2 = decimal_value_1 + 1
    binary_value_2 = dec2bin(decimal_value_2, 8)
    gray_code_2 = dec2gc(decimal_value_2, 8)
    to_print = f'{decimal_value_2} has '
    for gene in range(8):
        if binary_value_2[gene] != binary_value_1[gene]:
            to_print += f'[b][bright_magenta]{binary_value_2[gene]}[/][/]'
        else:
            to_print = to_print + str(binary_value_2[gene])
    to_print = to_print + ' as binary value\n   and '
    for gene in range(8):
        if gray_code_2[gene] != gray_code_1[gene]:
            to_print += f'[b][bright_magenta]{gray_code_2[gene]}[/][/]'
        else:
            to_print = to_print + str(gray_code_2[gene])
    to_print = to_print + ' as Gray code'
    console.print(to_print)

    # test inverse conversion
    back_to_decimal = gc2dec(gray_code_2)
    assert back_to_decimal == decimal_value_2
    back_to_decimal = bin2dec(binary_value_2)
    assert back_to_decimal == decimal_value_2
    return


@app.cell
def _():
    mo.md(r"""
    95 and 96 differ in **6 bits** in binary representation, but in only **1 bit** in Gray code.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Generation of the initial population

    The population is stored as a $N ~ \text{individuals} \times 8 ~ \text{genes}$ matrix. Here $N=20$.

    Let's populate the initial population with random individuals :
    """)
    return


@app.function
def generate_random_population(rng: np.random.Generator):
    return rng.integers(low=0, high=2, size=(20,8))


@app.cell
def _():
    # Create a random number generator with a specified seed.
    # `np.random.seed(value)` is considered a legacy function,
    # so let's use a `np.random.Generator`
    rng: np.random.Generator = np.random.default_rng(seed=112358)

    population = generate_random_population(rng)
    return population, rng


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Evaluate each individual of the initial population
    """)
    return


@app.function
# define how to evaluate a given population
def evaluate_population(population: np.ndarray) -> np.ndarray:
    scores = np.zeros((population.shape[0], 1))  # a column vector of size N, N = number of individuals
    for idx in range(0, population.shape[0]):
        scores[idx] = fitness_score(gc2dec(population[idx, :]))
    return scores


@app.function
# define how to display the scores
def display_scores(population: np.ndarray, scores: np.ndarray, console: Console):
    assert population.shape[0] == scores.shape[0]
    table = Table(title='Fitness scores')
    table.add_column('Index')
    table.add_column('Chromosome')
    table.add_column('Score')
    for idx in range(0, population.shape[0]):
        table.add_row(
            str(idx),
            chromosome2str(population[idx, :]),
            f'{scores[idx, 0]:0.3f}'
        )
    console.print(table)


@app.function
def compute_stats(values: np.ndarray) -> tuple[float, float, float]:
    return (values.mean(), values.max(), values.std())


@app.function
def print_generation_stats(console, mean: float, max: float, std_dev: float):
    console.print(
        Panel.fit(
            f'mean score = {mean:0.3f}\n' +
            f'best score = {max:0.3f}\n' +
            f' std. dev. = {std_dev:0.3f}',
            title=f'Score stats'
        )
    )


@app.cell
def _(console, population):
    # evaluate the initial population
    _scores = evaluate_population(population)

    # display the scores of the initial population
    display_scores(population, _scores, console)
    mean, max, std_dev = compute_stats(_scores)
    print_generation_stats(console, mean, max, std_dev)
    return max, mean, std_dev


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Plot the population
    """)
    return


@app.function
# define how to plot a given population
def plot_population(x, y, population, generation, score_mean, score_max, score_std_dev) -> go.Figure:
    fig = plot_objective_function(x, y)
    x_population = np.apply_along_axis(lambda x: float(gc2dec(x)), axis=1, arr=population)
    y_population = function_to_minimize(x_population)
    fig.add_trace(go.Scatter(
        x=x_population,
        y=y_population,
        mode='markers',
        marker_color='black',
        marker_size=10,
        name='individuals'
    ))
    fig.layout.update(showlegend=False)
    fig.update_layout(title_text=f'Generation {generation}   -   scores : avg = {score_mean:0.3f}, max = {score_max:0.3f}, sd = {score_std_dev:0.3f}')
    return fig


@app.cell
def _(max, mean, population, std_dev, x, y):
    _fig = plot_population(x, y, population, 0, mean, max, std_dev)
    _fig.show()
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Crossover

    In pairs, the chromosomes of two individuals are recombined according to a crossover point
    """)
    return


@app.function
def individuals_crossover(parent1: np.ndarray, parent2: np.ndarray, crossover_point: int) -> tuple[np.ndarray, np.ndarray]:
    assert(parent1.shape == (8,))
    assert(parent2.shape == (8,))
    # chromosomes of 8 genes
    # |0|1|2|3|4|5|6|7|
    #   ^ ^ ^ ^ ^ ^ ^
    #   0 1 2 3 4 5 6
    # -> 7 crossover points possible
    assert(crossover_point >= 0)
    assert(crossover_point <= 6)
    child1 = np.copy(parent1)
    child2 = np.copy(parent2)
    temp = np.copy(child1[0:crossover_point+1]) # copy the first chunk of child1
    child1[0:crossover_point+1] = np.copy(child2[0:crossover_point+1]) # replace the first chunk of child1 by the first chunk of child2
    child2[0:crossover_point+1] = np.copy(temp) # replace the first chunk of child2 by the saved first chunk of child1
    return child1, child2


@app.function
# define how to visualize a crossover between two chromosomes
def display_crossover(console, parent1: np.ndarray, parent2: np.ndarray, crossover_point: int, child1: np.ndarray, child2: np.ndarray):
    console.print(
        f'parent1 : [bright_red]{chromosome2str(parent1)}[/]\n'
        f'parent2 : [bright_cyan]{chromosome2str(parent2)}[/]\n'
        f' child1 : [bright_cyan]{chromosome2str(child1)[:crossover_point+1]}[/][bright_red]{chromosome2str(child1)[crossover_point+1:]}[/]\n'
        f' child2 : [bright_red]{chromosome2str(child2)[:crossover_point+1]}[/][bright_cyan]{chromosome2str(child2)[crossover_point+1:]}[/]'
    )


@app.cell
def _(console):
    # With predefined parents and crossover point

    parent1 = np.array([1,1,0,1,0,0,0,1], dtype=int)
    parent2 = np.array([0,0,0,0,1,1,0,1], dtype=int)
    crossover_point = 3
    child1, child2 = individuals_crossover(parent1,parent2,crossover_point)
    display_crossover(
        console,
        parent1,
        parent2,
        crossover_point,
        child1,
        child2
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    To apply the crossover on the whole population, we have to
    1. Shuffle the population
    1. Group parents in pairs
    1. For each pair, pick a random number according to a crossover probability<br/>If no crossover for the current pair, copy the parents chromosomes into the children ones
    """)
    return


@app.function
# if `console` is not None, crossovers are printed
def population_crossover(parents: np.ndarray, crossover_probability: float, rng: np.random.Generator, console: Console|None = None) -> np.ndarray:
    # shuffle the parents
    parents = parents[rng.permutation(parents.shape[0])]
    children = np.copy(parents)
    for i in np.arange(0,parents.shape[0],2): # two by two
        parent1 = parents[i,:]
        parent2 = parents[i+1,:]
        if rng.uniform(low=0.0, high=1.0) < crossover_probability:
            # generate a random crossover point
            crossover_point = rng.integers(low=0, high=6, size=1, dtype=int, endpoint=True)[0] 
            children[i,:], children[i+1,:] = individuals_crossover(
                parent1,
                parent2,
                crossover_point
            )
            if console is not None:
                display_crossover(
                    console,
                    parent1,
                    parent2,
                    crossover_point,
                    children[i,:],
                    children[i+1,:]
                )
        else:
            # no crossover, leave parents chromosomes in children[i,:] and children[i+1,:]
            if console is not None:
                # like display_crossover() but no colors
                console.print(
                    f'parent1 : {chromosome2str(parent1)}\n'
                    f'parent2 : {chromosome2str(parent2)}\n'
                    f' child1 : {chromosome2str(children[i,:])}\n'
                    f' child2 : {chromosome2str(children[i+1,:])}\n'
                )
        if console is not None:
            console.print('\n') # some space between pairs
    return children


@app.cell
def _(console, population, rng: np.random.Generator):
    crossover_probability = 0.8
    children = population_crossover(population,crossover_probability,rng,console)
    return (children,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Mutation

    At random, some gene the individuals are altered to express the other allele (bit flip)
    """)
    return


@app.function
# if `console` is not None, mutations are printed
def population_mutation(population: np.ndarray, mutation_probability: float, rng: np.random.Generator, console: Console|None = None):
    to_print = ''
    for idx in range(population.shape[0]):
        if console is not None:
            to_print = to_print + (chromosome2str(population[idx, :]) + ' -> ')
        for gene in range(8):
            if rng.uniform(low=0.0, high=1.0) < mutation_probability:
                population[idx, gene] = int(not bool(population[idx, gene]))  # binary complement
                if console is not None:
                    to_print = to_print + f'[orange1]{population[idx, gene]}[/]'
            elif console is not None:
                to_print = to_print + str(population[idx, gene])
        if console is not None:
            to_print = to_print + '\n'
    if console is not None:
        console.print(to_print)


@app.cell
def _(children, console, rng: np.random.Generator):
    mutation_probability = 0.08
    population_mutation(children, mutation_probability, rng, console)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Selection

    It consist of replacing some of the parents (individuals of the last generation) with some children (individuals go through crossover and mutations), to obtain a new generation.

    Several strategies are possible:
    - Keeping the overall $N$ best individuals
    - Replacing the $n$ worse parents by the $n$ best children
    - Roulette wheel
    - Tournament

    Here we will implement the 2nd one.
    """)
    return


@app.function
def indices_of_the_best(scores: np.ndarray, n: int) -> np.ndarray:
    original_indices = np.arange(scores.shape[0])
    sorting_indices = np.argsort(scores, 0) # the indices sorting the scores
    sorted_original_indices = original_indices[sorting_indices[::-1]] # apply them on the original indices
    return sorted_original_indices[:n] # keep only the n^th first


@app.function
def indices_of_the_worse(scores: np.ndarray, n: int) -> np.ndarray:
    original_indices = np.arange(scores.shape[0])
    sorting_indices = np.argsort(scores, 0) # the indices sorting the scores
    sorted_original_indices = original_indices[sorting_indices[::1]] # apply them on the original indices
    return sorted_original_indices[:n] # keep only the n^th first


@app.function
# if `console` is not None, populations before/after the selection are printed
def selection(parents: np.ndarray, children: np.ndarray, n, console: Console|None = None) -> np.ndarray:
    new_population = np.copy(parents)

    # recompute the scores of the parents, because the crossover shuffled them
    parents_scores = evaluate_population(parents)

    # compute the scores of the children
    children_scores = evaluate_population(children)
    
    indices_of_worse_parents = indices_of_the_worse(parents_scores, n)
    indices_of_best_children = indices_of_the_best(children_scores, n)
    
    if console is not None:
        # based on display_scores()
        
        table = Table(title='Parents')
        table.add_column('Index')
        table.add_column('Chromosome')
        table.add_column('Score')
        for idx in range(0, parents.shape[0]):
            score_str = f'{parents_scores[idx, 0]:0.3f}'
            table.add_row(
                str(idx),
                chromosome2str(parents[idx, :]),
                score_str if idx not in indices_of_worse_parents else '[bright_red]' + score_str + '[/]')
        console.print(table)
        
        table = Table(title='Children')
        table.add_column('Index')
        table.add_column('Chromosome')
        table.add_column('Score')
        for idx in range(0, children.shape[0]):
            score_str = f'{children_scores[idx, 0]:0.3f}'
            table.add_row(
                str(idx),
                chromosome2str(children[idx, :]),
                score_str if idx not in indices_of_best_children else '[bright_green]' + score_str + '[/]')
        console.print(table)

    # actual replacement
    for i in range(n):
        # replacement in the population of one of the worse parent
        # by one of the best child
        new_population[indices_of_worse_parents[i], :] = children[indices_of_best_children[i], :]
    return new_population # TODO also update & return scores?


@app.cell
def _(children, console, population):
    n = 10 # the 10 worse parents will be replaced by the 10 best children
    new_population = selection(population, children, n, console)
    return (new_population,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Evaluate and plot the new population (generation 1)
    """)
    return


@app.cell
def _(console, new_population, x, y):
    _scores = evaluate_population(new_population)
    display_scores(new_population, _scores, console)
    _fig = plot_population(x, y, new_population, 1, *compute_stats(_scores))
    _fig.show()
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
