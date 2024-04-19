Simple examples of genetic algorithms.
University project, originally with MATLAB, to be transcoded in Python.

There is a difference between the MATLAB script and the Python one :
In the MATALB script, half of the init population is selected, then crossover & mutations are applied on this half of the population.
An the end, a fixed number of the best children replaces the worse parents.
In the Python notebook, a copy of all the init population go through the crossover & mutation steps. At the end, the two set of N individuals are combined into a single set of N individuals, the new generation.

# Python environment

The project is managed with [Rye](https://rye-up.com/) ([`astral-sh/rye`](https://github.com/astral-sh/rye)).

1. Install Rye
1. `rye sync` creates the venv from `pyproject.toml`
1. If you use VSCode to edit/run Jupyter notebooks, select the Python kernel in `.venv`.

# MATLAB code for active contours

**Active Contour MATLAB algorithm from Ritwik Kumar**
Ritwik Kumar (2022). Snakes: Active Contour Models (https://www.mathworks.com/matlabcentral/fileexchange/28109-snakes-active-contour-models), MATLAB Central File Exchange. Retrieved April 7, 2022. 