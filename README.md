# CS 170 Project Spring 2020

Project spec:

Let *G = (V,E)* be a positive weighted, connected, undirected graph. We would like to find a subgraph *T* of *G* such that:

1. Every vertex is either in *T* or adjacent to a vertex in *T*.

2. *T* is a tree.

3. The average pairwise distance between all vertices in *T* is minimized.

Requirements:

Python 3.6+

You'll only need to install networkx to work with the starter code. For installation instructions, follow: https://networkx.github.io/documentation/stable/install.html

If using pip to download, run `python3 -m pip install networkx`

Files:
- `parse.py`: functions to read/write inputs and outputs
- `solver.py`: where you should be writing your code to solve inputs
- `utils.py`: contains functions to compute cost and validate NetworkX graphs

When writing inputs/outputs:
- Make sure you use the functions `write_input_file` and `write_output_file` provided
- Run the functions `read_input_file` and `read_output_file` to validate your files before submitting!
  - These are the functions run by the autograder to validate submissions
