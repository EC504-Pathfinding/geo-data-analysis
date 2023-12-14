# Geographic Analysis and Pathfinding System

This repository contains code and documentation for our final project in ENG EC 504: Advanced Data Structures at Boston University.  

## Project Team

- Kashyap Panda - kpanda@bu.edu
- Jordan Koseski - jkoseski@bu.edu  
- John Burke - jwburke@bu.edu
- Ritesh Rana - ritesh19@bu.edu

## Project Description  

This project implements two graph search algorithms, Dijkstra's algorithm and the A* search algorithm, to find efficient driving routes between cities in the contiguous United States.

The code reads in a dataset of US cities and constructs a graph representation where edges connect cities within a specified distance threshold. Dijkstra's algorithm iteratively explores the graph to find shortest paths from a start city to all possible destinations. In contrast, A* selectively explores the graph using a heuristic function to find a least-cost path from start to goal without considering all nodes.  

Both algorithms output an optimal path of cities to travel between a given start and end location, providing intermediate stopovers for long trips.

## Usage  

To use this project:

1. Clone the repository
2. Run `make` to compile the C++ code
3. Execute `./GeoDataAnalysis` to run the program 
4. Enter a source and destination city when prompted  
5. View the printed output path and timing comparisons  

Output files containing the generated routes are saved to the `output/` directory for each city pair analyzed.

## Files  

- `src/`: Source code for the project  
- `input/`: Input data files
- `output/`: Output files with generated routes
- `docs/`: Documentation including project report 
- `makefile`: Compilation scripts

## Dependencies  

- C++ compiler
- Makefile

## References  

See the project report for a list of references and data sources.