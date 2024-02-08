# Multilevel Kernighan-Lin

<p align="justify">Description of multiKL 



</p>

## Supported Platforms

The programs in this repository can be run on Linux, Windows, and macOS. Ensure you have a C++ compiler installed on your system.   
If you're using Windows, consider using MinGW or Visual Studio. On macOS and Linux, the default system compiler should suffice.


## Getting Started

To use the programs in this repository, follow these steps:

### Cloning the Repository

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/Eugenio8mola/GraphPartitioning/sequentialAlgorithms/multilevelKernighanLin.git

2. Navigate to the project directory:

   ```bash
   GraphPartitioning/sequentialAlgorithms/multilevelKernighanLin

   ```

### Compiling and Running

1. Compile the program using your preferred C++ compiler:

   ```bash
   g++ -o multilevelKL main.cpp
   ```

2. Run the program:

   ```bash
   ./multilevelKL num_nodes max_node_weight num_partitions
   ```

   Replace:  
   `num_nodes` with the desired number of nodes that should be used for graph generation,  
   `max_node_weight` with the maximum node weight that should be used for graph generation,  
   `num_partitions` with the desired number of partitions that should be applied.

# GraphHandler Class

The `GraphHandler` class provides functionality for handling graphs, including generating random graphs, printing graph information, saving and reading adjacency lists from files.

## Constructors

### First Constructor

## `GraphHandler`
Creates an `empty` GraphHandler object.   

### Second Constructor

## `GraphHandler(int n, int max_w)`

Creates a `GraphHandler` object, generates a random graph with:  
number of nodes `n`    
maximum node weight `max_w`.

## Public Methods

## `getGraph`  
Returns the `Graph` stored in the GraphHandler object.

## `printPairList`
Prints the pair list used for uncoarsing the `Graph g`.

## `print`

Prints the adjacency list and node weights of the `Graph g`

## `saveAdjacencyList`
Writes the `Adjacency List`of the graph on the specified output file.  
Receives three parameters:   
1. `g` the graph  
2. `n` number of nodes  
3. `filename` is the name of the file  

## `readAdjacencyList`

Reads the adjacency list from the specified input file and saves it in in a vector of vectors of the specified `Graph` object.  
Receives two parameters:  
the name of the file `filename`   
the `Graph g` in which the adjacency list should be saved 

## `computePartitionWeights`

The `computePartitionWeights` function calculates the total weight of each partition based on the nodes' weights.   
It takes two parameters:  
a vector of integers `partition` representing the partition number to which nodes belong to.  
a vector of integers `nodeWeights` representing the weights of nodes.  

## `computeCutCost`

The `computeCutCost` function calculates the cut cost for a given partition in a graph, which represents the sum of weights of edges crossing the partition boundary.  
It takes two parameters:  
a Graph data structure `myGraph`  
a partition vector `partition` specifying the current partition assigned to each node.  

## `kernighanLin`

The `kernighanLin` function implements the Kernighan-Lin algorithm, which iteratively improves the quality of the partition by swapping pairs of nodes to reduce the cut cost.  
It takes two parameters:  
a Graph data structure `myGraph`  
a partition vector `partition` specifying the current partition assigned to each node.  

## `isPowerOf2`

The `isPowerOf2` function checks if both the input integer `k` and the size of the graph are powers of 2.  
It utilizes bitwise operations to determine whether a number is a power of 2.  
It takes two parameters:   
the number of partitions to apply `k`.  
a Graph data structure `myGraph` used to obtain the size of the generated graph.


## `isContiguous`

The `isContiguous` function checks whether the generated graph is contiguous (connected) using Depth-First Search (DFS) traversal.  
It ensures that all nodes in the graph are visited, indicating connectivity.  
It takes as parameter a Graph data structure `graph`  

## `findNodeWithFewestAdjacents`

The `findNodeWithFewestAdjacents` function finds the index of a node in the graph that has the fewest number of adjacent   nodes, excluding locked nodes.  
It returns:  
the index of the found node 'nodeWithFewestAdjacents' if a node is found   
'-1' if all nodes are locked or the provided adjacency list is empty.   
It helps in selecting nodes for certain graph-based computations allowing to perform a smarter coarsening.  
It takes two parameters:  
a vector of vector of pairs of integers 'adjacencyList' which is the adjacency list of the graph  
a vector of boolean values 'locked' which is used to keep track of nodes already coarsed and that should not be considered.  

## `returnBestNode`

The `returnBestNode` function finds the best node connected to a node with the fewest adjacents, considering the weight of connected nodes and edges. It returns the index of the best node satisfying certain criteria.

## `removeElementsFromAdjacency`

The `removeElementsFromAdjacency` function selectively removes elements from the adjacency list and the corresponding marked matrix for a specific node, based on specified positions.

## `updateAdjacencyList`

The `updateAdjacencyList` function updates the adjacency list and related data structures of a graph based on specified pairs of nodes and their weights.

## `coarseGraph`

The `coarseGraph` function performs graph coarsening by merging nodes based on certain criteria, updating the graph structure, and returning a coarser-level graph.

## `uncoarseGraph`

The `uncoarseGraph` function reconstructs partitioning information for an uncoarsened graph after a graph coarsening operation.

## `assignPartition`

The `assignPartition` function initializes a partition vector with consecutive partition numbers starting from 1 up to `k`.

## `multilevelKL`

The `multilevelKL` function implements a multilevel graph partitioning algorithm based on the Kernighan-Lin heuristic. It recursively coarsens the graph, applies KL partitioning at each level, and then uncoarsens the graph.

For detailed usage and examples of each function, refer to the respective source code files in this repository.

