# Multilevel Kernighan-Lin

Description of multiKL

## Supported Platforms

The programs in this repository can be run on Linux, Windows, and macOS. Ensure you have a C++ compiler installed on your system. If you're using Windows, consider using MinGW or Visual Studio. On macOS and Linux, the default system compiler should suffice.


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
   g++ -o partition main.cpp
   ```

2. Run the program:

   ```bash
   ./programName num_nodes max_node_weight num_partitions
   ```

   Replace `num_nodes` with the desired number of nodes that should be used for graph generation, `max_node_weight` with the maximum node weight that should be used for graph generation and `num_partitions` with the desired number of partitions that should be applied.

## `computePartitionWeights(const vector<int>& partition, vector<int>& nodeWeights)`

The `computePartitionWeights` function calculates the total weight of each partition based on the nodes' weights. It takes two parameters: a vector of integers `partition` representing the partition IDs of nodes and a vector of integers `nodeWeights` representing the weights of nodes.

## `computeCutCost`

The `computeCutCost` function calculates the cut cost for a given partition in a graph, which represents the sum of weights of edges crossing the partition boundary. It takes a graph `myGraph` and a partition vector `partition` as input.

## `kernighanLin`

The `kernighanLin` function implements the Kernighan-Lin algorithm, which iteratively improves the quality of the partition by swapping pairs of nodes to reduce the cut cost. It takes a graph `myGraph` and a partition vector `partition` as input.

## `isPowerOf2`

The `isPowerOf2` function checks if both the input integer `k` and the size of the graph are powers of 2. It utilizes bitwise operations to determine whether a number is a power of 2.

## `isContiguous`

The `isContiguous` function checks whether a graph is contiguous (connected) using Depth-First Search (DFS) traversal. It ensures that all nodes in the graph are visited, indicating connectivity.

## `findNodeWithFewestAdjacents`

The `findNodeWithFewestAdjacents` function finds and returns the index of a node in the graph that has the fewest adjacent nodes, excluding locked nodes. It helps in selecting nodes for certain graph-based computations.

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

