# Multilevel Kernighan-Lin
The multilevel Kernighan-Lin algorithm is a graph partitioning technique that adopts an advanced stategy to efficiently divide a graph of `n` nodes into `k` partitions. (with `n` and `k` powers of 2)    
This technique completes the Kernighan-Lin algorithm for graph partitioning, which considers only the edge weights (node weights are considered to be all equal to 1), by considering also the node weights in order to  produce partitions whose total node weigth is balanced and whose edge cut size is minimized.   
It employs a strategy based on a divide and conquer approach.     
The initial graph is coarsed into a smaller one by recursive calls to the `coarseGraph` method.     
Coarsing involves merging multiple vertices into `supernodes` to reduce the graph's size.     
The process is iterated until a sufficiently small graph is obtained, which generally coincides with the size of the graph equal to the number of partitions `k`, which is the stopping condition.  
Once the coarsest graph is achieved, the initial partition is applied.   
Each node `n` is assigned to one of the `k` partitions.  
After obtaining a partition for the coarsest graph, the algorithm gradually unfolds the partitioning to finer levels of the original graph. 
This is done by recursively unfolding the graph calling the `uncoarseGraph` method and refining the partitioning obtained at each level using a graph partition algorithm used for the refinement.  
As the algorithm uncoarsens the graph, it refines the partitioning by unfoilding the merged node pairs instead of the `supernodes` and reconstructing the adjacency list of the graph, ensuring the partitioning is more   accurate.  
The process is iterated recursively until the original graph is fully reconstructed.   
At each uncoarsening step, the refinement process is applied.   
It involves refining the partitioning obtained from the previous coarsening level.     
This refinement can involve local optimization techniques, such as the `Kernighan-Lin algorithm`, to further improve the partition quality or other algorithms that can be used to perform the refinement.    
The application of the Kernighan-Lin algorithm involves iteratively moving vertices between the partitions to improve a predefined objective, such as minimizing edge cuts.   
By employing coarsening and uncoarsening strategies, the multilevel Kernighan-Lin algorithm effectively balances the trade-off between partition quality and computational efficiency, making it suitable for partitioning large graphs efficiently while maintaining partition quality and preserving the original graph structure.  

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
   cd GraphPartitioning/sequentialAlgorithms/multilevelKernighanLin

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
1. the name of the file `filename`   
2. the `Graph g` in which the adjacency list should be saved 

# Implemented methods

## `computePartitionWeights`

The `computePartitionWeights` function calculates the total weight of each partition based on the nodes' weights.   
It takes two parameters:  
1. the vector of integers `partition` representing the partition number to which nodes belong to.  
2. the vector of integers `nodeWeights` representing the weights of nodes.  

## `computeCutCost`

The `computeCutCost` function calculates the cut cost for a given partition in a graph, which represents the sum of weights of edges crossing the partition boundary.  
It takes two parameters:  
1. the Graph data structure `myGraph`  
2. the partition vector `partition` specifying the current partition assigned to each node.  

## `kernighanLin`

The `kernighanLin` function implements the Kernighan-Lin algorithm, which iteratively improves the quality of the partition by swapping pairs of nodes to reduce the cut cost.  
It takes two parameters:  
1. the Graph data structure `myGraph`  
2. the partition vector `partition` specifying the current partition assigned to each node.  

## `isPowerOf2`

The `isPowerOf2` function checks if both the input integer `k` and the size of the graph in terms of number of nodes `n` are powers of 2.  
It utilizes bitwise operations to determine whether a number is a power of 2.  
It takes two parameters:   
1. the number of partitions to apply `k`.  
2. the Graph data structure `myGraph` used to obtain the size of the generated graph.


## `isContiguous`

The `isContiguous` function checks whether the generated graph is contiguous (connected) using Depth-First Search (DFS) traversal.  
It ensures that all nodes in the graph are visited, indicating connectivity.  
It takes as parameter a Graph data structure `graph`  

## `findNodeWithFewestAdjacents`

The `findNodeWithFewestAdjacents` function finds the index of a node in the graph that has the fewest number of adjacent   nodes, excluding locked nodes.  
It returns:  
1.  `nodeWithFewestAdjacents` the index of the found node, if a node is found   
2. `-1` if all nodes are locked or the provided adjacency list is empty.   
It helps in selecting nodes for certain graph-based computations allowing to perform a smarter coarsening.  
It takes two parameters:  
1. the vector of vector of pairs of integers `adjacencyList` which is the adjacency list of the graph  
2. the vector of boolean values `locked` which is used to keep track of nodes already coarsed and that should not be considered.  

## `returnBestNode`

The `returnBestNode` function finds the best node to coarse in order to generate a `supernode`.
The `bestNode` is selected among the ones in the adjacency list of the node with the fewest number of adjacents returned by the function `findNodeWithFewestAdjacents`.  
This function selects appropriately the best matching node considering the sum of the node weights of the nodes to merge and compares it with a `threshold`computed on the entire node weigth set.  
Together with the `threshold` two intervals are also passed.  
This procedure turns the coarsening process into a `smart coarsening`.
The method returns:
the `pair<int, int>` = `(bestNode,edgeWeight)` in which:  
`bestNode` is the index of the best node found   
`edgeWeight` is the edge weight of the best node found.  
It takes six parameters:  
1. the Graph data structure `myGraph`  
2. the node for which the best pairing needs to be determined `nodeWithFewestAdjacents`   
3. the threshold parameter setting a condition for selecting the best pairing node `threshold`, nodes that meet or exceed this threshold are considered suitable candidates for pairing.   
4. the vector `maxInt` that constrains the selection of pairing nodes based on one of the corresponding maximum values specified in this vector used as correction factor to the `threshold` value.   
5. the vector `minInt` that constrains the selection of pairing nodes based on one of the corresponding minimum values specified in this vector used as correction factor to the `threshold` value. 
6. the boolean vector `locked` indicating whether each node in the graph is locked or not. 
   Locked nodes are excluded from consideration as potential pairing candidates, ensuring that they remain unchanged during the pairing process.  

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

