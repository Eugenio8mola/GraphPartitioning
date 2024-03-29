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
1. the Graph class `myGraph`  
2. the partition vector `partition` specifying the current partition assigned to each node.  

## `kernighanLin`

The `kernighanLin` function implements the Kernighan-Lin algorithm, which iteratively improves the quality of the partition by swapping pairs of nodes to reduce the cut cost. 

It takes two parameters:  
1. the Graph class `myGraph`  
2. the partition vector `partition` specifying the current partition assigned to each node.  

## `isPowerOf2`

The `isPowerOf2` function checks if both the input integer `k` and the size of the graph in terms of number of nodes `n` are powers of 2.  
It also checks that the number of partition `k` is doesn't exceed the number of nodes `n`.
It utilizes bitwise operations to determine whether `k` and `n` are powers of 2. 
If this condition is not satisfied it prints an `ERROR`.  

It takes two parameters:   
1. the number of partitions to apply `k`.  
2. the Graph data structure `myGraph` used to obtain the size of the generated graph.


## `isContiguous`

The `isContiguous` function checks whether the generated graph is contiguous (connected) using Depth-First Search (DFS) traversal.  
It ensures that all nodes in the graph are visited, indicating connectivity.  
If this condition is not satisfied it prints an `ERROR`.  
It takes as parameter a Graph data structure `graph`.    

## `findNodeWithFewestAdjacents`

The `findNodeWithFewestAdjacents` function finds the index of a node in the graph that has the fewest number of adjacent   nodes, excluding locked nodes. 
It helps in selecting nodes allowing to perform a smarter coarsening.  

It returns:  
+ `nodeWithFewestAdjacents` the index of the found node, if a node is found   
+ `-1` if all nodes are locked or the provided adjacency list is empty.  
   
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
+ the `pair<int, int>` = `(bestNode,edgeWeight)` in which:  
`bestNode` is the index of the best node found   
`edgeWeight` is the edge weight of the best node found.

It takes six parameters:  
1. the Graph class `myGraph`  
2. the node for which the best pairing needs to be determined `nodeWithFewestAdjacents`   
3. the threshold parameter setting a condition for selecting the best pairing node `threshold`, nodes that meet or exceed this threshold are considered suitable candidates for pairing.   
4. the vector `maxInt` that constrains the selection of pairing nodes based on one of the corresponding maximum values specified in this vector used as correction factor to the threshold value.   
5. the vector `minInt` that constrains the selection of pairing nodes based on one of the corresponding minimum values specified in this vector used as correction factor to the threshold value. 
6. the boolean vector `locked` indicating whether each node in the graph is locked or not. 
   Locked nodes are excluded from consideration as potential pairing candidates, ensuring that they remain unchanged during the pairing process.  

## `removeElementsFromAdjacency`

The `removeElementsFromAdjacency` function removes elements from the adjacency list.  
The corresponding element is removed also from the `marked` matrix in the same position.  

It takes four parameters:    
1. the index `i` of the row in the adjacency list from which elements need to be removed.  
2. the vector `keepPositionOfL` contains the shift to apply to the adacency list in order to compute the positions of elements that need to be removed.  
3. the vector of vectors `marked` representing the marked matrix associated with the adjacency list.  It is modified along with the adjacency list.   
4. the vector of vectors `NewAdjacencyList` representing the adjacency list from which elements have to be removed.  It is modified within the function.  

## `updateAdjacencyList`

The `updateAdjacencyList` function is responsible of the updating of the adjacency list of a graph right after it gets coarsed.
It squeezes the original adjacency list into a new smaller one resulting from the merging of two nodes into a `supernode`.  
Each `supernode` together with its index and edge weigth is put inside the new generated adjacency list.  
A double reference structure is used in order to keep track of already coarsed nodes added in the `adjacencyList`.  
The result is achieved by defining the vector of vectors of bolean values `marked` in order to mark nodes already coarsed and tell the program to don't check them again for a possible coarsing. 
The function `removeElementsFromAdjacency` is called to remove a pair from the adjacency list.    
The entry of the `marked` vector of vectors of bolean values associated with the removed pair is also removed due to the correspondence created between the two structures.   
The node pairs removed are all specified in the `vector<vector<pair<int,int>>>` `pairList`.  
Removed pairs are then replaced by the generated `supernode`.  

The function takes six parameters:  
1. the Graph class `myGraph`
2. the vector of vectors `pairList` containing pairs of nodes representing the pairings that need to be considered for updating the adjacency list.  
3. the vector of vectors `NewAdjacencyList` representing the updated adjacency list after the function execution.  
   It is modified within the function.  
4. the vector `allNodes` containing the list of nodes used to apply changes to the adjacency list.
5. the vector `newNodeWeights` storing the weights of nodes in the graph after the adjacency list is updated.  
6. the integer `sizeCounter` tracking the size of the adjacency list after modifications.
   
## `coarseGraph`

The `coarseGraph` function performs graph coarsening by merging nodes based on certain criteria, updating the graph   structure, adjacency list, node weights and returning a coarser graph.  
Computations performed involve multiple calls to the `updateAdjacencyList` function, the `returnBestNode` function and the `findNodeWithFewestAdjacents` function.  
Selected pairs of nodes are merged into `supernodes` and a `pairList` storing the pairs used for the merging process is generated.  

Returns:  
+ the Graph class `copyGraph` which is the graph at the `n-1` coarsening stage that will be used recursively in the 
  uncoarsening process.
  
It receives three parameters:   
1. the input graph to be coarsed `myGraph`.    
2. the vector `partition` storing the partition assignments for each node of the input graph to be coarsed.  
3. the reference to the size of the partition `partitionSize`.  

## `uncoarseGraph`

The `uncoarseGraph` function reconstructs partitioning information for a graph that has been coarsed, computes the total  
weight of a partition and calculates the cut cost for the reconstructed partition.   
The reconstruction is performed by exploiting the pairs contained in `pairList` and assigning to the pairs the same  
partition of the `supernodes` originated by their merging. 

It receives four parameters:  
1. the recursively uncoarsed graph `copyGraph`  
2. the coarsed graph `multilevelGraph`    
3. the vector `partition` storing the partition of the coarsed graph for each node.  
4. the value `partitionSize` referencing the size of the partition.   
5. the number of nodes `n` in the graph.    

## `assignPartition`

The `assignPartition` function assigns the partition to the graph, when the coarsest graph is obtained.    
`k` partitions are assigned, one for each element of the `partition` vector.    
Partition numbers are consecutive, they start from `1` and go up to `k`.    
Receives as parameter the vector `partition`.   

## `multilevelKL`

The `multilevelKL` function acts as a `wrapper` implementing a multilevel graph partitioning algorithm based on the Kernighan-Lin heuristic.  
It recursively coarsens the graph, uncoarsens the graph and applies KL partitioning at each level.  

It receives four paramters:  
1. the input graph to be partitioned `multilevelGraph`.  
2. the vector `partition` storing the initial partition assignments for each node.
3. the value `partitionSize` referencing the size of the partition.  
4. the number of nodes `n` in the graph.  

For detailed usage and implementation details, refer to the respective source code files in this repository.

