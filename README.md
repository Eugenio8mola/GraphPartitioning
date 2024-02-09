# Graph Partitioning

This project provides implementations of various graph partitioning algorithms in `C++` for efficient partitioning of large graphs.  
The implemented algorithms are:  
+ `Fiduccia-Mattheyses`,  
+ `multilevel Kernighan-Lin`,  
+ `Greedy Graph Growing Partitioning`,
    
  and their parallel versions:   
  + `Parallel Fiduccia-Mattheyses`,  
  + `Parallel multilevel Kernighan-Lin`,  
  + `Parallel Greedy Graph Growing Partitioning`.
  
Each of these algorithms aims at providing the best partition of a given graph `G` such that:     
1. node weights among partitions are `balanced`,  
2. edge cut size is `minimized`.

## Supported Platforms

The programs in this repository can be run on Linux, Windows, and macOS. Ensure you have a C++ compiler installed on your system. If you're using Windows, consider using MinGW or Visual Studio. On macOS and Linux, the default system compiler should suffice.

## Algorithms Implemented

### `Fiduccia-Mattheyses`
This algorithm is based on a local search approach and aims to balance the sizes of partitions while minimizing the edge cut.

### `Multilevel Kernighan-Lin`
The multilevel version of the Kernighan-Lin algorithm improves the partitioning quality by constructing a hierarchy of coarser graphs.
When the coarsest graph is reached,  it applies the partition and refines the edge cut at multiple levels by applying Kernighan-Lin and uncoarsening the graph till the original input graph is obtained.

### `Greedy Graph Growing Partitioning`
This algorithm iteratively grows a partition by greedily selecting nodes based on certain criteria, aiming to achieve a balanced partitioning with minimal communication overhead.

### Parallel Versions
Parallel versions of the above algorithms are implemented to leverage multicore processors for faster computation.

## Getting Started

To use the programs in this repository, follow these steps:

### Cloning the Repository

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/Eugenio8mola/GraphPartitioning.git

2. Navigate to the project directory:

   ```bash
   cd GraphPartitioning
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

   Replace:  
   `num_nodes` with the desired number of nodes that should be used for graph generation,   
   `max_node_weight` with the maximum node weight that should be used for graph generation,  
   `num_partitions` with the desired number of partitions that should be applied.  

## Acknowledgments
This `README`provides general information.  
Each algorithm has its own `README` inside its folder.  
Before running any of the algorithms implemented, read the `README` of the desired sequential or parallel algorithm that you are willing to execute.  
Further details on execution mode, parameters to pass to run the program and functions used are provided.  

