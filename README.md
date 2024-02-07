# Graph Partitioning

This project provides implementations of various graph partitioning algorithms in C++ for efficient partitioning of large graphs. The implemented algorithms include Fiduccia-Mattheyses, multilevel Kernighan-Lin, Greedy Graph Growing Partitioning, and their parallel versions: Parallel Fiduccia-Mattheyses, Parallel multilevel Kernighan-Lin, Parallel Greedy Graph Growing Partitioning.

## Supported Platforms

The programs in this repository can be run on Linux, Windows, and macOS. Ensure you have a C++ compiler installed on your system. If you're using Windows, consider using MinGW or Visual Studio. On macOS and Linux, the default system compiler should suffice.

## Algorithms Implemented

### Fiduccia-Mattheyses
This algorithm is based on a local search approach and aims to balance the sizes of partitions while minimizing the edge cut.

### Multilevel Kernighan-Lin
The multilevel version of the Kernighan-Lin algorithm improves the partitioning quality by constructing a hierarchy of coarser graphs.
When the coarsest graph is reached,  it applies the partition and refines the edge cut at multiple levels by applying Kernighan-Lin and uncoarsening the graph till the original input graph is obtained.

### Greedy Graph Growing Partitioning
This algorithm iteratively grows a partition by greedily selecting nodes based on certain criteria, aiming to achieve a balanced partitioning with minimal communication overhead.

### Parallel Versions
Parallel versions of the above algorithms are implemented to leverage multicore processors for faster computation.

## Getting Started

To use the programs in this repository, follow these steps:

### Cloning the Repository

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/Eugenio8mola/GraphPartitioning.git
