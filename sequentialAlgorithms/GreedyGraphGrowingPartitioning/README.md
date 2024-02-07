
## Greedy Graph Growing Partitioning Algorithm  

<p style="text-align: justify">The GGP algorithm consists of iteratively generating a set E that gathers half the
vertices of the graph G in terms of weight. This set E is initialized by randomly selecting
a vertex in the graph. During the execution of the algorithm, the vertices of the graph are
divided into three sets: the set E, the set of the vertices adjacent to E, called border of E
and denoted by F, and the set of the remaining vertices, denoted by R. At each iteration,
the vertices adjacent to the vertices of the set E are added to E, that is, E = E ∪F. The
process ends when E contains a set of vertices which represents half of the total weight
of the vertices of the graph. The algorithm then returns the bisection of the graph G
whose first part is E.</p>
<p style="text-align: justify">If this method is very easy to implement, the quality of the partition obtained
depends mainly on the starting vertex chosen. As we have seen, in the case of a
multilevel method, the number of vertices of the coarsened graph to partition is very
low, often less than 200. Thus, this algorithm can be executed several times – around
ten times for example – with a starting vertex which is always different. The bisection
of the lower cut will therefore be chosen as a partition for the refinement phase of the
multilevel method.</p>
<p style="text-align: justify">The GGGP algorithm starts with the initialization of the set E. E is initialized as a set that only contains a vertex randomly
selected in V. The sets R and F can then be initialized. Each iteration of the algorithm
consists of inserting in E the vertex of F that can reduce the most cut. This vertex,
denoted by v, is the vertex of maximal gain of F. Then each vertex of R adjacent to v is
moved to F and its gain is calculated. Similarly, each vertex of F adjacent to v has its
gain recalculated and the next iteration can therefore start. The algorithm ends when
the total weight of E is equal to half the weight of V.</p>


## Supported Platforms

The programs in this repository can be run on Linux, Windows, and macOS. Ensure you have a C++ compiler installed on your system. If you're using Windows, consider using MinGW or Visual Studio. On macOS and Linux, the default system compiler should suffice.

## Getting Started

To use the programs in this repository, follow these steps:

### Cloning the Repository

1. Clone this repository to your local machine:

   ```bash
   git clone https://github.com/Eugenio8mola/GraphPartitioning/sequentialAlgorithms/GreedyGraphGrowingPartitioning.git

2. Navigate to the project directory:

   ```bash
   cd GraphPartitioning/sequentialAlgorithms/GreedyGraphGrowingPartitioning
   ```

### Compiling and Running

1. Compile the program using your preferred C++ compiler:

   ```bash
   g++ -o gggp main.cpp
   ```

2. Run the program:

   ```bash
   ./programName num_nodes max_node_weight num_partitions
   ```

   Replace `num_nodes` with the desired number of nodes that should be used for graph generation, `max_node_weight` with the maximum node weight that should be used for graph generation, and `num_partitions` with the desired number of partitions that should be applied.

### Functions

## `graph`

In the class of the graph we have a constructor that receives two parameters:
1. `n` number of nodes
2. `max_w` maximum weight for nodes

This constructor randomly generates a completed graph approximately.

## `saveAdjacencyList`

This function receives three parameters and saves the graph in the format of the "Adjacency List" in the file. 
1. `g` the graph
2. `n` number of nodes
3. `filename` is the name of the file

## `readAdjacencyList`

This function receives two parameters and reads the graph from the file and puts it in a vector of vectors.
1. `filename` is the name of the file
2. `g` An object from the class of graph

## `get_vertex_weight`

This function receives the number of the node `vertex` and then returns the weight of the node.

## `total_vertex_weight`

This function returns the summation of the weights of all nodes.

## `max_vertex_weight`

This function returns the maximum weight of all nodes.

## `get_edge_weight`

This function receives a pair `edge` including two nodes and then returns the weight of edges:
if there exists an edge it returns the weight of the edge, otherwise it returns the 0.

## `get_adjacent_vertices`

This function receives the number of a node and then returns all nodes that are adjacent to it.
The output of the function is a `set<int>`.

## `make_new_graph`

After splitting the graph in each iteration we call this function and pass the partitions `vec` to this function to generate two new graphs.

## `calculate_gain`

This function receives three parameters and calculates the gain in each iteration to select the maximum gain node.
1. `graph` The graph
2. `vertex` A specific node
3. `E` is a set of nodes that are adjacent to that specific node.

## `GGGP`

This function implements "the Greedy Graph Growing Partitioning algorithm (GGGP)". It receives two parameters:
1. `g` The graph
2. `np` The number of partitions

## `print`

This function receives the grap `g` as an input and prints it in the output.

## `print_partitions`

After splitting the graph in each iteration we will have two partitions `result`, so we pass them to this function and it prints them in the output.
