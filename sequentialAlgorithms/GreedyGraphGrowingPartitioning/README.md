`graph(int n, int max_w)`

In the class of the graph we have a constructor that receives two parameters:
1. `n` number of nodes
2. `max_w` maximum weight for nodes

This constructor randomly generates a completed graph approximately.
#
`saveAdjacencyList(graph g, int n, const string &filename)`

This function receives three parameters and saves the graph in the format of the "Adjacency List" in the file. 
1. the graph
2. number of nodes
3. the name of the file
#
`readAdjacencyList(const string &filename, graph &g)`

This function receives two parameters and reads the graph from the file and puts it in a vector of vectors.
1. The name of the file
2. An object from the class of graph
#
`get_vertex_weight(int vertex)`

This function receives the number of the node and then returns the weight of the node
#
`total_vertex_weight()`

This function returns the summation of the weights of all nodes.
#
`max_vertex_weight()`

This function returns the maximum weight of all nodes.
#
`get_edge_weight(const pair<int, int>& edge)`

this function receives a pair including two nodes and then returns the weight of edges:
if there exists an edge it returns the weight of the edge, otherwise it returns the 0.
#
`get_adjacent_vertices(int vertex)`

This function receives the number of a node and then returns all nodes that are adjacent to it.
The output of the function is a set<int>.
#
`make_new_graph(vector<int>& vec)`

After splitting the graph in each iteration we call this function and pass the partitions to this function to generate two new graphs.
#
`calculate_gain(graph& graph, int vertex, unordered_set<int>& E)`

This function receives three parameters and calculates the gain in each iteration to select the maximum gain node.
1. The graph
2. A specific node
3. a set of nodes that are adjacent to that node.
#
`GGGP(graph& g, int np)`

This function implements "the Greedy Graph Growing Partitioning algorithm (GGGP)".
It receives two parameters:
1. The graph
2. The number of partitions
#
`print(graph g)`

This function receives the grap as an input and prints it in the output.
#
`print_partitions(vector<vector<int>> result)`

After splitting the graph in each iteration we will have two partitions, so we pass them to this function and it prints them in the output.


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

   Replace `num_nodes` with the desired number of nodes that should be used for graph generation, `max_node_weight` with the maximum node weight that should be used for graph generation and `num_partitions` with the desired number of partitions that should be applied.
