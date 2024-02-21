#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <math.h>
#include <numeric>
#include <stack>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <chrono>
#include <thread>
#include <unistd.h>

using namespace std;

int k = 8;
int n = 32;
int maxWeight = 6;

atomic<int> totalCutCost(0);
vector<int> Gain;
mutex lck,cmtx,umtx, emtx, adjmtx;

class Graph {

public:
    vector<vector<pair<int, int>>> graph;
    vector<int> nodeWeights;
    vector<int> partition;
    vector<vector<pair<int,int>>> pairList;

};

class GraphHandler {
private:
    Graph g;

public:
    //=================First Constructor========
    GraphHandler() { }

    Graph getGraph() const {
        return g;
    }
    //=================Second Constructor=======
    GraphHandler(int n, int max_w) {

        g.graph.resize(n,vector<pair<int, int>>(n, {0, 0}));
        random_device rd1;
        mt19937 gen1(rd1());
        uniform_int_distribution<int> distribution1(1, max_w);

        for (int i = 0; i < n; i++) {
            int elementsPerLine = rd1() % n + 1;
            for (int j = 0; j < elementsPerLine; j++) {
                if (i == j)
                    g.graph[i][j] = {0, 0};
                else {
                    int weight = distribution1(gen1);
                    g.graph[i][j] = {j+1, weight};
                    g.graph[j][i] = {i+1, weight};
                }
            }
        }

        // Remove pairs with (0, 0)
        for (int i = 0; i < n; i++) {
            g.graph[i].erase(remove_if(g.graph[i].begin(), g.graph[i].end(),
                                       [](const pair<int, int>& p) {
                                           return p.first == 0 && p.second == 0;
                                       }),
                             g.graph[i].end());
        }


        g.nodeWeights.resize(n);

        random_device rd2;
        mt19937 gen2(rd2());
        uniform_int_distribution<int> distribution2(1, max_w);

        for (int i = 0; i < n; i++)
            g.nodeWeights[i] = distribution2(gen2);
    }
    //
    static Graph extract_subgraph(const Graph g, int index_start, int index_end) {
        Graph subgraph;

        // Copy nodes and their weights within the specified range
        for (int i = index_start; i <= index_end; ++i) {
            subgraph.nodeWeights.push_back(g.nodeWeights[i]);

            // Filter graph based on the condition g.graph[i].first < (index_end + 1)
            vector<pair<int, int>> filtered_graph;
            for (const auto& node : g.graph[i]) {
                if (node.first < (index_end + 1)) {
                    filtered_graph.push_back(node);
                }
            }
            subgraph.graph.push_back(filtered_graph);
        }

        return subgraph;
    }


    //print Pairlist for uncoarsening
    static void printPairList(Graph g){
        cout << "PairList of the graph is: "<< endl;
        for (int i = 0; i <  g.pairList.size(); ++i) {
            cout << "Node " << i + 1 << ": ";
            for (const auto& pair : g.pairList[i]) {
                cout << "(" << pair.first << ", " << pair.second << ") ";
            }
            cout << endl;
        }
    }
    //print graph
    static void print(Graph g) {
        cout << "\nNumber of nodes = " << g.graph.size() << "\n";
        cout << "\nAdjacency List:\n";
        for (int i = 0; i < g.graph.size(); i++) {
            cout << i + 1;
            for (const auto &neighbor: g.graph[i]) {
                cout << " (" << neighbor.first << "," << neighbor.second << ")";
            }
            cout << endl;
        }
        cout << "\nNode weigths are: \n" << endl;
        cout << "[";
        for (int i = 0; i < g.graph.size(); i++) {
            cout << g.nodeWeights[i];
            if(i!= (g.graph.size()-1))
                cout << ", ";
        }
        cout << "]\n" << endl;
    }

    //write adjacency list on an output file
    static void saveAdjacencyList(Graph g, const string &filename) {
        ofstream outFile(filename);
        if (!outFile.is_open()) {
            cerr << "Error: Unable to open the file " << filename << endl;
            return;
        }
        for (int i = 0; i < g.graph.size(); i++) {
            outFile << i << " " << g.nodeWeights[i];
            for (const auto& neighbor : g.graph[i]) {
                outFile << " " << neighbor.first << ":" << neighbor.second;
            }
            outFile << endl;
        }
        outFile.close();
    }

    //read adjacency list and save it in a graph
    static void readAdjacencyList(const string &filename, Graph &g) {
        ifstream inFile(filename);
        if (!inFile.is_open()) {
            cerr << "Error: Unable to open the file " << filename << endl;
            return;
        }
        int counter = 0;
        string line;
        while (getline(inFile, line))
            counter++;
        g.graph.resize(counter);
        g.nodeWeights.resize(counter);
        inFile.clear();
        inFile.seekg(0);
        while (getline(inFile, line)) {
            istringstream iss(line);
            int index, weight;
            char ch;
            iss >> index >> weight;
            g.nodeWeights[index] = weight;
            int neighbor, edgeWeight;
            while (iss >> neighbor >> ch >> edgeWeight) {
                g.graph[index].emplace_back(neighbor, edgeWeight);
            }
        }
        inFile.close();
    }
};

void isValidConfiguration(int num_threads, const Graph& myGraph) {
    if (num_threads >= myGraph.graph.size() / 2) {
        std::cerr << "Error: Number of threads must be less than half the size of the graph." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if ((num_threads & (num_threads - 1)) != 0) {
        std::cerr << "Error: Number of threads must be a power of 2." << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void printPartition(const vector<int>& partition) {
    for (size_t i = 0; i < partition.size(); ++i) {
        std::cout << partition[i];
        if (i != partition.size() - 1) {
            cout << " ";
        }
    }
    cout << endl;
}

void computePartitionWeights(const vector<int>& partition, vector<int>& nodeWeights) {

    unordered_map<int, int> partitionWeights;

    for (int i = 0; i < partition.size(); ++i) {
        int currentNode = i + 1;  // start from 1

        // Find the partition to which the current node belongs
        int currentPartition = partition[i];

        // Add the weight of the current node to the corresponding partition's total weight
        partitionWeights[currentPartition] += nodeWeights[i];
    }
    cout <<" "<< endl;
    for (const auto& entry : partitionWeights) {
        cout << "Partition " << entry.first << ": Total Weight = " << entry.second << std::endl;
    }
    cout <<" "<< endl;

}

int computeCutCost(Graph& myGraph, vector<int>& partition) {
    vector<int> cutCost(k, 0);
    int intV = 0;
    int extV = 0;
    vector<int> intVs;
    vector<int> extVs;
    vector<int> myGain;

    for (int i = 0; i < myGraph.graph.size(); ++i) {
        for (const auto& edge : myGraph.graph[i]) {
            int j = edge.first - 1;
            int weight = edge.second;
            if (partition[i] != partition[j]) {
                cutCost[partition[i] - 1] += weight;
                cutCost[partition[j] - 1] += weight;
                extV += weight;
            } else {
                intV += weight;
            }
        }

        intVs.emplace_back(intV);
        extVs.emplace_back(extV);
        int G = extVs[i] - intVs[i];
        myGain.emplace_back(G);
        extV = 0;
        intV = 0;
    }

    for (int& cost : cutCost) {
        cost /= 2;
    }



    int totalCost = accumulate(cutCost.begin(), cutCost.end(), 0);

    // Update totalCutCost using atomic operation
    totalCutCost += totalCost;

    // Update Gain vector
    Gain = myGain;

    //CutCost for each partition
//    cout << "Computed cut cost for each partition is: " << endl;
//    for(size_t t = 0; t < cutCost.size(); t++)
//    {
//        cout << "Partition " << t + 1 << ":  " << cutCost[t]<< endl;
//    }
//
//    cout << "The total cut cost of this partition is: " << totalCutCost << endl;

    return totalCost;
}

void kernighanLin(Graph& myGraph, vector<int>& partition) {
    int nodes = partition.size();
    int n_iterations = nodes / 2;
    vector<bool> locked(nodes, false);

    bool improvement = true;
    vector<pair<int, int>> nodePair;
    pair<int, int> maxPositionIndedeces;
    vector<int> Glist;


    int initialCutCost = computeCutCost(myGraph, partition);

    int count = 0;
    while (improvement && count != (n_iterations - 1)) {
        improvement = false;
        for (size_t iter = 0; iter < n_iterations; iter++) {
            count = iter;
            int maxGain = 0;
            for (int i = 0; i < nodes; i++) {
                if (!locked[i]) {
                    for (int j = i + 1; j < nodes; j++) {
                        if (!locked[j] && partition[i] != partition[j]) {
                            bool neighbor = false;
                            int weight;

                            for (const auto& edge : myGraph.graph[i]) {
                                int node_num = edge.first - 1;
                                weight = edge.second;
                                if (node_num == j) {
                                    neighbor = true;
                                }
                            }

                            if (neighbor) {
                                int GainV1_V2 = Gain[i] + Gain[j] - 2 * weight;
                                if (abs(GainV1_V2) > maxGain) {
                                    maxGain = GainV1_V2;
                                    maxPositionIndedeces.first = i;
                                    maxPositionIndedeces.second = j;
                                }
                            } else {
                                int GainV1_V2 = Gain[i] + Gain[j];
                                if (abs(GainV1_V2) > maxGain) {
                                    maxGain = GainV1_V2;
                                    maxPositionIndedeces.first = i;
                                    maxPositionIndedeces.second = j;
                                }
                            }
                        }
                    }
                }
            }


            Glist.emplace_back(maxGain);

            swap(partition[maxPositionIndedeces.first], partition[maxPositionIndedeces.second]);
            locked[maxPositionIndedeces.first] = true;
            locked[maxPositionIndedeces.second] = true;
            nodePair.emplace_back(maxPositionIndedeces.first, maxPositionIndedeces.second);
            improvement = true;

            int newCutCost = computeCutCost(myGraph, partition);

            if (newCutCost < initialCutCost) {
                initialCutCost = newCutCost;
                improvement = true;
            } else {
                Glist.pop_back();
                swap(partition[maxPositionIndedeces.first], partition[maxPositionIndedeces.second]);
                locked[maxPositionIndedeces.first] = false;
                locked[maxPositionIndedeces.second] = false;
                nodePair.pop_back();
                improvement = false;
            }
        }
    }


    //    cout <<"*---------------------------------------------------------------------------------------------------------------------------------------------------------*" << endl;
//    cout << "Partition after KL:  ";
//    printPartition(partition);
//    cout <<"*---------------------------------------------------------------------------------------------------------------------------------------------------------*" << endl;
    int finalCutCost = computeCutCost(myGraph, partition);
    cout <<"CutSize after KL execution is: " << finalCutCost << endl;
}

void isPowerOf2(int k, const Graph& graph) {
    // Check if k is a power of 2
    bool isKPowerOf2 = (k > 0) && ((k & (k - 1)) == 0);

    // Check if graph size is a power of 2
    bool isGraphSizePowerOf2 = (graph.graph.size() > 0) && ((graph.graph.size() & (graph.graph.size() - 1)) == 0);

    bool isKlessOrEqThanSize = k <= graph.graph.size();
    if (isKPowerOf2 && isGraphSizePowerOf2 && isKlessOrEqThanSize) {
        return;
    } else {
        cerr << "Error: ";
        if (!isKPowerOf2) {
            cerr << "k is not a power of 2. ";
        }
        if (!isGraphSizePowerOf2) {
            cerr << "Graph size is not a power of 2.";
        }
        if (!isKlessOrEqThanSize) {
            cerr << "K is not lower than Graph size.";
        }
        cerr << endl;

        exit(EXIT_FAILURE);
    }
}

void isContiguous(const Graph& graph) {
    if (graph.graph.empty()) {
        cerr << "Error: Graph is empty." << endl;
        exit(EXIT_FAILURE);
    }

    int n = graph.graph.size();
    unordered_set<int> visited;

    stack<int> dfsStack;
    dfsStack.push(0);  // Start DFS from node 0

    while (!dfsStack.empty()) {
        int current = dfsStack.top();
        dfsStack.pop();

        if (visited.count(current) == 0) {
            visited.insert(current);
            for (const auto& neighbor : graph.graph[current]) {
                dfsStack.push(neighbor.first - 1);  // Adjust for 1-based indexing
            }
        }
    }

    if ( visited.size() == n) {
        cout << "The graph is contiguous." << endl;
    } else {
        cerr << "Error: The graph is not contiguous." << endl;
        exit(EXIT_FAILURE);
    }

}

int findNodeWithFewestAdjacents(const std::vector<std::vector<pair<int,int>>>& adjacencyList, vector<bool>& locked, const int index_start, const int index_end) {

    if (adjacencyList.empty()) {
        // Handle the case where the adjacency list is empty
        cout << "List is empty" << endl;
        return -1;
    }

    bool allTrue = std::all_of(locked.begin() + index_start, locked.begin() + index_end + 1, [](bool value) {
        return value;
    });

    if(allTrue)
        return -1; // stop if there are no more nodes

    int minAdjacentNodes = std::numeric_limits<int>::max();
    int nodeWithFewestAdjacents;

    for (int i = 0; i < adjacencyList.size(); ++i) {
        if(locked[index_start + i] == true)
            continue;
        if (adjacencyList[i].size() <= minAdjacentNodes) {
            minAdjacentNodes = adjacencyList[i].size();
            nodeWithFewestAdjacents = i;
        }

    }

    nodeWithFewestAdjacents = nodeWithFewestAdjacents + index_start;
    locked[nodeWithFewestAdjacents] = true;
    return nodeWithFewestAdjacents;

}

pair<int,int> returnBestNode(const Graph& myGraph, const Graph& subgraph,int nodeWithFewestAdjacents, int treshold, vector<int> maxInt,vector<int> minInt, vector<bool> locked, const int index_start, const int index_end) {



    int nearTreshold = std::numeric_limits<int>::max();
    int bestNode = -1;
    int edgeWeight = 0;
    for (const auto &edge: subgraph.graph[nodeWithFewestAdjacents-index_start]) {

        int nodeWithFewestAdjacentsWeight = myGraph.nodeWeights[nodeWithFewestAdjacents];
        int node_num = edge.first - 1;
        if(locked[node_num] == true)
            continue;
        int node_numWeight = myGraph.nodeWeights[node_num];
        int totalNodeWeight = nodeWithFewestAdjacentsWeight + node_numWeight;

        auto it = std::find(maxInt.begin(), maxInt.end(), totalNodeWeight);
        if (it != maxInt.end()) {
            if (totalNodeWeight < nearTreshold) {
                nearTreshold = totalNodeWeight;
                bestNode = node_num;
                edgeWeight = edge.second;
            }


        } else {
            if (totalNodeWeight < treshold) {
                auto it = std::find(minInt.begin(), minInt.end(), totalNodeWeight);
                if (it != minInt.end()) {

                    if (totalNodeWeight < nearTreshold) {
                        nearTreshold = totalNodeWeight;
                        bestNode = node_num;
                        edgeWeight = edge.second;
                    }
                }
            }
        }
    }
    return {bestNode,edgeWeight};
}

void removeElementsFromAdjacency(int i, vector<int> keepPositionOfL, vector<vector<bool>>&marked, vector<vector<pair<int,int>>>& NewAdjacencyList)
{
    for(size_t index = 0; index < keepPositionOfL.size(); index++)
    {
        int value = keepPositionOfL[index];
        if(index!=0)
            value--;
        NewAdjacencyList[i].erase(NewAdjacencyList[i].begin() + value);
        marked[i].erase(marked[i].begin() + value);

    }

}

void updateAdjacencyList(Graph &myGraph,  vector<vector<pair<int,int>>>& pairList, vector<vector<pair<int,int>>>& NewAdjacencyList, vector<int>& allNodes,vector<int>& newNodeWeights,int& sizeCounter) {


    lck.lock();
    Graph copy = myGraph;
    pairList.resize(myGraph.graph.size() / 2);
    NewAdjacencyList.resize(myGraph.graph.size() / 2);
    int nodeWithFewestAdjacents;
    int z;

    for (int index = 0; index < allNodes.size() ; index+=2) {

        nodeWithFewestAdjacents = allNodes[index];
        z = allNodes[index + 1];


        int indexToWrite = z + 1;
        //build and update new Adjacency List

        for (const auto &pair: copy.graph[nodeWithFewestAdjacents]) {
            int nodeToWrite = pair.first;
            int weightToWrite = pair.second;

            auto iterat = std::find_if(copy.graph[z].begin(), copy.graph[z].end(),
                                       [nodeToWrite](const std::pair<int, int> &p) {
                                           return p.first == (nodeToWrite);
                                       });

            if (iterat != copy.graph[z].end()) {
                //REPLACE THE PAIR with the updated one
                // If present, increase the second value
                weightToWrite += iterat->second;

                copy.graph[z].erase(
                        remove_if(copy.graph[z].begin(), copy.graph[z].end(),
                                  [nodeToWrite](std::pair<int, int> p) {
                                      return p.first == nodeToWrite;
                                  }),
                        copy.graph[z].end());

                copy.graph[z].emplace_back(nodeToWrite, weightToWrite);

            } else {
                // If not present, add it to the list
                copy.graph[z].emplace_back(nodeToWrite, weightToWrite);
            }

        }

        lck.unlock();
        int indexedNodeWithFewest = nodeWithFewestAdjacents + 1;
        // Remove pairs with pair.first equal to ( indexToWrite | indexedNodeWithFewest )
        lck.lock();
        copy.graph[z].erase(
                remove_if(copy.graph[z].begin(), copy.graph[z].end(),
                          [indexToWrite, indexedNodeWithFewest](std::pair<int, int> &p) {
                              return (p.first == indexToWrite | p.first == indexedNodeWithFewest);
                          }),
                copy.graph[z].end());
        lck.unlock();

        //create new list starting from node 0
        for (const auto &pair: copy.graph[z]) {
            lck.lock();
            NewAdjacencyList[sizeCounter].emplace_back(pair);
            lck.unlock();
        }


        //push the pair into the pair List:
        lck.lock();
        pairList[sizeCounter].emplace_back(indexToWrite, indexedNodeWithFewest);
        sizeCounter++;
        lck.unlock();

    }

    lck.lock();
    //UPDATE new adj list
    vector<vector<bool>> marked(NewAdjacencyList.size());
    for (size_t i = 0; i < NewAdjacencyList.size(); ++i) {
        marked[i].resize(NewAdjacencyList[i].size(), false);
    }
    lck.unlock();

    lck.lock();
    bool update = false;

    for (int j = 0; j < pairList.size(); j++) {
        int totWeight = 0;

        for (const auto &p: pairList[j]) {

            vector<int> keepPositionOfL;
            for (int i = 0; i < NewAdjacencyList.size(); i++) {
                for (int l = 0; l < NewAdjacencyList[i].size(); l++ ) {
                    auto pair = NewAdjacencyList[i][l];
                    int value1 = pair.first;
                    int value2 = pair.second;

                    if ((pair.first == p.first | pair.first == p.second) && marked[i][l] != true) {
                        totWeight += pair.second;
                        update = true;
                        keepPositionOfL.push_back(l);
                    }
                }
                //END OF LINE
                if (update) {

                    //remove all the pairs
                    removeElementsFromAdjacency( i, keepPositionOfL, marked, NewAdjacencyList);
                    //ADD THE UPDATED PAIR
                    marked[i].push_back(true);
                    NewAdjacencyList[i].emplace_back(j+1, totWeight);
                    totWeight = 0;
                    update = false;
                    keepPositionOfL.clear();

                }

            }

        }
    }

    //UPDATE adjList & weights
    myGraph.graph = NewAdjacencyList;
    myGraph.nodeWeights = newNodeWeights;
    lck.unlock();

    //Print New Coarsed graph
    GraphHandler::print(myGraph);

    myGraph.pairList = pairList;
    GraphHandler::printPairList(myGraph);
    cout << "\n" << endl;

}

void coarsenNodes(Graph &myGraph, vector<int>& partition, int &partitionSize, int avgNodeWeight, vector<int> maxInterval,vector<int> minInterval, vector<bool>&locked,vector<int>& allNodes, vector<int>& newNodeWeights, int sizeCounter, vector<vector<pair<int,int>>> NewAdjacencyList,vector<vector<pair<int,int>>> pairList,  int thread_id, int num_threads)
{


    int index_start = thread_id*(myGraph.graph.size()/num_threads);
    int index_end = (myGraph.graph.size()/num_threads + thread_id*(myGraph.graph.size()/num_threads));
    int counter = index_start;
    vector<int> localNodes;
    Graph subgraph = GraphHandler::extract_subgraph(myGraph, index_start, index_end - 1);
    while(counter < index_end) {
        int nodeWithFewestAdjacents = findNodeWithFewestAdjacents(subgraph.graph, locked, index_start, index_end);
        adjmtx.lock();
        //cout << "Node With Fewest : " << nodeWithFewestAdjacents << "  thread_id : " << thread_id << endl;
        adjmtx.unlock();
        if (nodeWithFewestAdjacents == -1) {
            break;
        }

        localNodes.emplace_back(nodeWithFewestAdjacents);
        adjmtx.lock();
        int nodeWeight =0;
        nodeWeight+= myGraph.nodeWeights[nodeWithFewestAdjacents];
        adjmtx.unlock();

        adjmtx.lock();
        newNodeWeights.emplace_back(nodeWeight);
        adjmtx.unlock();

        partitionSize--;
        counter++;
    }
    adjmtx.lock();
    allNodes.insert(allNodes.end(), localNodes.begin(), localNodes.end());
    adjmtx.unlock();
}

Graph coarseGraph(Graph &myGraph,vector<int>& partition, int &partitionSize, int num_threads) {


    Graph copyGraph = myGraph;
    vector<int> allNodes;
    int maxNode = *std::max_element(myGraph.nodeWeights.begin(), myGraph.nodeWeights.end());
    int minNode = *std::min_element(myGraph.nodeWeights.begin(), myGraph.nodeWeights.end());

    vector<int> maxInterval;
    vector<int> minInterval;
    vector<bool> locked(myGraph.graph.size(), false);

    vector<int> newNodeWeights;
    vector<vector<pair<int,int>>> NewAdjacencyList;
    vector<vector<pair<int,int>>> pairList;

    int sizeCounter = 0;
    NewAdjacencyList.resize(partitionSize);

    int totNodeWeight = std::accumulate(myGraph.nodeWeights.begin(), myGraph.nodeWeights.end(), 0);
    int avgNodeWeight = ceil((2.0*totNodeWeight/myGraph.graph.size())); // 8

    emtx.lock();
    for(int i = avgNodeWeight; i <= (2*maxNode); i++) {
        maxInterval.emplace_back(i); // 8 9 10
    }

    for(int i = (avgNodeWeight - 2*minNode - 1); i < (avgNodeWeight); i++) {
        minInterval.emplace_back(i); //5 6 7
    }
    emtx.unlock();

    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++) {
        threads.emplace_back(coarsenNodes, std::ref(myGraph), std::ref(partition), std::ref(partitionSize), std::ref(avgNodeWeight),  std::ref(maxInterval),std::ref(minInterval), std::ref(locked),std::ref(allNodes), std::ref(newNodeWeights), std::ref(sizeCounter), std::ref(NewAdjacencyList),std::ref(pairList), i, num_threads);
    }

    for (auto& thread : threads) {
        thread.join();
    }

//    cout <<"AllNodes: "<< endl;
//    printPartition(allNodes);


    updateAdjacencyList(myGraph, pairList, NewAdjacencyList, allNodes, newNodeWeights, sizeCounter);


    partitionSize = newNodeWeights.size();
    copyGraph.pairList = pairList;

    return copyGraph;
}

void uncoarseGraph(Graph copyGraph, Graph multilevelGraph, vector<int>& partition, int &partitionSize, int n)
{

    //Print uncoarsed graph
    //GraphHandler::print(copyGraph);

    //GraphHandler::printPairList(copyGraph);

    int size = copyGraph.graph.size();
    vector<int> newPartition(size);


    //reconstructs the partition
    int numberOfPartition = 0;
    for(const auto& vector :copyGraph.pairList)
    {
        for(const auto& pair : vector)
        {
            int index1 = pair.first - 1;
            int index2 = pair.second - 1;

            newPartition[index1] = partition[numberOfPartition];
            newPartition[index2] = partition[numberOfPartition];

        }
        numberOfPartition++;
    }

    partition.resize(size);
    partition = newPartition;

    computePartitionWeights(partition,copyGraph.nodeWeights);
    int cutSize = computeCutCost(copyGraph,partition);
    cout <<"CutSize of the uncoarsed graph: " << cutSize<< endl;

}

void assignPartition(vector<int>& partition)
{
    partition.resize(k);
    int partitionNum = 1;
    for(size_t i = 0; i < partition.size(); i++)
        partition[i] = partitionNum++;
}

void multilevelKL(Graph &multilevelGraph, vector<int>& partition, int &partitionSize, int n, int num_threads) {


    static int unfold = 1;
    static int fold  = 1;
    static auto effective_time = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - chrono::high_resolution_clock::now());
    static auto un_effective_time = effective_time;
    static auto kl_effective_time = effective_time;


    if (n == k) {
        // Base case: Stop recursion when there are only two nodes
        std::cout << "\nCoarsest Graph reached.\nInitial partition is applied.\nNumber of nodes: " << n << ", Number of partitions: " << k << "\n";
        assignPartition(partition);
        printPartition(partition);
        return;
    }else
    {

        auto coarse_start_time = chrono::high_resolution_clock::now();
        cout <<"*---------------------------------------------------------------------------------------------------------------------------------------------------------*" << endl;
        cout << "running coarseGraph() pass = " << fold << endl;
        cmtx.lock();
        Graph copyGraph = coarseGraph(multilevelGraph, partition, partitionSize,num_threads);
        cmtx.unlock();
        auto coarse_end_time = chrono::high_resolution_clock::now();
        auto coarse_duration = chrono::duration_cast<chrono::microseconds>(coarse_end_time - coarse_start_time);
        effective_time = chrono::duration_cast<chrono::microseconds>( effective_time + coarse_duration);
        //cout << "Coarsening Time: " << effective_time.count() << " microseconds" << endl;
        fold++;
        partitionSize = partitionSize/2;
        multilevelKL(multilevelGraph, partition, partitionSize, n/2, num_threads);
        cout <<"*---------------------------------------------------------------------------------------------------------------------------------------------------------*" << endl;
        cout << "running uncoarseGraph() pass = " << unfold<< endl;
        unfold++;
        auto uncoarse_start_time = chrono::high_resolution_clock::now();
        umtx.lock();
        uncoarseGraph(copyGraph, multilevelGraph, partition, partitionSize, n);
        //Uncomment to print partitition
        //cout << "\nPartition of the uncoarsed graph: ";
        //printPartition(partition);
        umtx.unlock();
        auto uncoarse_end_time = chrono::high_resolution_clock::now();
        auto uncoarse_duration = chrono::duration_cast<chrono::microseconds>(uncoarse_end_time - uncoarse_start_time);
        un_effective_time = chrono::duration_cast<chrono::microseconds>( un_effective_time + uncoarse_duration);
        //cout << "Uncoarsening Time: " << un_effective_time.count() << " microseconds" << endl;
        auto kl_start_time = chrono::high_resolution_clock::now();
        kernighanLin(copyGraph, partition);
        auto kl_end_time = chrono::high_resolution_clock::now();
        auto kl_duration = chrono::duration_cast<chrono::microseconds>(kl_end_time - kl_start_time);
        kl_effective_time = chrono::duration_cast<chrono::microseconds>( kl_effective_time + kl_duration);
        //cout << "Kernighan-Lin Time: " << kl_effective_time.count() << " microseconds" << endl;

    }

}

int main()  {

    /*int argc, char *argv[]*/
    /*
    // Check command-line parameters
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <num_nodes>  <max_node_weight> <num_partitions> <num_threads>" << std::endl;
        return EXIT_FAILURE;;
    }

    n = std::stoi(argv[1]);
    maxWeight = std::stoi(argv[2]);
    k = std::stoi(argv[3]);
    int num_threads = std::stoi(argv[4]);
     */

    int num_threads = 4;
    auto start_time = chrono::high_resolution_clock::now();
    vector<int> partition(n, -1);
    //generate Graph
    auto generation_start_time = chrono::high_resolution_clock::now();
    GraphHandler myGraphHandler(n, maxWeight);
    auto generation_end_time = chrono::high_resolution_clock::now();
    auto gen_duration = chrono::duration_cast<chrono::microseconds>(generation_end_time - generation_start_time);
    //cout << "Generation Time: " << gen_duration.count() << " microseconds" << endl;
    // Save the graph
    GraphHandler::saveAdjacencyList(myGraphHandler.getGraph(), "graph.txt");
    // Read the graph from file
    Graph myGraph;
    GraphHandler::readAdjacencyList("graph.txt", myGraph);
    cout << "\n\nGraph Generated."<< endl;

    //Uncomment to print the graph
    GraphHandler::print(myGraph);
    isValidConfiguration(num_threads, myGraph);
    isPowerOf2( k, myGraph);
    isContiguous(myGraph);
    int partitionSize = myGraph.graph.size()/2;
    multilevelKL(myGraph, partition, partitionSize, n, num_threads);
    cout <<"Final partition is: "<< endl;
    printPartition(partition);
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    //cout << "Execution Time: " << duration.count() << " microseconds" << endl;
    return 0;
}
