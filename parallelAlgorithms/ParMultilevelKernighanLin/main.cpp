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
#include <future>
#include <unistd.h>

using namespace std;

int k = 2;
int n = 16;
int maxWeight = 6;
sem_t s;

atomic<int> totalCutCost(0);
vector<int> Gain;
mutex mtx,lck,glock,cmtx,umtx, emtx, adjmtx, fewestmtx, bestmtx;

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

    //print Pairlist for uncoarsening
    static void printPairList(Graph g){
        cout << "PairList unfolded for Uncoarsing the graph is: "<< endl;
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
                cout << " (" << neighbor.first << "," << neighbor.second << ") ";
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

void computePartitionWeights(const vector<int>& partition, vector<int>& nodeWeights) {

    unordered_map<int, int> partitionWeights;

    for (int i = 0; i < partition.size(); ++i) {
        int currentNode = i + 1;  // start from 1

        // Find the partition to which the current node belongs
        int currentPartition = partition[i];

        // Add the weight of the current node to the corresponding partition's total weight
        partitionWeights[currentPartition] += nodeWeights[i];
    }

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

    computePartitionWeights(partition, myGraph.nodeWeights);
    computeCutCost(myGraph, partition);
}

void isPowerOf2(int k, const Graph& graph) {
    // Check if k is a power of 2
    bool isKPowerOf2 = (k > 0) && ((k & (k - 1)) == 0);

    // Check if graph size is a power of 2
    bool isGraphSizePowerOf2 = (graph.graph.size() > 0) && ((graph.graph.size() & (graph.graph.size() - 1)) == 0);

    if (isKPowerOf2 && isGraphSizePowerOf2) {
        return;
    } else {
        cerr << "Error: ";
        if (!isKPowerOf2) {
            cerr << "k is not a power of 2. ";
        }
        if (!isGraphSizePowerOf2) {
            cerr << "Graph size is not a power of 2.";
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

int findNodeWithFewestAdjacents(const std::vector<std::vector<pair<int,int>>>& adjacencyList, vector<bool>& locked) {

    if (adjacencyList.empty()) {
        // Handle the case where the adjacency list is empty
        cout << "List is empty" << endl;
        return -1;
    }

    bool allTrue = std::all_of(locked.begin(), locked.end(), [](bool value) {
        return value;
    });

    if(allTrue)
        return -1; // stop if there are no more nodes

    int minAdjacentNodes = std::numeric_limits<int>::max();
    int nodeWithFewestAdjacents;

    for (int i = 0; i < adjacencyList.size(); ++i) {
        if(locked[i] == true)
            continue;
        if (adjacencyList[i].size() <= minAdjacentNodes) {
            minAdjacentNodes = adjacencyList[i].size();
            nodeWithFewestAdjacents = i;
        }
    }

    return nodeWithFewestAdjacents;

}

pair<int,int> returnBestNode(const Graph& myGraph, int nodeWithFewestAdjacents, int treshold, vector<int> maxInt,vector<int> minInt, vector<bool> locked) {


    int nearTreshold = std::numeric_limits<int>::max();
    int bestNode = -1;
    int edgeWeight = 0;
    for (const auto &edge: myGraph.graph[nodeWithFewestAdjacents]) {

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


    cout <<"-----------------------------"<< endl;
    cout <<"START UPDATE ADJ LIST" << endl;
    cout <<"-----------------------------"<< endl;
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
        cout <<"ERROR HERE line 470"<< endl;
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
        cout <<"Passed here line 501"<< endl;

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
cout <<"not passed"<< endl;


        //create new list starting from node 0
        for (const auto &pair: copy.graph[z]) {
            NewAdjacencyList[sizeCounter].emplace_back(pair);
        }


cout <<"Problem line 522"<< endl;
        //push the pair into the pair List:
        //unique_lock<mutex> myLck{mtx};
        lck.lock();
        pairList[sizeCounter].emplace_back(indexToWrite, indexedNodeWithFewest);
        sizeCounter++;
        lck.unlock();
        //myLck.release();
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
    cout <<"update weights and newAdj"<< endl;
    //UPDATE adjList & weights
    myGraph.graph = NewAdjacencyList;
    myGraph.nodeWeights = newNodeWeights;
    lck.unlock();

    int number = 1;
    cout << "Adjacency List of the coarsed Graph is:" << endl;
    for (const auto &vector: NewAdjacencyList) {
        cout << number++ << " ";
        for (const auto &pair: vector) {
            cout << "(" << pair.first << ", " << pair.second << ") ";
        }
        cout << " " << endl;
    }

    myGraph.pairList = pairList;
    GraphHandler::printPairList(myGraph);
    cout << "\n" << endl;

    for (size_t i = 0; i < marked.size(); ++i) {
        std::cout << "Row " << i << ": ";
        for (size_t j = 0; j < marked[i].size(); ++j) {
            std::cout << marked[i][j] << " ";
        }
        std::cout << std::endl;
    }
cout << "END reached"<< endl;
}

void coarsenNodes(Graph &myGraph, vector<int>& partition, int &partitionSize, int avgNodeWeight, vector<int> maxInterval,vector<int> minInterval, vector<bool>&locked,vector<int>& allNodes, vector<int>& newNodeWeights, int sizeCounter, vector<vector<pair<int,int>>> NewAdjacencyList,vector<vector<pair<int,int>>> pairList,  int thread_id, int num_threads) {




    //n_th = 4;
    //n = 16 -> 8  0 -> 2 1
    sem_wait(&s);
    int num_nodesHalved = myGraph.graph.size();
    int counter = thread_id*(num_nodesHalved/num_threads);
        while(counter < (num_nodesHalved/num_threads + thread_id*(num_nodesHalved/num_threads))) {
            int nodeWithFewestAdjacents = findNodeWithFewestAdjacents(myGraph.graph, locked);
            //cout << "Node With Fewest : " << nodeWithFewestAdjacents << endl;
            if (nodeWithFewestAdjacents == -1)
                break;

            int weight;
            int nodeWeight = myGraph.nodeWeights[nodeWithFewestAdjacents];
            int totEdgeWeight = 0;
            int count_node = 1;

            vector<int> indexes;
            if (locked[nodeWithFewestAdjacents] == false) {
                indexes.emplace_back(nodeWithFewestAdjacents);
                allNodes.emplace_back(nodeWithFewestAdjacents);
                //select the node in the adjacency list whose summed weight is closest to the AVG node weight for the merging
                pair<int, int> p = returnBestNode(myGraph, nodeWithFewestAdjacents, avgNodeWeight, maxInterval,
                                                  minInterval, locked); // index shift in vector
                int node_num = p.first;
                if (locked[node_num] == true)
                    continue;
                weight = p.second;
                nodeWeight += myGraph.nodeWeights[node_num];
                totEdgeWeight += weight;
                count_node++;



                auto it = std::find(maxInterval.begin(), maxInterval.end(), nodeWeight);
                if (it != maxInterval.end()) {
                    indexes.emplace_back(node_num);
                    allNodes.emplace_back(node_num);

                    for (int index = 0; index < indexes.size(); index++) {
                        locked[indexes[index]] = true;
                        partition[indexes[index]] = partitionSize;

                    }

                    newNodeWeights.emplace_back(nodeWeight);

                } else {
                    if (nodeWeight < avgNodeWeight && count_node != 1) {
                        auto it = std::find(minInterval.begin(), minInterval.end(), nodeWeight);
                        if (it != minInterval.end()) {
                            indexes.emplace_back(node_num);
                            allNodes.emplace_back(node_num);

                            for (int index = 0; index < indexes.size(); index++) {
                                locked[indexes[index]] = true;
                                partition[indexes[index]] = partitionSize;

                            }

                            newNodeWeights.emplace_back(nodeWeight);


                        } else {
                            indexes.emplace_back(node_num);
                        }

                    } else {
                        nodeWeight -= myGraph.nodeWeights[node_num];
                    }
                }

            }
            partitionSize--;
            counter++;
        }
    cout << "allNodes Are: ";
    for (int index = 0; index < allNodes.size() ; index++) {

        cout << allNodes[index] << " ";
    }
    cout << ""<< endl;
    sem_post(&s);
}

Graph coarseGraph(Graph &myGraph,vector<int>& partition, int &partitionSize, int num_threads) {


    Graph copyGraph = myGraph;
    vector<int> allNodes;
    int maxNode = *std::max_element(myGraph.nodeWeights.begin(), myGraph.nodeWeights.end());
    int minNode = *std::min_element(myGraph.nodeWeights.begin(), myGraph.nodeWeights.end());

    vector<int> maxInterval;
    vector<int> minInterval;
    vector<bool> locked(myGraph.graph.size(), false);

    vector<pair<int,int>> PairnewNodeWeights;

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

    updateAdjacencyList(myGraph, pairList, NewAdjacencyList, allNodes, newNodeWeights, sizeCounter);

    cout <<"va anche qui?"<< endl;
    partitionSize = newNodeWeights.size();
    copyGraph.pairList = pairList;
/*
    GraphHandler::print(copyGraph);

    GraphHandler::printPairList(copyGraph);

    cout << "Pairlist of coarsed graph \n"<< endl;

    GraphHandler::print(myGraph);

    GraphHandler::printPairList(myGraph);
*/

    return copyGraph;
}

void uncoarseGraph(Graph copyGraph, Graph multilevelGraph, vector<int>& partition, int &partitionSize, int n)
{

    GraphHandler::print(copyGraph);

    GraphHandler::printPairList(copyGraph);

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
    computeCutCost(copyGraph,partition);

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
         cout << "Coarsening Time: " << effective_time.count() << " microseconds" << endl;
        fold++;
        partitionSize = partitionSize/2;
        multilevelKL(multilevelGraph, partition, partitionSize, n/2, num_threads);
        cout <<"*---------------------------------------------------------------------------------------------------------------------------------------------------------*" << endl;
        cout << "running uncoarseGraph() pass = " << unfold<< endl;
        unfold++;
        auto uncoarse_start_time = chrono::high_resolution_clock::now();
        umtx.lock();
        uncoarseGraph(copyGraph, multilevelGraph, partition, partitionSize, n);
        umtx.unlock();
        auto uncoarse_end_time = chrono::high_resolution_clock::now();
        auto uncoarse_duration = chrono::duration_cast<chrono::microseconds>(uncoarse_end_time - uncoarse_start_time);
        un_effective_time = chrono::duration_cast<chrono::microseconds>( un_effective_time + uncoarse_duration);
        cout << "Uncoarsening Time: " << un_effective_time.count() << " microseconds" << endl;
        auto kl_start_time = chrono::high_resolution_clock::now();
        kernighanLin(copyGraph, partition);
        auto kl_end_time = chrono::high_resolution_clock::now();
        auto kl_duration = chrono::duration_cast<chrono::microseconds>(kl_end_time - kl_start_time);
        kl_effective_time = chrono::duration_cast<chrono::microseconds>( kl_effective_time + kl_duration);
        cout << "Kernighan-Lin Time: " << kl_effective_time.count() << " microseconds" << endl;

    }

}

int main() {


    //int main(int argc, char *argv[]) {
        // Check for command-line parameter specifying the number of threads
        /*
        if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <num_threads>" << std::endl;
            return 1;
        }


        int num_threads = std::stoi(argv[1]);
         */


    int num_threads = 4;
    sem_init(&s, 0, num_threads);
    auto start_time = chrono::high_resolution_clock::now();
    vector<int> partition(n, -1);
    //generate Graph
    auto generation_start_time = chrono::high_resolution_clock::now();
    //GraphHandler myGraphHandler(n, maxWeight);
    auto generation_end_time = chrono::high_resolution_clock::now();
    auto gen_duration = chrono::duration_cast<chrono::microseconds>(generation_end_time - generation_start_time);
    cout << "Generation Time: " << gen_duration.count() << " microseconds" << endl;
    // Save the graph
    //GraphHandler::saveAdjacencyList(myGraphHandler.getGraph(), "graph.txt");
    // Read the graph from file
    Graph myGraph;
    GraphHandler::readAdjacencyList("graph.txt", myGraph);
    cout << "\n\nGraph Generated."<< endl;

    GraphHandler::print(myGraph);


    isPowerOf2( k, myGraph);
    isContiguous(myGraph);
    int partitionSize = myGraph.graph.size()/2;
    multilevelKL(myGraph, partition, partitionSize, n, num_threads);
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
    cout << "Execution Time: " << duration.count() << " microseconds" << endl;
    return 0;
}

#if 0
//-----------GRAPH model example-------------
//    Graph myGraph;
    myGraph.graph.resize(8);
    myGraph.nodeWeights = {1, 4, 5, 2, 3, 6, 5, 3};

    myGraph.graph[0].emplace_back(5, 1);
    myGraph.graph[0].emplace_back(3, 2);
    myGraph.graph[0].emplace_back(2, 1);
    myGraph.graph[1].emplace_back(1, 1);
    myGraph.graph[1].emplace_back(3, 2);
    myGraph.graph[1].emplace_back(4, 1);
    myGraph.graph[2].emplace_back(5, 3);
    myGraph.graph[2].emplace_back(4, 2);
    myGraph.graph[2].emplace_back(2, 2);
    myGraph.graph[2].emplace_back(1, 2);
    myGraph.graph[3].emplace_back(2, 1);
    myGraph.graph[3].emplace_back(3, 2);
    myGraph.graph[3].emplace_back(6, 2);
    myGraph.graph[3].emplace_back(7, 5);
    myGraph.graph[4].emplace_back(1, 1);
    myGraph.graph[4].emplace_back(3, 3);
    myGraph.graph[4].emplace_back(6, 2);
    myGraph.graph[5].emplace_back(5, 2);
    myGraph.graph[5].emplace_back(4, 2);
    myGraph.graph[5].emplace_back(7, 6);
    myGraph.graph[6].emplace_back(6, 6);
    myGraph.graph[6].emplace_back(4, 5);
    myGraph.graph[6].emplace_back(8, 3);
    myGraph.graph[7].emplace_back(7, 3);
#endif
