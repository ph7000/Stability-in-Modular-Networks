/**
* Network.cpp
*
* PHYS-302 Fall 2018
* Modular Network Project
* 
* Benjamin Stein-Lubrano
* Professor David Egolf
*
**/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>

using namespace std;

namespace patch
{
	template < typename T > std::string to_string(const T& n)
	{
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

/****************************************************************************
**                                                                         **
**                               constants                                 **
**                                                                         **
*****************************************************************************/

const int NUM_MODULES = 8,
NODES_PER_MODULE = 8,
NUM_NODES = (NUM_MODULES * NODES_PER_MODULE),
NUM_NETWORKS = 100;

const int MAX_DISTANCE = NUM_NODES - 1;

const double
PROB_INT_CONN = (2.0 / 7.0), // (3.2 / 7.0),
PROB_EXT_CONN = (2.0 / 56.0); // (0.8 / 56.0);

/****************************************************************************
**                                                                         **
**                        class Node declaration                           **
**                                                                         **
*****************************************************************************/

class Node
{
public:
	// constructors.
	// Default
	// Convert: population, then population and lyapunov vector contribution
	// Copy
	Node();
	Node(double);
	Node(double, double);
	Node(Node&);

	// accessors and mutators
	double get_startPop() { return startPop; }
	void set_startPop(double pop) { startPop = pop; }
	double get_population() { return population; }
	void set_population(double pop) { population = pop; }
	int get_numIntConns() { return numIntConns; }
	void set_numIntConns(int Conn) { numIntConns = Conn; }
	int get_numExtConns() { return numExtConns; }
	void set_numExtConns(int Conn) { numExtConns = Conn; }

private:
	double startPop;
	double population;
	int numIntConns;
	int numExtConns;
}; // END Node class declaration

/****************************************************************************
**                                                                         **
**                         function declarations                           **
**                                                                         **
*****************************************************************************/

// functions
void condenseConnections();

void showConnectionsBoolean(fstream&);
void showConnectionsVec(fstream&);
void listNumConnections(fstream&);
void showPopulations(fstream&);

void createNodesRand();
void createInternalConnections(int, double);
void createExternalConnections(double);

double modCalc();

void findNodesDistance(int, int,
	vector<int> condConn[NUM_NODES], vector<int> &ultDistNodes);

/****************************************************************************
**                                                                         **
**                                 setup                                   **
**                                                                         **
*****************************************************************************/

// array of Nodes
Node arrNodes[NUM_NODES];

// large boolean array of connections between pairs of nodes
int connection[NUM_NODES][NUM_NODES] = {};

// smaller, condensed array of vectors of Nodes connected to a particular Node
vector<int> condConn[NUM_NODES] = {};

/****************************************************************************
**                                                                         **
**                             function main                               **
**                                                                         **
*****************************************************************************/

int main(int argc, char *argv[])/**/
//int main()
{
	/*if (argc != 2)
	{
		cout << "Intended command line inputs: .\\name, Rinterval"
			<< endl;
		return -1;
	}*/

    int randSeed = atof(argv[1]); // + atof(argv[3]); /**/
    //for (int networkNum = 0; networkNum < NUM_NETWORKS; networkNum++)
    //{ // for loop to make 100 networks (Paul)
    //int randSeed = networkNum + 1; // Paul


    // sets random seed
    srand(randSeed);

    // creates netork. If not all nodes are connected to all other nodes,
    // restart.
    bool connected = false;
    int networksTried = 1;
    while (connected == false)
    {
        // clears arrays
        for (int nodeNum = 0; nodeNum < NUM_NODES; nodeNum++)
        {
            vector<int> clearedConn;
            condConn[nodeNum] = clearedConn;

            for (int secondNodeNum = 0; secondNodeNum < NUM_NODES; secondNodeNum++)
            {
                connection[nodeNum][secondNodeNum] = 0;
            }
        }

        // creates nodes with random populations from 0 to 1. Replaces
        // nodes if they already exist
        createNodesRand();

        // creates internal connections between nodes in every module
        for (int modNum = 0; modNum < NUM_MODULES; modNum++)
        {
            createInternalConnections(modNum, PROB_INT_CONN);
        }

        cout << "Internal Connections Created" << endl;

        // creates external connections between nodes in different modules
        createExternalConnections(PROB_EXT_CONN);

        cout << "External Connections Created" << endl;

        // creates simpler version of connection array, condConn
        condenseConnections();

        // all nodes connected to centerNodeNum are added to connectedNodes
        int centerNodeNum = 0;
        vector<int> connectedNodes;
        int distance = 1;
        while (distance < MAX_DISTANCE)
        {
            vector<int> ultDistNodes;

            findNodesDistance(distance, centerNodeNum,
                condConn, ultDistNodes);

            if (ultDistNodes.size() == 0)
            {
                distance = MAX_DISTANCE;
            }
            else
            {
                connectedNodes.insert(connectedNodes.end(), ultDistNodes.begin(), ultDistNodes.end());

                distance++;
            }
        }

        // if connectedNodes contains all other Nodes in the network, use that
        // network. If not, try again.
        if (connectedNodes.size() == (NUM_NODES - 1))
        {
            cout << "Network " << networksTried << " is connected." << endl;
            connected = true;
        }
        else
        {
            networksTried++;
            cout << networksTried << " Networks Tried" << endl;
        }
    }


    cout << modCalc() << endl;

    cout << "done" << endl;

    fstream outFileStruct;
    string structFileName = //"inputstruct" + patch::to_string(randSeed)
        "../NetworkSimFixed/data/inputstruct" + patch::to_string(randSeed)
        + ".txt";
    outFileStruct.open(structFileName.c_str(), fstream::out);

    showConnectionsBoolean(outFileStruct);
    
    //} // end of for loop (Paul)
    
	getchar();
	getchar();

	return 0;
} // END function main

/****************************************************************************
**                                                                         **
**                      class Node member functions                        **
**                                                                         **
*****************************************************************************/

Node::Node()
{
	startPop = -1;
	population = 0;
	numIntConns = 0;
	numExtConns = 0;
}

Node::Node(double pop)
{
	startPop = pop;
	population = pop;
	numIntConns = 0;
	numExtConns = 0;
}

Node::Node(double pop, double firstLyap)
{
	startPop = pop;
	population = pop;
	numIntConns = 0;
	numExtConns = 0;
}

Node::Node(Node& oldNode)
{
	startPop = oldNode.get_startPop();
	population = oldNode.get_population();
	numIntConns = oldNode.get_numIntConns();
	numExtConns = oldNode.get_numExtConns();
}

/****************************************************************************
**                                                                         **
**                      functions currently in use                         **
**                                                                         **
*****************************************************************************/

// creates all the nodes in the simulation, numbered 0 to (numNodes - 1),
// grouped by module in order of module number, and stored in arrNodes. Assigns
// each node a random population 0.XX
void createNodesRand()
{
	for (int i = 0; i < NUM_NODES; i++)
	{
		double population = 0;

		Node Node1(population);
		arrNodes[i] = Node1;
	}
} // END function createNodes()

// sets up internal connections in a module. Connects any two different
// nodes in the module with probability PROB_INT_CONN
void createInternalConnections(int modNum, double probIntConn)
{
	int firstModNode = modNum * NODES_PER_MODULE;

	for (int activeNode = firstModNode; 
		activeNode < (firstModNode + NODES_PER_MODULE ); activeNode++)
	{
		for (int dependentNode = activeNode + 1; 
			dependentNode < (firstModNode + NODES_PER_MODULE); dependentNode++)
		{
			double randNum = rand();
			randNum = randNum / RAND_MAX;

			if (randNum < probIntConn)
			{
				connection[activeNode][dependentNode] = 1;
				connection[dependentNode][activeNode] = 1;
				arrNodes[activeNode].set_numIntConns(
					arrNodes[activeNode].get_numIntConns() + 1);
				arrNodes[dependentNode].set_numIntConns(
					arrNodes[dependentNode].get_numIntConns() + 1);
			}
		}
	}
} //END function createInternalConnections

// sets up all external connections in the network. Connects any two nodes in
// different modules with probability PROB_EXT_CONN
void createExternalConnections(double probExtConn)
{
	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		for (int dependentNode = activeNode + 1;
			dependentNode < NUM_NODES; dependentNode++)
		{
			int activeModule = activeNode / NODES_PER_MODULE;
			int dependentModule = dependentNode / NODES_PER_MODULE;

			if (activeModule != dependentModule)
			{
				double randNum = rand();
				randNum = randNum / RAND_MAX;

				if (randNum < probExtConn)
				{
					connection[activeNode][dependentNode] = 1;
					connection[dependentNode][activeNode] = 1;
					arrNodes[activeNode].set_numExtConns(
						arrNodes[activeNode].get_numExtConns() + 1);
					arrNodes[dependentNode].set_numExtConns(
						arrNodes[dependentNode].get_numExtConns() + 1);
				}
			}
		}
	}
} // END function createExternalConnections

// creates condConn[]<> out of connections[][]
void condenseConnections()
{
	// loops over all the nodes in the array
	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		int numIntConns = arrNodes[activeNode].get_numIntConns();
		int numExtConns = arrNodes[activeNode].get_numExtConns();
		int numConns = numIntConns + numExtConns;
		int connNum = 0;

		// loops through all the positions in activeNode's row of connections[][],
		// and wherever there's a connection, adds the connected Node's number to
		// activeNode's condConn vector
		for (int activeConn = 0; activeConn < numConns; activeConn++)
		{
			bool connected = false;

			// loops until it reaches a connection.
			while (connected == false)
			{
				if (connection[activeNode][connNum] == 1)
				{
					connected = true;

					// adds the connected node to activeNode's vector
					condConn[activeNode].push_back(connNum);
				}

				connNum++;
			}
		}
	}
} // END function condenseConnections

// displays the matrix of connections in a table, with the number of each node
// across the top and down the side. 0 = not connected, 1 = connected
void showConnectionsBoolean(fstream& outFile)
{
	outFile << "NumNodes = " << NUM_NODES << endl;
	outFile << "   ";
	for (int i = 0; i < NUM_NODES; i++)
	{
		outFile << setw(2) << i << " ";
	}
	outFile << endl;

	for (int i = 0; i < NUM_NODES; i++)
	{
		outFile << setw(2) << i << " ";
		for (int j = 0; j < NUM_NODES; j++)
		{
			outFile << setw(2) << connection[i][j] << " ";
		}
		outFile << endl;
	}
	outFile << endl;
} // END function showConnectionsBoolean()

// lists the Nodes each Node is connected to. Will only work after
// condenseConnections() has been called
void showConnectionsVec(fstream& outFile)
{
	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		int numIntConns = arrNodes[activeNode].get_numIntConns();
		int numExtConns = arrNodes[activeNode].get_numExtConns();
		int numConns = numIntConns + numExtConns;
		
		outFile << "Node " << setw(2) << activeNode << " connections: ";

		for (int conn = 0; conn < numConns; conn++)
		{
			outFile << setw(2) << condConn[activeNode].at(conn) << ", ";
		}

		outFile << endl;
	}

	outFile << endl;
} // END function showConnectionsVec

// displays the number of connections each Node has
void listNumConnections(fstream& outFile)
{
	outFile << "Internal, External, then Total Connections per Node:" << endl;

	for (int i = 0; i < NUM_NODES; i++)
	{
		int intConns = arrNodes[i].get_numIntConns();
		int extConns = arrNodes[i].get_numExtConns();
		int totConns = intConns + extConns;
		outFile << "Node " << setw(2) << i << ": " << setw(2) << intConns << ", "
			<< setw(2) << extConns << ", " << setw(2) << totConns << endl;
	}

	outFile << endl;
} // END listNumConnections

// calculates network modularity
double modCalc()
{
	double modularity = 0;

	int numConns = 0;

	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		numConns = numConns + arrNodes[activeNode].get_numIntConns()
			+ arrNodes[activeNode].get_numExtConns();
	}

	numConns = numConns / 2;

	// loop over the modules in the network
	for (int modNum = 0; modNum < NUM_MODULES; modNum++)
	{
		double numIntEdges = 0;
		double nodeDegSum = 0;
		double modComponent = 0;

		// loop over nodes in module
		for (int nodeInMod = 0; nodeInMod < NODES_PER_MODULE; nodeInMod++)
		{
			int nodeNum = (modNum * NODES_PER_MODULE) + nodeInMod;
			int nodeIntConns = arrNodes[nodeNum].get_numIntConns();
			int nodeExtConns = arrNodes[nodeNum].get_numExtConns();

			// find total internal edges and total degree of nodes
			numIntEdges = numIntEdges + nodeIntConns;
			nodeDegSum = nodeDegSum + nodeIntConns + nodeExtConns;
		}

		numIntEdges = numIntEdges / 2;

		// calculate modularity contribution for the current module
		modComponent = (numIntEdges / numConns) -
			(nodeDegSum * nodeDegSum / (4 * numConns * numConns));

		// add modularity component to overall modularity
		modularity = modularity + modComponent;
	}

	//cout << numConns << endl;

	return modularity;
}

// finds all nodes a given distance from a given node
void findNodesDistance(int ultDistance, int centerNodeNum,
	vector<int> condConn[NUM_NODES], vector<int> &ultDistNodes)
{
	// creates vector of nodes closer to or at desired distance from center node.
	vector<int> closerNodes;

	// puts center node in closer nodes vector
	closerNodes.push_back(centerNodeNum);

	// creates vector of just reached nodes, to be cleared and updated during
	// every distance loop
	vector<int> justReachedNodes;

	// puts center node in initial justReachedNodes vector
	justReachedNodes.push_back(centerNodeNum);

	int closerNodesSize = closerNodes.size();
	int numJustReachedNodes = justReachedNodes.size();

	for (int distance = 1; distance <= ultDistance; distance++)
	{
		closerNodesSize = closerNodes.size();
		numJustReachedNodes = justReachedNodes.size();
		vector<int> newReachedNodes;

		for (int numNodeInJustReached = 0; numNodeInJustReached < numJustReachedNodes; numNodeInJustReached++)
		{
			int nodeInJustReached = justReachedNodes[numNodeInJustReached];
			int numConnectedNodes = condConn[nodeInJustReached].size();

			for (int numConnectedNode = 0; numConnectedNode < numConnectedNodes; numConnectedNode++)
			{
				int connectedNode = condConn[nodeInJustReached][numConnectedNode];

				bool nodeIsCloser = false;

				for (int numCloserNode = 0; numCloserNode < closerNodesSize; numCloserNode++)
				{
					int closerNode = closerNodes[numCloserNode];

					if (connectedNode == closerNode)
					{
						nodeIsCloser = true;
					}
				}

				if (nodeIsCloser == false)
				{
					newReachedNodes.push_back(connectedNode);
					closerNodes.push_back(connectedNode);
					closerNodesSize++;
				}
			}
		}

		justReachedNodes = newReachedNodes;
	}

	ultDistNodes = justReachedNodes;
}

/****************************************************************************
**                                                                         **
**                    functions not currently in use                       **
**                                                                         **
*****************************************************************************/
