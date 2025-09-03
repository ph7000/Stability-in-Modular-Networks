/**
* Network.cpp
*
* PHYS-301 Fall 2020/Spring 2021
* Modular Network Project
* 
* Paul Heyden
* Original Author: Benjamin Stein-Lubrano
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
#include <time.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdexcept>
#include <exception>
#include <errno.h>
#include <bits/stdc++.h>

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

const int
NUM_NODES = 64,
NUM_LYAPS = 1; //NUM_NODES; //(NUM_NODES / 2);

// only used in calculating coarse-grained participation ratios. Hopefully
// we'll be able to get rid of it if/when we add the code that figures out
// what nodes belong to what modules.
const int NUM_MODULES = 8,
  NODES_PER_MODULE = 8,
  DEGREE = 4;


const int
NUM_RAND_NETWORKS = 1, // 100,
NUM_D_VALUES = 100, //= 100,
NUM_TIMESTEPS = 100000,
OUTPUT_DELAY = 50000;

// WARNING: switching the time series on causes massive amounts of data to be 
// outputted, so it should only be true if the code is run for a single 
// connection strength value. 
const bool DISPLAY_TIME_SERIES = true;

/*const double
  MU_FACTOR = 3.55; /**/

/*const double
NUM_TOTAL_CONN = 4.0,
NUM_INT_CONN = 3.2,
NUM_EXT_CONN = NUM_TOTAL_CONN - NUM_INT_CONN,
PROB_INT_CONN = (NUM_INT_CONN / 7.0), // (3.2 / 7.0),                                                                                                            
PROB_EXT_CONN = (NUM_EXT_CONN / 56.0), // (0.8 / 56.0),                                                                                                          
PERCENT_INT_CONN = (NUM_INT_CONN / NUM_TOTAL_CONN); /**/

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
	double get_oldPopulation() { return oldPopulation; }
	void set_oldPopulation(double pop) { oldPopulation = pop; }
	double get_population() { return population; }
	void set_population(double pop) { population = pop; }
	double get_oldLyapVecContr(int lyapNum) { return oldLyapVecContr[lyapNum]; }
	void set_oldLyapVecContr(int lyapNum, double pop) { oldLyapVecContr[lyapNum] = pop; }
	double get_lyapVecContr(int lyapNum) { return lyapVecContr[lyapNum]; }
	void set_lyapVecContr(int lyapNum, double pop) { lyapVecContr[lyapNum] = pop; }

	int get_numConns() { return numConns; }
	void set_numConns(int Conn) { numConns = Conn; }

private:
	double startPop;
	double oldPopulation;
	double population;
	double oldLyapVecContr[NUM_LYAPS];
	double lyapVecContr[NUM_LYAPS];
	int numConns;
}; // END Node class declaration

/****************************************************************************
**                                                                         **
**                         function declarations                           **
**                                                                         **
*****************************************************************************/

// functions
int loadNetworkStructure(string);
void condenseConnections();

void showConnectionsBoolean(fstream&);
void showConnectionsVec(fstream&);
void listNumConnections(fstream&);
void showPopulations(fstream&);
void showShortTermExp(int, fstream&);

void nodeGrowth(double);
void nodeLyapGrowth(double);
void nodeInteraction(double);
void nodeLyapInteraction(double);
void modifiedGramSchmidt(int, double*, double*);
void normalizeLyap();
void modulePrevalence(double*);
void displayLyapContr();

double modCalc();

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

// array of sum of short-time lyapunov exponents. Used to find long-time
// lyapunov exponents
double lyapExponents[NUM_LYAPS] = {};

// arrays of sums of short-time inverse participation ratios. Used to find long-
// time inverse participation ratios
double partRatioArr[NUM_LYAPS] = {};

// array of sum of short-time inverse coarse participation ratios. Used to
// find long-time inverse coarse participation ratios
double partRatioArrCoarse[NUM_LYAPS] = {};

//array of percentages of how much of the lyapunov contribuions are in the  module with the biggest lyapunov contribution
double largestContrModArr[NUM_LYAPS] = {};

//array of mu values
double muValues[10] = {3.1, 3.3, 3.5, 3.55, 3.555, 3.559, 3.56, 3.561, 3.565, 3.57};
//double muValues[3] = {3.6, 3.8, 4.0};

//array of percent internal
double PValues[1] = {0.4};  
//double PValues[4] = {0.2, 0.4, 0.6, 0.8};

/****************************************************************************
**                                                                         **
**                             function main                               **
**                                                                         **
*****************************************************************************/

int main(int argc, char *argv[]) /**/
{
    double startTime = clock();

    if (argc != 2)
	{
		cout << "Intended command line inputs: .\\name, Dinterval"
			<< endl;
		return -1;
	}

    //The following code allows for job arrays of 0-999, with the hundreds...
    //decimal place giving the mu value and the other two places giving ...
    //(DFactor - 1). For example a command argument of 142 gives mu3.3, DFactor = 43.
    //Can change to doing just one mu value if the coded-out DFactor line is used.

    string commandArgString = string(argv[1]);
    int commandArg = (int)atof(argv[1]);  
    int muPosition = 0;
    int onesDigit = commandArg % 10;
    int tensDigit = ((commandArg - onesDigit)/10) % 10;
    int DValue = onesDigit + 10*tensDigit;

    int argSize = commandArgString.size();
      
    if (argSize == 3)
    {
      muPosition = (commandArg - DValue)/100;
    }

    double DFactor = 0.01 * ((double)DValue + 1.0);
    double muFactor = muValues[muPosition];

    cout << "The size of the command argument is " << argSize << endl;
    cout << muFactor << " <- mu factor, DFactor -> " << DFactor << endl; /**/

    /*double DFactor = 0.01 * atof(argv[1]); /**/
    
    for (int PPosition = 0; PPosition < sizeof(PValues) / sizeof(PValues[0]); PPosition++)
    { // runs the entire code for each percent internal (P) value
        double percentIntConn = PValues[PPosition]; /**/
	
	string trialFilePathMu = "./../Consolidate/data/mu" + patch::to_string(muFactor);     
      
	if (mkdir(trialFilePathMu.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)        
	  {                                                                                       
	    if( errno == EEXIST ) {                                                             
	      // already exists                                                               
	    } else                                                                              
	      {                                                                                   
		// something else                                                               
		cout << "error in mu file";                                                     
	      }                                                                                   
	  }                                                                                       
                                                         
	string trialFilePathFull = trialFilePathMu + "/" + "mu" + patch::to_string(muFactor) +  "P" + patch::to_string(percentIntConn);  
	
	if (mkdir(trialFilePathFull.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)      
	  {                                                                                       
	    if( errno == EEXIST ) {                                                             
	      // already exists                                                                   
	    } else {                                                                            
	      // something else                                                               
	      cout << "error in mu file";                                                     
	    }                                                                                   
	  }                                                                                       
                                                                                                    
	if (NUM_NODES != 64)                                                                    
	  {                                                                                       
	    trialFilePathFull += "/Mo" + patch::to_string(NUM_MODULES) + "/" + "mu" + patch::to_string(muFactor) + "P" + patch::to_string(percentIntConn) + "Mo" + patch::to_string(NUM_MODULES) + "No" + patch::to_string(NODES_PER_MODULE);                                   
                                                                                                    
	    if (mkdir(trialFilePathFull.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)  
	      {                                                                                   
		if( errno == EEXIST ) {                                                         
		  // already exists                                                           
		} else {                                                                        
		  // something else                                                           
		  cout << "error in mu file";                                                 
                }                                                                                   
                                                                                                    
	      }                                                                                   
                                                                                                    
	  } 

	if (DEGREE != 4)
          {
            trialFilePathFull += "/mu" + patch::to_string(muFactor) + "P" + patch::to_string(percentIntConn) + "k" + patch::to_string(DEGREE);

            if (mkdir(trialFilePathFull.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
              {
                if( errno == EEXIST ) {
                  // already exists                                                                                                                                                                     
                } else {
                  // something else                                                                                                                                                                     
                  cout << "error in mu file";
                }

              }

          }

	string trialFileName = trialFilePathFull
	  + "/" + "D" + patch::to_string(DFactor);

	fstream outFilePart;
	string partFileName = trialFileName + "part.txt";
	outFilePart.open(partFileName.c_str(), fstream::out); /**/

	/*fstream outFilePartSquared;                                                                                                                                                   
            string partSquaredFileName = trialFileName + "partSq.txt";                                                                                                                      
            outFilePartSquared.open(partSquaredFileName.c_str(), fstream::out); /**/

	fstream outFilePartCoarse;
	string partCoarseFileName = trialFileName + "coarse.txt";
	outFilePartCoarse.open(partCoarseFileName.c_str(), fstream::out); /**/

	/*fstream outFilePartCoarseSquared;                                                                                                                                             
            string partCoarseSquaredFileName = trialFileName + "coarseSq.txt";                                                                                                              
            outFilePartCoarseSquared.open(partCoarseSquaredFileName.c_str(), fstream::out); /**/

	fstream outFileExp;
	string expFileName = trialFileName + "exp.txt";
	outFileExp.open(expFileName.c_str(), fstream::out); /**/

	/*fstream outFileExpSquared;                                                                                                                                                    
            string expSquaredFileName = trialFileName + "expSq.txt";                                                                                                                        
            outFileExpSquared.open(expSquaredFileName.c_str(), fstream::out); /**/

	fstream outFileDim;
	string dimFileName = trialFileName + "dim.txt";
	outFileDim.open(dimFileName.c_str(), fstream::out); /**/

	/*fstream outFileDimSquared;                                                                                                                                                    
            string dimSquaredFileName = trialFileName + "dimSq.txt";                                                                                                                        
            outFileDimSquared.open(dimSquaredFileName.c_str(), fstream::out); /**/

	fstream outFileMod;
	string modFileName = trialFileName + "mod.txt";
	outFileMod.open(modFileName.c_str(), fstream::out); /**/

	/*fstream outFileModSquared;                                                                                                                                                    
            string modSquaredFileName = trialFileName + "modSq.txt";                                                                                                                        
            outFileModSquared.open(modSquaredFileName.c_str(), fstream::out); /**/

	fstream outFileTimeSeries;
	fstream outFileTimeSeriesPopulation;
	//The following is only necessary for time series files
	if(DISPLAY_TIME_SERIES)
	  {
	    string timeSeriesFilePath = "./../Consolidate/data/timeSeries";

	    if (mkdir(timeSeriesFilePath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
	      {            
		if( errno == EEXIST ) {   
		  // already exists                     
		} else                         
		  {    
		    // something else                                      
		    cout << "error in mu file";                 
		  }         
	      } 

	    string timeSeriesFileName = timeSeriesFilePath
	      + "/" + "mu" + patch::to_string(muFactor)
	      + "P" + patch::to_string(percentIntConn)
	      + "D" + patch::to_string(DFactor)
	      + ".txt";

	    outFileTimeSeries.open(timeSeriesFileName.c_str(), fstream::out); /**/

	    string timeSeriesPopulationFilePath = "./../Consolidate/data/timeSeriesPopulation";

	    if (mkdir(timeSeriesPopulationFilePath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
	      {
		if( errno == EEXIST ) {
		  // already exists                                                                   
		} else
		  {
		    // something else                                                                 
		    cout << "error in mu file";
		  }
	      }

	    string timeSeriesPopulationFileName = timeSeriesPopulationFilePath
	      + "/" + "mu" + patch::to_string(muFactor)
	      + "P" + patch::to_string(percentIntConn)
	      + "D" + patch::to_string(DFactor)
	      + ".txt";

	    outFileTimeSeriesPopulation.open(timeSeriesPopulationFileName.c_str(), fstream::out); /**/
	  }

        for (int randNetwork = 1; randNetwork < NUM_RAND_NETWORKS + 1; randNetwork++)
        { // basically runs the entire code for each network structure
            int randSeed = (2 * randNetwork + 1);

            // sets random seed
            srand(randSeed);

            // resets lyapunov exponent and participation ratio arrays to 0
            for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
            {
                lyapExponents[lyapNum] = 0;
                partRatioArr[lyapNum] = 0;
                partRatioArrCoarse[lyapNum] = 0;
                largestContrModArr[lyapNum] = 0;
            }

            for (int nodeNum = 0; nodeNum < NUM_NODES; nodeNum++)
            {
                condConn[nodeNum].clear();
            }

            string inputStructFileName = "./data/P" + patch::to_string(percentIntConn);
	    
	    if (NUM_NODES != 64 || DEGREE != 4)
	      {
		inputStructFileName += "/Mo" + patch::to_string(NUM_MODULES) 
		  + "/" + "No" + patch::to_string(NODES_PER_MODULE) 
		  + "/k" + patch::to_string(DEGREE) 
		  + "P" + patch::to_string(percentIntConn) 
		  + "Mo" + patch::to_string(NUM_MODULES) 
		  + "No" + patch::to_string(NODES_PER_MODULE) 
		  + "k" + patch::to_string(DEGREE) + "/";
	      }
	    
	    inputStructFileName += "inputstruct" + patch::to_string(randNetwork) + ".txt";

            // loads Network structure from inputStructFileName.
            int numNodes = loadNetworkStructure(inputStructFileName);

            if (numNodes != NUM_NODES)
            {
                cout << "Incorrect Number of Nodes in Network" << endl;

                getchar();
                getchar();

                return 1;
            }

	    /*string trialFileName = trialFilePathFull
	      + "/" + "D" + patch::to_string(DFactor);/**/

            // creates simpler version of connection array, condConn
            condenseConnections();
        
            /*string trialFilePathMu = "./../Consolidate/data/mu" + patch::to_string(muFactor);

            if (mkdir(trialFilePathMu.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
            {
                if( errno == EEXIST ) {
                    // already exists
                } else
                {
                    // something else
                    cout << "error in mu file";
                }
            }
        
            string trialFilePathFull = trialFilePathMu + "/" + "mu" + patch::to_string(muFactor) + "P" + patch::to_string(percentIntConn);

            if (mkdir(trialFilePathFull.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
            {
                if( errno == EEXIST ) {
                // already exists
                } else {
                    // something else
                    cout << "error in mu file";
                }
            }
        
            if (NUM_NODES != 64)
            {
                trialFilePathFull += "/Mo" + patch::to_string(NUM_MODULES) + "/" + "mu" + patch::to_string(muFactor) + "P" + patch::to_string(percentIntConn) + "Mo" + patch::to_string(NUM_MODULES) + "No" + patch::to_string(NODES_PER_MODULE);
            
                if (mkdir(trialFilePathFull.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
                {
                    if( errno == EEXIST ) {
                        // already exists
                    } else {
                        // something else
                        cout << "error in mu file";
                }
                
                }

            }

            string trialFileName = trialFilePathFull
                + "/" + "D" + patch::to_string(DFactor)
                + "R" + patch::to_string(randSeed)
	        + "S" + patch::to_string(randNetwork); /**/

            /* fstream outFileStruct;
            string outStructFileName = "outputstruct"
                + patch::to_string(randNetwork) + ".txt";
            outFileStruct.open(outStructFileName.c_str(), fstream::out);

            showConnectionsBoolean(outFileStruct);
            showConnectionsVec(outFileStruct); /**/

            // actually running simlation starts here

            normalizeLyap();

            /*fstream outFilePart;
            string partFileName = trialFileName + "part.txt";
            outFilePart.open(partFileName.c_str(), fstream::out); /**/

            /*fstream outFilePartSquared;
            string partSquaredFileName = trialFileName + "partSq.txt";
            outFilePartSquared.open(partSquaredFileName.c_str(), fstream::out); /**/

            /*fstream outFilePartCoarse;
            string partCoarseFileName = trialFileName + "coarse.txt";
            outFilePartCoarse.open(partCoarseFileName.c_str(), fstream::out); /**/

            /*fstream outFilePartCoarseSquared;
            string partCoarseSquaredFileName = trialFileName + "coarseSq.txt";
            outFilePartCoarseSquared.open(partCoarseSquaredFileName.c_str(), fstream::out); /**/

            /*fstream outFileExp;
            string expFileName = trialFileName + "exp.txt";
            outFileExp.open(expFileName.c_str(), fstream::out); /**/

            /*fstream outFileExpSquared;
            string expSquaredFileName = trialFileName + "expSq.txt";
            outFileExpSquared.open(expSquaredFileName.c_str(), fstream::out); /**/

            /*fstream outFileDim;
            string dimFileName = trialFileName + "dim.txt";
            outFileDim.open(dimFileName.c_str(), fstream::out); /**/

            /*fstream outFileDimSquared;
            string dimSquaredFileName = trialFileName + "dimSq.txt";
            outFileDimSquared.open(dimSquaredFileName.c_str(), fstream::out); /**/

            /*fstream outFileMod;
            string modFileName = trialFileName + "mod.txt";
            outFileMod.open(modFileName.c_str(), fstream::out); /**/

            /*fstream outFileModSquared;
            string modSquaredFileName = trialFileName + "modSq.txt";
            outFileModSquared.open(modSquaredFileName.c_str(), fstream::out); /**/

            // sets up output files for population and lyapunov vector components
            /*fstream outFilePop;
            string popFileName = trialFileName + "pop.pgm";
            outFilePop.open(popFileName.c_str(), fstream::out);
            outFilePop << "P2" << endl
            << NUM_NODES << " " << (NUM_TIMESTEPS - OUTPUT_DELAY) << endl
            << "250" << endl; /**/

            /*fstream outFileLyap[NUM_LYAPS];
            string lyapFileName = trialFileName + "lyap";
            for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
            {
                string fileName = lyapFileName + patch::to_string(lyapNum)
                    + ".pgm";
                outFileLyap[lyapNum].open(fileName.c_str(), fstream::out);
            } /**/


            /*for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
            {
            outFileLyap[lyapNum] << "P2" << endl
            << NUM_NODES << " " << (NUM_TIMESTEPS - OUTPUT_DELAY) << endl
            << "250" << endl;
            } /**/

            // array of participation ratios, and pointer to first element
            // of array
            double partRatio[NUM_LYAPS] = {};
            double *partRatioPt = partRatio;

            // array of coarse-grained participation ratios (module
            // participation, rather than node participation), and pointer
            // to first element of array
            double partRatioCoarse[NUM_LYAPS] = {};
            double *partRatioCoarsePt = partRatioCoarse;

            //same as above but with largest mod
            double largestContrMod[NUM_LYAPS] = {};
            double *largestContrModPt = largestContrMod;

            /*double largestMod[NUM_LYAPS] = {};
            double *largestModPt = largestMod;/**/

            for (int timeStep = 0; timeStep < NUM_TIMESTEPS; timeStep++)
            {
                int checkStep = timeStep % (NUM_TIMESTEPS / 10);

                if (checkStep == 0)
                {
                    cout << "interactions " << (timeStep / (NUM_TIMESTEPS / 10))
                        << "/10 complete" << endl;
                }

                nodeGrowth(muFactor);

                nodeLyapGrowth(muFactor);

                nodeInteraction(DFactor);

                nodeLyapInteraction(DFactor);

		if (DISPLAY_TIME_SERIES)
		  {

		    if (timeStep >= OUTPUT_DELAY)
		      {
			outFileTimeSeries << timeStep << " ";
		      }

		    showShortTermExp(timeStep, outFileTimeSeries);
		  }

                modifiedGramSchmidt(timeStep, partRatioPt, partRatioCoarsePt);

                normalizeLyap();

                modulePrevalence(largestContrModPt);

                if (timeStep >= OUTPUT_DELAY)
                {
                    /*for (int j = 0; j < NUM_NODES; j++)
                    {
                        double color1 = arrNodes[j].get_population() * 250.0;
                        int color2 = color1;
                        outFilePop << setw(4) << color2;
                    }
                    outFilePop << endl; /**/

                    /*for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
                    {
                        for (int j = 0; j < NUM_NODES; j++)
                        {
                            double color1 = arrNodes[j].get_lyapVecContr(lyapNum) * 125.0 + 125.0;
                            int color2 = color1;
                            outFileLyap[lyapNum] << setw(4) << color2;
                        }
                        outFileLyap[lyapNum] << endl;
                    } /**/

                    /* outFilePart << timeStep << ": ";
                    for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
                    {
                        outFilePart << partRatio[lyapNum] << " ";
                    }
                    outFilePart << endl; /**/

                    for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
                    {
                        double oneTimePartRatio = partRatio[lyapNum];
                        partRatioArr[lyapNum] += oneTimePartRatio;

                        partRatioArrCoarse[lyapNum] += partRatioCoarse[lyapNum];
                
                        largestContrModArr[lyapNum] += largestContrMod[lyapNum];
                    } /**/
		    
		    //timeSeriesPopulation output
		    if(DISPLAY_TIME_SERIES)
		      {
			outFileTimeSeriesPopulation << timeStep << " ";
			showPopulations(outFileTimeSeriesPopulation);
		    
			//timeSeries output
			for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
			  {
			    outFileTimeSeries << partRatio[lyapNum] << " " << partRatio[lyapNum]*partRatio[lyapNum] << " ";
			
			    outFileTimeSeries << partRatioCoarse[lyapNum]<< " " << partRatioCoarse[lyapNum]*partRatioCoarse[lyapNum] << " ";
			
			    outFileTimeSeries << largestContrMod[lyapNum]<< " " << largestContrMod[lyapNum]*largestContrMod[lyapNum] << endl;
			  } /**/
		      }
                }
            }
	    
	    // Now that the time steps are complete, we can display the lyapContr's to see the largest one
	    // This is relevant any time we want to check the behavior of a time series
	    displayLyapContr();

            double lyapDimension = 0;
            double lyapExpSum = 0;
            for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
            {
                double lyapExponent = lyapExponents[lyapNum] /
                    (NUM_TIMESTEPS - OUTPUT_DELAY);
                double lyapExponentSquared = lyapExponent * lyapExponent;
                outFileExp << lyapExponent << " " << lyapExponentSquared << " ";


                if (lyapExpSum >= 0)
                {
                    if ((lyapExpSum + lyapExponent) >= 0)
                        lyapDimension++;

                    if ((lyapExpSum + lyapExponent) < 0)
                    {
                        lyapDimension -= (lyapExpSum / lyapExponent);
                    }
                }
                lyapExpSum += lyapExponent;
            }

            double lyapDimensionSquared = lyapDimension * lyapDimension;

            outFileDim << lyapDimension << " " << lyapDimensionSquared; /**/

            for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
            {
                double partRatioLong = partRatioArr[lyapNum] /
                    (NUM_TIMESTEPS - OUTPUT_DELAY);
                double partRatioLongSquared = partRatioLong * partRatioLong;

                double partRatioCoarseLong = partRatioArrCoarse[lyapNum] /
                    (NUM_TIMESTEPS - OUTPUT_DELAY);
                double partRatioCoarseLongSquared = partRatioCoarseLong * partRatioCoarseLong;

                double largestContrModLong = largestContrModArr[lyapNum] /
                    (NUM_TIMESTEPS - OUTPUT_DELAY);
                double largestContrModLongSquared = largestContrModLong * largestContrModLong;

                outFilePart << partRatioLong << " " << partRatioLongSquared << " ";
                outFilePartCoarse << partRatioCoarseLong << " " << partRatioCoarseLongSquared << " ";

                outFileMod << largestContrModLong << " " << largestContrModLongSquared << " ";
            } /**/

	    outFileExp << "\n";
	    outFileDim << "\n";
	    outFilePart << "\n";
	    outFilePartCoarse << "\n";
	    outFileMod << "\n";

            cout << "network " << randNetwork << " done" << endl;

            outFileExp.flush();
            //outFileExpSquared.flush();
            outFileDim.flush();
            //outFileDimSquared.flush();
            outFilePart.flush();
            //outFilePartSquared.flush();
            outFilePartCoarse.flush();
            //outFilePartCoarseSquared.flush();
            outFileMod.flush();
            //outFileModSquared.flush();
	    if(DISPLAY_TIME_SERIES)
	      {
		outFileTimeSeries.flush();
		outFileTimeSeriesPopulation.flush();
	      }

            /*outFileExp.close();
            //outFileExpSquared.close();
            outFileDim.close();
            //outFileDimSquared.close();
            outFilePart.close();
            //outFilePartSquared.close();
            outFilePartCoarse.close();
            //outFilePartCoarseSquared.close();
            outFileMod.close();
            //outFileModSquared.close(); /**/

        } // End loop over network structures

	outFileExp.close();
	//outFileExpSquared.close();                                                                                                                                                    
	outFileDim.close();
	//outFileDimSquared.close();                                                                                                                                                    
	outFilePart.close();
	//outFilePartSquared.close();                                                                                                                                                   
	outFilePartCoarse.close();
	//outFilePartCoarseSquared.close();                                                                                                                                             
	outFileMod.close();
	//outFileModSquared.close(); 
	if (DISPLAY_TIME_SERIES)
	  {
	    outFileTimeSeriesPopulation.close();
	    outFileTimeSeries.close();/**/ 
	  }

    } // End loop over percent internal values
    
    cout << "done" << endl;

    double endTime = clock();
    double runTime = (endTime - startTime) / CLOCKS_PER_SEC;

    cout << "runtTime = " << runTime << " seconds";

    // actually running simulation ends here

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
	oldPopulation = 0;
	population = 0;
	numConns = 0;

	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		oldLyapVecContr[lyapNum] = 0;
		lyapVecContr[lyapNum] = 0;
	}
}

Node::Node(double pop)
{
	startPop = pop;
	oldPopulation = pop;
	population = pop;
	numConns = 0;
	
	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		oldLyapVecContr[lyapNum] = 0;
		lyapVecContr[lyapNum] = 0;
	}
}

Node::Node(double pop, double firstLyap)
{
	startPop = pop;
	oldPopulation = pop;
	population = pop;
	numConns = 0;

	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		oldLyapVecContr[lyapNum] = firstLyap;
		lyapVecContr[lyapNum] = firstLyap;
	}
}

Node::Node(Node& oldNode)
{
	startPop = oldNode.get_startPop();
	oldPopulation = oldNode.get_oldPopulation();
	population = oldNode.get_population();
	numConns = oldNode.get_numConns();

	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
			oldLyapVecContr[activeNode] = oldNode.get_oldLyapVecContr(activeNode);

		for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
			lyapVecContr[activeNode] = oldNode.get_lyapVecContr(activeNode);
	}
}

/****************************************************************************
**                                                                         **
**                      functions currently in use                         **
**                                                                         **
*****************************************************************************/

int loadNetworkStructure(string inputFileName)
{
	fstream inFileStruct;
	string inStructFileName = inputFileName;
	inFileStruct.open(inStructFileName.c_str(), fstream::in);

	// loads network number of nodes
	int numNodes = 0;
	string emptyString = "";
	inFileStruct >> emptyString;
	inFileStruct >> emptyString;
	inFileStruct >> numNodes;
	cout << endl;

	// creates nodes according to number of nodes, with random populations
	// and lyapunov vector contributions
	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		double randNum1 = rand();
		double population = randNum1 / RAND_MAX;

		Node Node1(population);
		arrNodes[activeNode] = Node1;

		for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
		{
			double randNum2 = rand();
			double lyapVecContr = randNum2 / RAND_MAX;

			arrNodes[activeNode].set_lyapVecContr(lyapNum, lyapVecContr);
		}
	}

	// bypass first row of input file
	for (int emptyIndex = 0; emptyIndex < numNodes; emptyIndex++)
	{
		inFileStruct >> emptyString;
	}

	// imports connection array
	int activeNode = 0;
	for (activeNode = 0; activeNode < numNodes; activeNode++)
	{
		inFileStruct >> emptyString;

		int connExist = 0;
		for (int dependentNode = 0; dependentNode < numNodes; dependentNode++)
		{
			inFileStruct >> connExist;
			connection[activeNode][dependentNode] = connExist;

			if (connExist == 1)
			{
				arrNodes[activeNode].set_numConns(arrNodes[activeNode].get_numConns() + 1);
			}
		}
	}

	return numNodes;
} // END function loadNetworkStructure

// creates condConn[]<> out of connection[][]
void condenseConnections()
{
	// loops over all the nodes in the array
	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		int numConns = arrNodes[activeNode].get_numConns();
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
	outFile << "Boolean Connections Matrix" << endl;
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
		int numConns = arrNodes[activeNode].get_numConns();
		
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
		int totConns = arrNodes[i].get_numConns();
		outFile << "Node " << setw(2) << i << ": " << setw(2) << totConns << endl;
	}

	outFile << endl;
} // END listNumConnections

// displays the population of each Node
void showPopulations(fstream& outFile)
{
	for (int i = 0; i < NUM_NODES; i++)
	{
		//outFile << "Node " << setw(2) << i << ": ";
	  outFile << setprecision(12) << arrNodes[i].get_population() << " ";
	}
	outFile << endl;
} // END function showPopulations()

// lists the Nodes each Node is connected to. Will only work after
// condenseConnections() has been called                                       
void showShortTermExp(int timeStep, fstream& outFile)
{
  for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
    {
      // calculates lyapunov vector magnitudes and then normalizes lyapunov     
      // vector. Calculat                      
      double lyapMagSqr = 0;
      for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
        {
          double lyapContr = arrNodes[activeNode].get_lyapVecContr(lyapNum);
          lyapMagSqr += (lyapContr * lyapContr);
        }

      double lyapMag = sqrt(lyapMagSqr);
      if (timeStep >= OUTPUT_DELAY)
        {
          double lyapExp = (log(lyapMagSqr)) / 2;
          double lyapExpSq = lyapExp*lyapExp;
          outFile << lyapExp << " " << lyapExpSq << " ";
        }
    }
} // END function showShortTermExp 

// calculates new Node population based on the Node's inherent growth
void nodeGrowth(double muFactor)
{
	for (int numActiveNode = 0; numActiveNode < NUM_NODES; numActiveNode++)
	{
		double activeNodePop = arrNodes[numActiveNode].get_population();
		double newActiveNodePop = muFactor * activeNodePop * (1 - activeNodePop);

		arrNodes[numActiveNode].set_oldPopulation(newActiveNodePop);
		arrNodes[numActiveNode].set_population(newActiveNodePop);
	}
} // END function nodeGrowth

// calculates each node's lyapunov vector component ignoring D
void nodeLyapGrowth(double muFactor)
{
	for (int numActiveNode = 0; numActiveNode < NUM_NODES; numActiveNode++)
	{
		double activeNodePop = arrNodes[numActiveNode].get_population();

		for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
		{
			double activeNodeLyap = arrNodes[numActiveNode].get_lyapVecContr(lyapNum);
			double newActiveNodeLyap = muFactor * activeNodeLyap * (1 - (2 * activeNodePop));
			arrNodes[numActiveNode].set_oldLyapVecContr(lyapNum, newActiveNodeLyap);
			arrNodes[numActiveNode].set_lyapVecContr(lyapNum, newActiveNodeLyap);
		}
	}
} // END function nodeLyapGrowth

// calculates new Node population from contributions from connected Nodes
void nodeInteraction(double DFactor)
{
	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		int numActiveConn = arrNodes[activeNode].get_numConns();
		double oldPop = arrNodes[activeNode].get_oldPopulation();
		double newPop = (1 - DFactor) * oldPop;

		// adds all the connection contributions to the Node's population
		for (int nthConnection = 0; nthConnection < numActiveConn; nthConnection++)
		{
			int connectedNode = condConn[activeNode][nthConnection];
			newPop += (DFactor / numActiveConn)
				* arrNodes[connectedNode].get_oldPopulation();
		}

		// sets new population equal to the one calculated
		arrNodes[activeNode].set_population(newPop);
	}

	for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
	{
		double newPop = arrNodes[activeNode].get_population();
		arrNodes[activeNode].set_oldPopulation(newPop);
	}
} // END function nodeInteraction

// calculates each node's lyapunov vector component after D interaction
void nodeLyapInteraction(double DFactor)
{
	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		for (int i = 0; i < NUM_NODES; i++)
		{
			int numConn = arrNodes[i].get_numConns();
			double oldLyap = arrNodes[i].get_oldLyapVecContr(lyapNum);
			double newLyap = (1 - DFactor) * oldLyap;

			// adds all the connection contributions to the Node's lyapunov
			// component
			for (int j = 0; j < numConn; j++)
			{
				int connectedNode = condConn[i][j];
				newLyap += (DFactor / numConn)
					* arrNodes[connectedNode].get_oldLyapVecContr(lyapNum);
			}

			// sets new lyapunov component equal to the one calculated
			arrNodes[i].set_lyapVecContr(lyapNum, newLyap);
		}

		for (int i = 0; i < NUM_NODES; i++)
		{
			double newLyap = arrNodes[i].get_lyapVecContr(lyapNum);
			arrNodes[i].set_oldLyapVecContr(lyapNum, newLyap);
		}
	}
} // END function nodeInteraction

// uses Gram-Schmidt process, subtracting the projection of each lyapunov
// vector onto each more dominant vectors, to get the less-dominant lyapunov
// vectors. Modified: substracts each vector from later vectors instead of one
// at a time. Normalizes lyapunov vectors and calculates short-time exponents
// as it goes. Also calculates participation ratios.
void modifiedGramSchmidt(int timeStep, double *partRatioPt, double *partRatioCoarsePt)
{
	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		// calculates lyapunov vector magnitudes and then normalizes lyapunov
		// vector. Calculates lyapunov exponent
		double lyapMagSqr = 0;
		for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
		{
			double lyapContr = arrNodes[activeNode].get_lyapVecContr(lyapNum);
			lyapMagSqr += (lyapContr * lyapContr);
		}

		double lyapMag = sqrt(lyapMagSqr);
		if (timeStep >= OUTPUT_DELAY)
		{
			double lyapExp = (log(lyapMagSqr)) / 2;
			lyapExponents[lyapNum] += lyapExp;
		}

		for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
		{
			double oldLyapContr = arrNodes[activeNode].get_lyapVecContr(lyapNum);
			double newLyapContr = oldLyapContr / lyapMag;
			arrNodes[activeNode].set_lyapVecContr(lyapNum, newLyapContr);
		}

		// calculates participation ratio
		double *activePartRatio = partRatioPt + lyapNum;
		*activePartRatio = 0;
		for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
		{
			double lyapContr = arrNodes[activeNode].get_lyapVecContr(lyapNum);
			double lyapContrSqr = lyapContr * lyapContr;
			*activePartRatio += lyapContrSqr * lyapContrSqr;
		}
		*activePartRatio = 1.0 / (*activePartRatio);

		// this part calculates the coarse-grained participation ratios
		double *activePartRatioCoarse = partRatioCoarsePt + lyapNum;
		*activePartRatioCoarse = 0;
		for (int activeMod = 0; activeMod < NUM_MODULES; activeMod++)
		{
			// number of first node in module
			int modFirstNode = activeMod * NODES_PER_MODULE;

			double lyapContrCoarse = 0;

			// loop over nodes in module
			for (int activeNode = 0; activeNode < NODES_PER_MODULE; activeNode++)
			{
				double lyapContrNode = arrNodes[modFirstNode + activeNode].get_lyapVecContr(lyapNum);
				lyapContrCoarse += (lyapContrNode * lyapContrNode);
			}

			*activePartRatioCoarse += lyapContrCoarse * lyapContrCoarse;
		}
		*activePartRatioCoarse = 1.0 / (*activePartRatioCoarse); /**/

		// modified gram-Schmidts the vectors: one lyapunov vector at a time, it
		// subtracts the projections of later vectors onto that vector from the
		// later vectors
		for (int lyapSubNum = lyapNum + 1; lyapSubNum < NUM_LYAPS; lyapSubNum++)
		{
			double dotProduct = 0;
			// calculates the dot product of the active lyapunov vector with
			// a later vector
			for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
			{
				double activeLyapContr =
					arrNodes[activeNode].get_lyapVecContr(lyapNum);
				double subLyapContr =
					arrNodes[activeNode].get_lyapVecContr(lyapSubNum);
				dotProduct += (activeLyapContr * subLyapContr);
			}

			// subtracts the projection of the later lyapunov vector onto the
			// active vector from the later vector (componentwise)
			for (int activeNode = 0; activeNode < NUM_NODES; activeNode++)
			{
				double lyapContr =
					arrNodes[activeNode].get_lyapVecContr(lyapNum);
				double subLyapContr =
					arrNodes[activeNode].get_lyapVecContr(lyapSubNum);
				double projection = (lyapContr * dotProduct);
				double newLyapContr = subLyapContr - projection;
				arrNodes[activeNode].set_lyapVecContr(lyapSubNum, newLyapContr);
			}
		}
	}
} // END function modifiedGramSchmidt

// normalizes the Lyapunov vector to magnitude 1, so the next round of lyapunov
// vector calculations are 
void normalizeLyap()
{
	for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
	{
		double lyapMag = 0;
		for (int i = 0; i < NUM_NODES; i++)
		{
			double lyapContr = arrNodes[i].get_lyapVecContr(lyapNum);
			lyapMag += (lyapContr * lyapContr);
		}
		lyapMag = sqrt(lyapMag);

		for (int i = 0; i < NUM_NODES; i++)
		{
			double oldLyapContr = arrNodes[i].get_lyapVecContr(lyapNum);
			double newLyapContr = oldLyapContr / lyapMag;
			arrNodes[i].set_lyapVecContr(lyapNum, newLyapContr);
		}
	}
} // END function normalizeLyap

// this part calculates the coarse-grained participation ratios        
void modulePrevalence(double *largestContrMod)
{
        for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
        {
                double *activeLargestContrMod = largestContrMod + lyapNum;
		*activeLargestContrMod = 0;

		double maxModTotalContr = 0;
	
		double maxModTotalCheck = 0;

		for (int activeMod = 0; activeMod < NUM_MODULES; activeMod++)
		{
			// number of first node in module                              
                        int modFirstNode = activeMod * NODES_PER_MODULE;

			double modTotalContr = 0;

			// loop over nodes in module                                   
			for (int activeNode = 0; activeNode < NODES_PER_MODULE; activeNode++)
			{
                                double lyapContrNode = arrNodes[modFirstNode + activeNode].get_lyapVecContr(lyapNum);
                                modTotalContr += (lyapContrNode * lyapContrNode);
			}

			//modTotalContr = sqrt(modTotalContr);

			if (modTotalContr > maxModTotalContr)
			{
			          maxModTotalContr = modTotalContr;
			}
			maxModTotalCheck += modTotalContr;
		}
		*activeLargestContrMod = maxModTotalContr;
	}
}


// outputs the lyapContr for each node. Finds the biggest one and outputs it. 
void displayLyapContr()
{
  for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
    {
      /*for (int i = 0; i < NUM_NODES; i++)
	{
	  int numConn = arrNodes[i].get_numConns();
	  double oldLyap = arrNodes[i].get_oldLyapVecContr(lyapNum);
	  double newLyap = (1 - DFactor) * oldLyap;

	  // adds all the connection contributions to the Node's lyapunov                                        
	  // component                                                                                           
	  for (int j = 0; j < numConn; j++)
	    {
	      int connectedNode = condConn[i][j];
	      newLyap += (DFactor / numConn)
		* arrNodes[connectedNode].get_oldLyapVecContr(lyapNum);
	      /**/

      for (int i = 0; i < NUM_NODES; i++)
	{
	  double newLyap = arrNodes[i].get_lyapVecContr(lyapNum);
	  cout << i << " " << newLyap << ", ";
	}
      cout << endl;

      double lyapHolder = 0.0;
      int nodeHolder = 0;
      for (int i = 0; i < NUM_NODES; i++)
        {
          double newLyap = arrNodes[i].get_lyapVecContr(lyapNum);
	  if (abs(newLyap) > abs(lyapHolder))
	    {
	      lyapHolder = newLyap;
	      nodeHolder = i;
	    }
	}
      cout << nodeHolder << " " << lyapHolder << ", ";
    
    }
} // END function displayLyap 
/****************************************************************************
**                                                                         **
**                    functions not currently in use                       **
**                                                                         **
*****************************************************************************/
