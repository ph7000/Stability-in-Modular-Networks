/**
* Consolidate.cpp
*
* PHYS-302 Spring 2021
* Modular Network Project
* 
* Updated by Paul Heyden
* Original author: Benjamin Stein-Lubrano
* Professor David Egolf
*
**/

/*
* Code purpose: consolidate output data from multiple files into one file
*/

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

const int
NUM_MODULES = 8,
NODES_PER_MODULE = 16,
NUM_NODES = NUM_MODULES * NODES_PER_MODULE,
NUM_LYAPS = 1, // NUM_NODES / 2
NUM_RAND_NETWORKS = 100, // = 100,
NUM_D_VALUES = 100; // = 100;

const int NUM_FILE_TYPES = 5; //8;
const int NUM_FILE_TYPES_AVERAGES = NUM_FILE_TYPES;
const int NUM_FILE_TYPES_CON_AVGS = NUM_FILE_TYPES;

/*double MU_FACTOR = 3.56; // 3.58;

double
NUM_TOTAL_CONN = 4.0,
NUM_INT_CONN = 0.8,
NUM_EXT_CONN = NUM_TOTAL_CONN - NUM_INT_CONN,
PROB_INT_CONN = (NUM_INT_CONN / 7.0), // (3.2 / 7.0),                                                                                                            
PROB_EXT_CONN = (NUM_EXT_CONN / 56.0), // (0.8 / 56.0),                                                                                                          
PERCENT_INT_CONN = (NUM_INT_CONN / NUM_TOTAL_CONN);/**/

/****************************************************************************
**                                                                         **
**                             function main                               **
**                                                                         **
*****************************************************************************/
int main(int argc, char *argv[])
{
    //array of mu values. This is used in the for loop to automate the consolidation of network data instead of running manually one by one
    double muValues[10] = {3.1, 3.3, 3.5, 3.55, 3.555, 3.559, 3.56, 3.561, 3.565, 3.58};
    // double muValues[1] = {3.8};
    // array of P values
    double PValues[4] = {0.2,0.4,0.6,0.8};
    // double PValues[1] = {0.5};
    // double kValues[4] = {4,6,8,12};
    // double kValues[1] = {4};
    
    int cmd = int(atof(argv[1]));
    int muPosition = (cmd) % 10;
    int kOrPPosition = ((cmd - muPosition) / 10);
    // the above would select one mu and one P value from the list of 10 mu and 4 k values.
    // Since the command line argument allows for multiple processes to run in parallel, typing in 
    // 0-39 would run all 40 combinations of mu and P. For example, 
    // 31 would give P value 3 in the P array (0.8) and mu value 1 in the array (3.3).

    double MU_FACTOR = muValues[muPosition];
    
    double PERCENT_INT_CONN = 0.0;
    double DEGREE = 0.0;
    // now we can get the value of k or P based on the command line argument and the length of each vector.
    int kSize = sizeof(kValues) / sizeof(kValues[0]);
    int PSize = sizeof(PValues) / sizeof(PValues[0]);
    if (kSize == 1)
    { // then we iterate over P values
        DEGREE = kValues[0];
        PERCENT_INT_CONN = PValues[kOrPPosition];
    }
    if (PSize == 1)
    { // then we iterate over k values
        PERCENT_INT_CONN = PValues[0];
        DEGREE = kValues[kOrPPosition];
    }
    cout << "degree: " << patch::to_string(DEGREE) << " percent internal: " << patch::to_string(PERCENT_INT_CONN) << patch::to_string(MU_FACTOR) << endl;

    // Input files for the specific mu, P, M, N, and k values of this network.
    string directoryName = "./data/mu" + patch::to_string(MU_FACTOR)
        + "/P" + patch::to_string(PERCENT_INT_CONN)
        + "/M" + patch::to_string(NUM_MODULES)
        + "/N" + patch::to_string(NODES_PER_MODULE)
        + "/mu" + patch::to_string(MU_FACTOR)
        + "P" + patch::to_string(PERCENT_INT_CONN)
        + "M" + patch::to_string(NUM_MODULES)
        + "N" + patch::to_string(NODES_PER_MODULE)
        + "k" + patch::to_string(DEGREE);
    
    cout << directoryName << endl;
	/************************************************************************
	**                          consolidation start                        **
	*************************************************************************/

        /*for (int fileDNum = 0; fileDNum < NUM_D_VALUES; fileDNum++)
	{
		//double Dfactor = (fileDNum / 100.0) + 1 / NUM_D_VALUES;
	        double Dfactor = (1.0 + fileDNum) / NUM_D_VALUES;
		
		for (int numFileType = 0; numFileType < NUM_FILE_TYPES; numFileType++)
		{
			fstream output;

			string fileType;
			switch (numFileType)
			{
		            case 0: fileType = "exp";
		       	        break;
	         	    case 1: fileType = "expSq";
       		  		break;
       		  	    case 2: fileType = "part";
	       	  		break;
	       	  	    case 3: fileType = "partSq";
	       	  		break;
	       	     	    case 4: fileType = "coarse";
	       	  		break;
	                    case 5: fileType = "coarseSq";
		      		break;
			    case 6: fileType = "dim";
		       		break;
			    case 7: fileType = "dimSq";
		       		break;
			    case 8: fileType = "mod";
				break;
			    case 9: fileType = "modSq";
		       	        break;
			}
			
			string outputFileName = directoryName + "/D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ ".txt";

			output.open(outputFileName.c_str(), fstream::out);

			// loop over random networks
                        for (int randNetwork = 1; randNetwork < NUM_RAND_NETWORKS + 1; randNetwork++)
			{
				int randSeed = (2 * randNetwork) + 1;
				// one input stream for each network
				fstream input;

				string inputFileName = directoryName + "/D" + patch::to_string(Dfactor)
					+ "R" + patch::to_string(randSeed)
					+ "S" + patch::to_string(randNetwork)
					+ fileType.c_str()
					+ ".txt";

				input.open(inputFileName.c_str());

				if (!input)
				{
					cout << "not in file " << fileType << " consolidation, Network # = "
						<< randNetwork << ", Dfactor = " << Dfactor << endl;
				}

				double inputNumber = 0;
				while (input >> inputNumber)
				{
					output << inputNumber << " ";
				}

				output << endl;

				input.close();

				remove(inputFileName.c_str());
			}

			output.flush();
			output.close();
		}
	}

	cout << "consolidation done" << endl;

	/************************************************************************
	**                           consolidation end                         **
	*************************************************************************/

	/************************************************************************
	**                             averaging start                         **
	*************************************************************************/

	// exp and part averaging
	for (int fileDNum = 0; fileDNum < NUM_D_VALUES; fileDNum++)
	{
        //double Dfactor = (fileDNum / 100.0) + 1 / NUM_D_VALUES;
        double Dfactor = (1.0 + fileDNum) / NUM_D_VALUES;

		for (int numFileType = 0; numFileType < NUM_FILE_TYPES_AVERAGES; numFileType++)
		{
			// one input stream
			fstream input;

			// one output stream
			fstream output;

			string fileType;
			switch (numFileType)
			{
		      	    case 0: fileType = "exp";
		                break;
				//case 1: fileType = "expSq";
			  	//break;
			    case 1: fileType = "part";
		       		break;
				//case 3: fileType = "partSq";
			  	//break;
			    case 2: fileType = "coarse";
			  	break;
				//case 5: fileType = "coarseSq";
			  	//break;
			    case 3: fileType = "dim";
			  	break;
				//case 7: fileType = "dimSq";
			      //break;
			    case 4: fileType = "mod";
			        break;
				//case 9: fileType = "modSq";
			        //break;
			}

			string inputFileName = directoryName + "/D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ ".txt";

			input.open(inputFileName.c_str());

			if (!input)
				cout << "not in file " << fileDNum << endl;

			string outputFileName = directoryName + "/D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ "avg" + ".txt";

			output.open(outputFileName.c_str(), fstream::out);

			// for exp, part, partSq, and coarse, (NUM_LYAPS * 2) elements per row.
			// for dim, only two elements to be averaged per row
			double fileTypeSum[NUM_LYAPS] = {};
			double fileTypeSqSum[NUM_LYAPS] = {};
			double inputNumber = 0;
			double inputNumberSq = 0;
			for (int numNetwork = 1; numNetwork < NUM_RAND_NETWORKS + 1; numNetwork++)
			{
			  if ( fileType == "dim" )//|| fileType == "dimSq" )
				{
					input >> inputNumber;
					fileTypeSum[0] = fileTypeSum[0] + inputNumber;
					input >> inputNumberSq;
					fileTypeSqSum[0] = fileTypeSqSum[0] + inputNumberSq;
				}
				else
				{
					for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
					{
						input >> inputNumber;
						fileTypeSum[lyapNum] = fileTypeSum[lyapNum] + inputNumber;
						input >> inputNumberSq;
						fileTypeSqSum[lyapNum] = fileTypeSqSum[lyapNum] + inputNumberSq;
					}
				}
			}

			if (fileType == "dim" )//|| fileType == "dimSq" )
			{
				double valueAvg = fileTypeSum[0] / NUM_RAND_NETWORKS;
				output << valueAvg << " ";
				double valueSqAvg = fileTypeSqSum[0] / NUM_RAND_NETWORKS;
				output << valueSqAvg << " ";
			}
			else
			{
				for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
				{
					double valueAvg = fileTypeSum[lyapNum] / NUM_RAND_NETWORKS;
					output << valueAvg << " ";
					double valueSqAvg = fileTypeSqSum[lyapNum] / NUM_RAND_NETWORKS;
					output << valueSqAvg << " ";
				}
			}

			output.flush();

			remove(inputFileName.c_str());
		}
	}

	cout << "averaging done" << endl;

	/************************************************************************
	**                              averaging end                          **
	*************************************************************************/

	/************************************************************************
	**                       consolidate averages start                    **
	*************************************************************************/

	// exp part starts
	for (int numFileType = 0; numFileType < NUM_FILE_TYPES_CON_AVGS; numFileType++)
	{

		//fstream input[NUM_D_VALUES];

		fstream output;

		string fileType;
		switch (numFileType)
		{
		    case 0: fileType = "exp";
	       		break;
			//case 1: fileType = "expSq";
		  	//break;
	            case 1: fileType = "part";
		  	break;
			//case 3: fileType = "partSq";
		  	//break;
		    case 2: fileType = "coarse";
		 	break;
			//case 5: fileType = "coarseSq";
		  	//break;
		    case 3: fileType = "dim";
		        break;
			//case 7: fileType = "dimSq";
		  	//break;
		    case 4: fileType = "mod";
		        break;
			//case 9: fileType = "modSq";
		        //break;
		}

		string outputFileName = directoryName + "/" + fileType.c_str()
			+ "avg" + ".txt";

		output.open(outputFileName.c_str(), fstream::out);

		for (int fileNum = 0; fileNum < NUM_D_VALUES; fileNum++)
		{
			fstream input;

			//double Dfactor = (fileNum / 100.0) + 1 / NUM_D_VALUES;
                        double Dfactor = (1.0 + fileNum) / NUM_D_VALUES;
			string inputFileName = directoryName + "/D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ "avg" + ".txt";

			input.open(inputFileName.c_str());

			if (!input)
				cout << "not in file " << fileNum << endl;

			double inputNumber = 0;

			while (input >> inputNumber)
			{
				output << inputNumber << " ";
			}

			output << endl;

			input.close();

			remove(inputFileName.c_str());
		}

		output.flush();
		output.close();
	}

	cout << "consolidate averages done" << endl;

	/************************************************************************
	**                       consolidate averages start                    **
	*************************************************************************/
    
	cout << "done done" << endl;

	return 0;
}
