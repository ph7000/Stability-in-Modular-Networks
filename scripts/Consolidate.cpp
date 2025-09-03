/**
* Consolidate.cpp
*
* PHYS-301 Spring 2018
* Modular Network Project
* 
* Benjamin Stein-Lubrano
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

const int NUM_NODES = 64,
NUM_LYAPS = NUM_NODES, // NUM_NODES / 2
NUM_RAND_NETWORKS = 100, // = 100,
NUM_D_VALUES = 100; // = 100;

const int NUM_FILE_TYPES = 8;
const int NUM_FILE_TYPES_AVERAGES = NUM_FILE_TYPES;
const int NUM_FILE_TYPES_CON_AVGS = NUM_FILE_TYPES;

const double MU_FACTOR = 3.56; // 3.58;

/****************************************************************************
**                                                                         **
**                             function main                               **
**                                                                         **
*****************************************************************************/

int main()
{
	/************************************************************************
	**                          consolidation start                        **
	*************************************************************************/

	for (int fileDNum = 0; fileDNum < NUM_D_VALUES; fileDNum++)
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
			}

			string outputFileName = "./data/mu" + patch::to_string(MU_FACTOR)
				+ "D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ ".txt";

			output.open(outputFileName.c_str(), fstream::out);

			// loop over random networks
            for (int randNetwork = 1; randNetwork < NUM_RAND_NETWORKS + 1; randNetwork++)
			{
				int randSeed = (2 * randNetwork) + 1;
				// one input stream for each network
				fstream input;

				string inputFileName = "./data/mu" + patch::to_string(MU_FACTOR)
					+ "D" + patch::to_string(Dfactor)
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
			}

			string inputFileName = "./data/mu" + patch::to_string(MU_FACTOR)
				+ "D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ ".txt";

			input.open(inputFileName.c_str());

			if (!input)
				cout << "not in file " << fileDNum << endl;

			string outputFileName = "./data/mu" + patch::to_string(MU_FACTOR)
				+ "D" + patch::to_string(Dfactor)
				+ fileType.c_str()
				+ "avg" + ".txt";

			output.open(outputFileName.c_str(), fstream::out);

			// for exp, part, partSq, and coarse, NUM_LYAPS elements per row.
			// for dim, only one element to be averaged per row
			double fileTypeSum[NUM_LYAPS] = {};
			double inputNumber = 0;
			for (int numNetwork = 1; numNetwork < NUM_RAND_NETWORKS + 1; numNetwork++)
			{
				if ( fileType == "dim" || fileType == "dimSq" )
				{
					input >> inputNumber;
					fileTypeSum[0] = fileTypeSum[0] + inputNumber;
				}
				else
				{
					for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
					{
						input >> inputNumber;
						fileTypeSum[lyapNum] = fileTypeSum[lyapNum] + inputNumber;
					}
				}
			}

			if (fileType == "dim" || fileType == "dimSq" )
			{
				double lyapExpAvg = fileTypeSum[0] / NUM_RAND_NETWORKS;
				output << lyapExpAvg << " ";
			}
			else
			{
				for (int lyapNum = 0; lyapNum < NUM_LYAPS; lyapNum++)
				{
					double lyapExpAvg = fileTypeSum[lyapNum] / NUM_RAND_NETWORKS;
					output << lyapExpAvg << " ";
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
		}

		string outputFileName = "./data/mu" + patch::to_string(MU_FACTOR)
			+ fileType.c_str()
			+ "avg" + ".txt";

		output.open(outputFileName.c_str(), fstream::out);

		for (int fileNum = 0; fileNum < NUM_D_VALUES; fileNum++)
		{
			fstream input;

			//double Dfactor = (fileNum / 100.0) + 1 / NUM_D_VALUES;
            double Dfactor = (1.0 + fileNum) / NUM_D_VALUES;
			string inputFileName = "./data/mu" + patch::to_string(MU_FACTOR) + "D"
				+ patch::to_string(Dfactor)
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
