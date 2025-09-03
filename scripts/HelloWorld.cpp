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
#include <sys/stat.h>
#include <sys/types.h>
#include <stdexcept>
#include <exception>

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

const int
NUM_NODES = 64,
NUM_LYAPS = NUM_NODES; //(NUM_NODES / 2);

// only used in calculating coarse-grained participation ratios. Hopefully
// we'll be able to get rid of it if/when we add the code that figures out
// what nodes belong to what modules.
const int NUM_MODULES = 8,
NODES_PER_MODULE = 8;

const int
NUM_RAND_NETWORKS = 100, // 100,
NUM_D_VALUES = 100, //= 100,
NUM_TIMESTEPS = 100000,
OUTPUT_DELAY = 50000;

double
MU_FACTOR = 3.56; // 3.58;

/*class Node
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
}; // END Node class declaration*/


//Node arrNodes[64];

//void listNumConnections(fstream&);

int main(int argc, char *argv[])
{
    double mu_int_array[12][2] = {{3.56, 0.8},{3.56, 1.6},{3.56, 2.4},{3.56, 3.2},{3.561, 1.6},{3.565, 1.6},{3.57, 1.6},{3.6, 1.6},{3.559, 1.6},{3.555, 1.6},{3.55, 1.6},{4.0, 1.6}};
    /*for(int i = 0; i<12; i++)
    {
        cout << mu_int_array[i][0] << " " << mu_int_array[i][1] << endl;
    }*/
    
    int i = int(atof(argv[1]));
    
    double
    MU_FACTOR = mu_int_array[i][0],
    INT_FACTOR = mu_int_array[i][1];
    
    cout << patch::to_string(MU_FACTOR) << " " << INT_FACTOR << endl;
    
    /*int check;
    string mu = patch::to_string(MU_FACTOR);
    
    int n = mu.length();
  
    // declaring character array
    char dirname[n + 1];
  
    // copying the contents of the
    // string to char array
    strcpy(dirname, mu.c_str());
    //clrscr();*/
  
    /*if (mkdir("geeksforgeeks", 0777) == -1)
        cerr << "Error :  " << strerror(errno) << endl;
  
    else
        cout << "Directory created";*/
    
    //int status;
    //status = mkdir("./HI/hi", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    string finalpath = patch::to_string(MU_FACTOR);
    
    string helloworld = "hello";
    helloworld += "world";
    
    cout << helloworld << endl;
    
    if (mkdir(finalpath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if( errno == EEXIST ) {
            // already exists
        } else {
       // something else
        //std::cout << "cannot create sessionnamefolder error:" << strerror(errno) << std::endl;
        //throw std::runtime_exception( strerror(errno) );
            cout << "error in mu file";
        }
    }
    cout << "still going" << endl;
    
    return 0;
}


/*int main ()
{
    /*string name = "John";
    string family = "Smith";
    name += " K. ";         // c-string
    name += family;         // string
    name += '\n';           // character

    cout << name;*/
    
    /*string trialFileName = "mu";
    trialFileName += "DRS";
    
    cout << trialFileName;
    
    /*fstream outFileExp;
    string expFileName = trialFileName + "exp.txt";
    outFileExp.open(expFileName.c_str(), fstream::out); /**/
    
    //outFileExp << "wassup ";
/*}

/*std::vector<int> myvector (3,100);
std::vector<int>::iterator it;

it = myvector.begin();
it = myvector.insert ( it , 200 );

myvector.insert (it,2,300);

// "it" no longer valid, get a new one:
it = myvector.begin();

std::vector<int> anothervector (2,400);
myvector.insert (it+2,anothervector.begin(),anothervector.end());

int myarray [] = { 501,502,503 };
myvector.insert (myvector.begin()+1, myarray, myarray+3);

std::cout << "myvector contains:";
for (it=myvector.begin(); it<myvector.end(); it++)
  std::cout << ' ' << *it;
std::cout << '\n';

return 0;/*

/*int main()
{
    vector<int> &a;
    vector<int> &b;
    for (int count; count < 10; count++)
    {
        a.push_back(count);
        b.push_back(count*2+1);
    }
    //for (int i; i < 10; i++)
    //{
    //    cout << a[i] << " " << b[i] << endl;
    //}
    
    
    //a.insert(a.end(), b.begin(), b.end());
    for (int i; i < 10; i++)
    {
        int c = a[i];
        cout << a[i] << endl;
    }
    
    
    
    return -1;
}
    /*const string ex = "exam.txt";
    fstream outF;
    outF.open(ex);
    listNumConnections(outF);*/
    /*const string test = "test.txt";
    std::ofstream outfile (test);

    outfile << "my text here!" << endl;

    outfile.close();
    
    
    fstream PaulsFile;
    string PaulsFileName = "PaulsFile.txt";
    PaulsFile.open(PaulsFileName.c_str(), fstream::out);
    PaulsFile << "Paul's Text" << endl;
    for (int i=0; i<100; i++)
    {
        PaulsFile << "new line" << endl;
    }
    PaulsFile << endl;*/
    
    /*fstream inFileStruct;
    string inStructFileName = "inputstruct1.txt";
    inFileStruct.open(inStructFileName.c_str(), fstream::in);

    // loads network number of nodes
    int numNodes = 0;
    string emptyString = "";
    inFileStruct >> emptyString;
    inFileStruct >> emptyString;
    inFileStruct >> numNodes;
    cout << numNodes;
    cout << endl;*/
    /*arrNodes = {1, 2, 3}
    cout << arrNodes << endl;
    delete[] arrNodes;
    cout << arrNodes << endl;
    
    
    return 0;
}

/*int main(int argc, char *argv[])
//int main()
{
    cout << "Hello World"
            << endl;
    
    fstream inFileStruct;

    // loads network number of nodes
    int numNodes = 0;
    string emptyString = "";
    inFileStruct >> emptyString;
    inFileStruct >> emptyString;
    inFileStruct >> numNodes;
    cout << emptyString << emptyString << numNodes << "h" << endl;
    
    
    char a = getchar();
    cout << a;

}*/

/*// writing on a text file
#include <iostream>
#include <fstream>

int main () {
  /*ofstream myfile ("example.txt");
  if (myfile.is_open())
  {
    myfile << "This is a line.\n";
    myfile << "This is another line.\n";
    myfile.close();
  }
  else cout << "Unable to open file";
    
  int randNetwork = 1;
  string inputStructFileName = "inputstruct" + patch::to_string(randNetwork)
  + ".txt";

  fstream inFileStruct;
  string inStructFileName = inputStructFileName;
  inFileStruct.open(inStructFileName.c_str(), fstream::in);

  // loads network number of nodes
  int numNodes = 0;
  string emptyString = "";
  inFileStruct >> emptyString;
  cout << emptyString;
  inFileStruct >> emptyString;
  inFileStruct >> numNodes;
  cout << emptyString << numNodes << endl;
    
    
  string fileName = "example.txt";
  string line;
  ifstream myfile (fileName);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      cout << line << '\n';
    }
    myfile.close();
  }
  else cout << "Unable to open file";

  return 0;
}*/

/*void listNumConnections(fstream& outFile)
{
    outFile << "Internal, External, then Total Connections per Node:" << endl;

    for (int i = 0; i < 64; i++)
    {
        //int intConns = arrNodes[i].get_numIntConns();
        //int extConns = arrNodes[i].get_numExtConns();
        //int totConns = intConns + extConns;
        outFile << "Node " << setw(2) << i << ": " << setw(2) << endl;//intConns << ", "
            //<< setw(2) << extConns << ", " << setw(2) << totConns << endl;
    }

    outFile << endl;
} // END listNumConnections */
