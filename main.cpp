#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <string>
#include <sstream>

using namespace std;

int n=32; //Grid size
double T=2.35;

//Declare matrix for the spins
vector<vector<int> > grid; //Stores -1 or 1

//Mod function to set boundary conditions
double mod(int a, int N)
{
    return (a%N +N)%N;
}

//Flips the spin
void flip(int x, int y)
{
    grid.at(x).at(y)*=-1;
}

//Calculates the system's energy
double calculateE()
{
    int E=0;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            E+= -1*(grid.at(i).at(j)*grid.at(i).at(mod(j+1,n)) + grid.at(i).at(j)*grid.at(i).at(mod(j-1,n)) + grid.at(i).at(j)*grid.at(mod(i+1,n)).at(j) + grid.at(i).at(j)*grid.at(mod(i-1,n)).at(j));
        }
    }

    return E/2; //Divide so as not to count each bond two times
}

//Calculates the system's magnetization
double calculateM()
{
    int M=0;

    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            M+= grid.at(i).at(j);
        }
    }

    return M;
}

//Calculates the Boltzmann factor for the spin
double boltzmann(int x, int y)
{
    double dE = (grid.at(x).at(y)*grid.at(x).at(mod(y+1,n)) + grid.at(x).at(y)*grid.at(x).at(mod(y-1,n)) + grid.at(x).at(y)*grid.at(mod(x+1,n)).at(y) + grid.at(x).at(y)*grid.at(mod(x-1,n)).at(y)) + (grid.at(x).at(y)*grid.at(x).at(mod(y+1,n)) + grid.at(x).at(y)*grid.at(x).at(mod(y-1,n)) + grid.at(x).at(y)*grid.at(mod(x+1,n)).at(y) + grid.at(x).at(y)*grid.at(mod(x-1,n)).at(y));
    double p=exp(-dE/T);

    return p;
}

int main()
{
    srand(time(NULL));

    vector<double> ExSpin;
    vector<double> MxSpin;

    //Initialize grid with random spins
    grid.resize(n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            if(rand()%2==0)
            {
                grid.at(i).push_back(-1);
            }
            else
            {
                grid.at(i).push_back(1);
            }
        }
    }

    //Save initial energy and magnetization
    ExSpin.push_back(calculateE()/(n*n));
    MxSpin.push_back(calculateM()/(n*n));

    //Time cycle
    for(int i=1; i<n*n*500; i++)
    {
        //If this iteration completes a time step (n*n iterations), save energy and magnetization
        if(i%(n*n)==0)
        {
            ExSpin.push_back(calculateE()/(n*n));
            MxSpin.push_back(calculateM()/(n*n));
        }

        //Choose random spin and w between 0 and 1
        int x=rand()%n;
        int y=rand()%n;
        double w=1.0*rand()/RAND_MAX;

        double p=boltzmann(x,y);

        //Flips the spin if it corresponds
        if(w<=p)
        {
            flip(x,y);
        }
    }

    //Open text file
    ofstream oFile;
    oFile.open("T2,35.dat");
    oFile.is_open();
    oFile<< "t " << "E "<<"M"<< endl;

    //Save data
    for(int i=0; i<ExSpin.size(); i++)
    {
        oFile<<i<<" "<<ExSpin.at(i)<<" "<<MxSpin.at(i)<<'\n';
    }

    return 0;
}
