#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#include <armadillo>




using namespace std;
using namespace arma;

void Heisen(void)
{
    mat eigV;
    eigV.load("eigvalues.dat");

   double tra=trace(eigV);
   tra=tra*tra;

   double tra2=norm(eigV);


    ofstream OUT("HEISTIME");
    OUT<<tra2<<" "<<tra<<endl;

OUT.close();


}

int main()
{
    Heisen();
}
