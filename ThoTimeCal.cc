#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <numeric>
#include<vector>
#include<array>
#include<string>
#include <complex>
#include <algorithm>



using namespace std;



void ThouTimeCal(string file)
{



        ifstream dataFile(file);
        ofstream OutFile("TIME"+file);
        double valx,valy,valz;

        double pi=3.14159265358979323846;
        while (dataFile>>valx>>valy>>valz)
        {
            valx=valx/2/pi;
            if(abs(log10(valy/(2.0*valx - valx*log(1.0 + 2.0*valx ))))<0.05)
            {
                OutFile<<valx<<endl;
                break;
            }
        }




OutFile.close();
dataFile.close();

}

int main(int argc, char *argv[])
{


    std::string str(argv[1]);



        ThouTimeCal(str);


}
