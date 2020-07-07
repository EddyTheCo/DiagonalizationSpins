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



void aveNCol(string file,const size_t NCol)
{



        ifstream dataFile(file);
        ofstream OutFile("AVE"+file);
        vector<vector<double>> arr;
        double val;
        size_t step=0;
        while (dataFile>>val)
        {
            arr.at(step%NCol).at(step/NCol)=val;
            step++;
        }


        for(size_t i=0;i<NCol;i++)
        {
            vector<double> li_r_a=arr.at(i);

            double sum = std::accumulate(li_r_a.begin(), li_r_a.end(), 0.0);

            double mean = sum / li_r_a.size();

            std::vector<double> diff(li_r_a.size());

            std::transform(li_r_a.begin(), li_r_a.end(), diff.begin(), [mean](double x) { return x - mean; });

            double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

            double stdev = std::sqrt(sq_sum /li_r_a.size());

            OutFile<<mean<<" "<<stdev<<endl;
        }



OutFile.close();
dataFile.close();

}

int main(int argc, char *argv[])
{

    int column = atoi(argv[2]);
    std::string str(argv[1]);
    aveNCol(str,column);

}
