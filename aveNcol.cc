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



void aveNCol(string file,const size_t NCol, const size_t colErr=0)
{



        ifstream dataFile(file);
        ofstream OutFile("AVE"+file);
        vector<vector<double>> arr;
        vector<double> SUMERR;
        double val;
        size_t step=0;

        while (dataFile>>val)
        {
            if(colErr!=0&&(step%NCol==colErr-1))
            {
                SUMERR.push_back(val);
            }
            if(step/NCol==0)
            {
                vector<double> var;
                var.push_back(val);
                arr.push_back(var);

            }
            else
            {
                arr.at(step%NCol).push_back(val);
            }

            step++;
        }

        double sq_sum1 = std::inner_product(SUMERR.begin(), SUMERR.end(), SUMERR.begin(), 0.0);

        double stdev1 = std::sqrt(sq_sum1 )/SUMERR.size();


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

OutFile<<stdev1<<endl;

OutFile.close();
dataFile.close();

}

int main(int argc, char *argv[])
{

    int column = atoi(argv[2]);
    std::string str(argv[1]);


    if(argc>3)
    {
        int column2 = atoi(argv[3]);
        aveNCol(str,column,column2);
    }
    else
    {
        aveNCol(str,column);
    }

}
