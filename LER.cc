#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#include <armadillo>
#include <complex>
#include <cmath>
#include <chrono>
#include<vector>
#include<array>
#include<string>
#include<queue>


using namespace std;
using namespace arma;


void LER(const int L)
{

        queue <double>  var;
        double varia;

        vector<double> li_r_a;
        mat eigV;
        eigV.load("eigvalues.dat");

        for(size_t i=0;i<eigV.size();i++)
        {
                varia=eigV.at(i);

                if(abs(((var.size())?var.back():0.0)-varia)>0.)
                {
                 var.push(varia);
                 if(var.size()==3)
                 {
                        double var0=var.front();
                        var.pop();
                        double var1=var.front();
                        var.pop();
                        var.push(var1);
                        double var2=var.front();
                        var.pop();
                        var.push(var2);
                        double r_a=(var2-var1)/(var1-var0);
                        if(r_a<1./r_a)
                        {
                                li_r_a.push_back(r_a);
                        }
                        else
                        {
                                li_r_a.push_back(1./r_a);
                        }
                 //cout<<li_r_a.back()<<endl;

                }
                }


        }

        double sum = std::accumulate(li_r_a.begin(), li_r_a.end(), 0.0);

double mean = sum / li_r_a.size();

std::vector<double> diff(li_r_a.size());

std::transform(li_r_a.begin(), li_r_a.end(), diff.begin(), [mean](double x) { return x - mean; });

double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

double stdev = std::sqrt(sq_sum /li_r_a.size());

ofstream LER("LER");

LER<<" "<<mean<<" "<<stdev/sqrt(li_r_a.size())<<endl;
LER.close();



}

void TEO(const string str, const int column=1)
{
        int status = system(("awk '{print $"+to_string(column)+"}' "+ str+"P*" +"|tail -n +2 > .collapsed.dat").c_str());

        ifstream histdata(".collapsed.dat");
        double var,var1,var2;
        vector<double> vec;
        size_t step=0;
        while(histdata>>var)
        {
            var2=var;
            if(step)
            {
                //cout<<abs(var2-var1)<<endl;
                vec.push_back(abs(var2-var1));
            }
            else
            {
                step++;
            }
            var1=var2;
        }
        histdata.close();
        double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
        
        
double mean = sum / vec.size();

std::vector<double> diff(vec.size());

std::transform(vec.begin(), vec.end(), diff.begin(), [mean](double x) { return x - mean; });

double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

double stdev = std::sqrt(sq_sum /vec.size());

ofstream histdata2("TEO");
histdata2<<mean<<" "<<stdev/std::sqrt(vec.size())<<endl;

histdata2.close();

}

int main(int argc, char *argv[])
{
std::string str(argv[1]);
LER(1);
//int column = atoi(argv[2]);
//TEO(str, column);
}




