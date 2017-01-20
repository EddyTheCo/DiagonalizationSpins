#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<vector>
#include <random>
#include"constants.hh"
#include "site.hh"
#include"myTGraph.hh"
#include<stack>
#include<string>
using namespace std;


long int Number_equi, N_sweeps, N_flips;

int main()
{
void crearMuestra(int , Site**,double);
double makeMeasur(Site**,double);
double calculateStd( stack<double> &,double);
Site** lattice;
int l,binsconst,lfinal,lstep,lini;
double aver,T;
vector <double> theMeas,theLeng,theMeasErr;
stack <double> errorSamp;


string INIMacro="INIMacro";
ifstream theMacro (INIMacro.c_str());

    theMacro>>T;
    theMacro>>Number_equi;
    theMacro>>N_sweeps;
    theMacro>>N_flips;
    theMacro>>binsconst;
    theMacro>>lini;
    theMacro>>lfinal;
    theMacro>>lstep;

    char title[20];
    sprintf(title,"Temp= %.2f ",T);

    string output="OUTPUT";
    ofstream theoutput (output.c_str());

    for ( l=lini;l>lfinal;l-=lstep)
    {
        aver=0;
        int cont=0;
        for(int bins=0;bins<binsconst;bins++)
        {
            lattice=new Site*[l];
            for(int i=0;i<l;i++)
            {
                lattice[i]=new Site[l];
                for(int j=0;j<l;j++)
                {
                    lattice[i][j]=Site(i,j);
                }
            }

            crearMuestra(l,lattice,T);

            double samp=makeMeasur(lattice,l);

        }

            theMeas.push_back(aver/cont);
            double therror=calculateStd(errorSamp,aver/cont);
            theMeasErr.push_back(therror);
            theLeng.push_back(l);

            theoutput<<l<<" "<<aver/cont<<" "<<therror<<endl;


    }
    MyTGraph * he=new MyTGraph(theMeas.size(),&theLeng[0],&theMeas[0],0,&theMeasErr[0],title,T);







}


void crearMuestra(int l,Site** lattice,double T)
{
    uniform_int_distribution<int> dist2(0,l-1);

    lattice[0][0].setLattice(lattice,l,T);

    lattice[0][0].TEnergy=0;
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<l;j++)
        {
            lattice[i][j].setPtr();
            lattice[i][j].totalEnergy();

        }
    }
    //cout<<lattice[0][0].TEnergy<<endl;
    for(int steps=0;steps<Number_equi;steps++)
    {
        lattice[dist2(Constants::mt)][dist2(Constants::mt)].changeConf();
    }

}

double makeMeasur(Site** lattice,double l)
{
    vector <double> energy;
    uniform_int_distribution<int> dist2(0,l-1);
     int N_meas=N_sweeps/N_flips;
    double E=0;


    for( long int step=0;step<N_sweeps;step++)
    {
        lattice[dist2(Constants::mt)][dist2(Constants::mt)].changeConf();
        if(step%N_flips==0)
        {
            E+=lattice[0][0].TEnergy;
            energy.push_back(lattice[0][0].TEnergy);
        }
    }
    E=E/(N_meas);
//cout<<E<<endl;



    for(int k=0;k<N_meas;k++)
    {
        energy[k]-=E;
    }
    double gamma=0;
    double gamma2=1234;
    for(int j=0;j<N_meas;j++)
    {
        double varA2=0,varA=0;
        for(int i=0;i<N_meas-j;i++)
        {
            varA2+=pow(energy[i],2);
            varA+=energy[i]*energy[i+j];
        }
        gamma+=varA/varA2;

        if(abs(gamma2-gamma)>0.01)
        {
            gamma2=gamma;
        }
        else
        {
            //cout<<gamma2<<endl;
            break;
        }
     //cout<<j<<" "<<varA/varA2<<" "<<gamma<<endl;

    }


    return gamma+0.5;

}
double calculateStd(stack<double> & vect,double ave)
{
    double sum=0;
    int k=vect.size();
    while(!vect.empty())
    {
        sum+=pow((vect.top()-ave),2);

        vect.pop();
    }

    return sqrt(sum/k);

}
