#include<iostream>
#include "site.hh"
#include<random>
#include"constants.hh"
using namespace std;



uniform_int_distribution<int> dist(-1,0);
uniform_real_distribution<double> distd(0,1);

double Site::TEnergy=0;
double Site::totalMag=0;
int Site::numberOfSites=0;
int Site::L=0;
double Site::T=0;
Site** Site::theLattice=NULL;
Site::Site()
{

}

Site::Site(int i,int j):i(i),j(j)
{
    numberOfSites+=1;
    spin=dist(Constants::mt);
    if(spin==0)
        spin=1;

}
Site::~Site()
{
	
	
}
void Site::calculateEnergy()
{
    //cout<<left->spin<<" "<<right->spin<<" "<<down->spin<<" "<<up->spin<<endl;
   energy= -1*spin*(left->spin+right->spin+down->spin+up->spin);
   dE=-2*energy;
}
void Site::totalEnergy()
{

    calculateEnergy();
   // cout<<energy<<endl;
    TEnergy+=energy/2;

}
void Site::setPtr(void)
{

    left=&theLattice[(i-1+L)%L][j];
    right=&theLattice[(i+1+L)%L][j];
    down=&theLattice[i][(j-1+L)%L];
    up=&theLattice[i][(j+1+L)%L];
}
void Site::changeConf()
{


calculateEnergy();
    if(dE<=0.)
    {
        spin*=-1;
        TEnergy+=dE;


    }
    else
    {
        if(exp(-1./T*dE)>distd(Constants::mt))
        {
            spin*=-1;
            TEnergy+=dE;


        }
    }
}
