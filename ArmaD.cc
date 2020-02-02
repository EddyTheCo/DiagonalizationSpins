#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<TH2D.h>
#include<TFile.h>
#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#include <armadillo>
#include <complex>
#include <cmath>
#include <chrono>
#include<vector>

#include<bits/stdc++.h>

using namespace std;
using namespace arma;


vector<int> BasisMat;
vector<size_t> Rperiodicity;
int Z2_index;
int Z2_state=0;
int NumberOfstatesLessHalfChain=0;

size_t  L;
size_t kmomentum;
int refRotate=0;
int leftRotate(int n, unsigned int d)
{
    return ((n << d)|(n >> (L - d)))&(refRotate);
}

int reflect(int number)
{
    int reflect=0;
    for(size_t j=0;j<L;j++)
    {
        if(number&(1<<j))reflect=reflect|(1<<(L-1-j));
    }

    return reflect;
}
void createZ2()
{

    for(size_t i=0;i<L-1;i++)
    {
        Z2_state+=((i%2)?0:1)*static_cast<int>(pow(2,i));
    }


}



void checkstate(const int &number)
{

    int shift=number;
    for(size_t j=0;j<L;j++)
    {
        shift=leftRotate(shift,1);

        if(shift<number)
        {
            return;
        }
        else
        {
            if(shift==number)
            {
                if((kmomentum%L/(j+1))!=0)return;
                BasisMat.push_back(number);
                Rperiodicity.push_back(j+1);

                return;
            }
        }

    }
}
void CreateBasis()
{
    createZ2();
    for(size_t i=0;i<L;i++)
    {
        refRotate+=pow(2,i);
    }

    for (int i = 0; i < pow(2,L); i++)
    {

        if(!(i&(leftRotate(i,1))))
        {
            BasisMat.push_back(i);
            //checkstate(i);
        }

    }


}


int main()
{


    cout<<"Enter L and k"<<endl;
    cin>>L>>kmomentum;



    auto start0 = chrono::high_resolution_clock::now();
    CreateBasis();
/*cout<<"*************"<<endl;
    for (size_t i=0;i<BasisMat.size();i++)
    {
        cout<<BasisMat.at(i)<<endl;
    }
cout<<"*************"<<endl;*/

    auto start1 = chrono::high_resolution_clock::now();
    cout<<" Elapsed time Creating the Basis " << chrono::duration<double>(start1 - start0).count()<<"s"<<endl;
    cout<<"Creating Hamiltonian matrix"<<endl;

    const size_t BASISSIZE=BasisMat.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    mat Hamil(BASISSIZE,BASISSIZE,fill::zeros);
    for(size_t i=0;i<BASISSIZE;i++)
    {
        for(size_t j=0;j<i;j++)
        {

            if(__builtin_popcount(BasisMat.at(i)^BasisMat.at(j))==1)
            {
                Hamil(i,j)=1.;
                Hamil(j,i)=1.;
            }


        }
    }

    auto start2 = chrono::high_resolution_clock::now();

    cout<<" Elapsed time Creating Hamiltonian Matrix: " << chrono::duration<double>(start2 - start1).count()<<"s"<<endl;




   cout<<"Start the diagonalization"<<endl;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, Hamil);


    auto start3 = chrono::high_resolution_clock::now();

    cout<<" Elapsed time Diagonalizing Hamiltonian Matrix: " << chrono::duration<double>(start3 - start2).count()<<"s"<<endl;



    TFile *RootFile = new TFile("RootFile.root","RECREATE");

    TH2D* Prof=new TH2D("E_Vs_Log10(overlap)","E_Vs_Log10(overlap)",1000,-20.,20.,1000,-10.,0.);
      TH2D* Prof2=new TH2D("E_Vs_Z1","E_Vs_Z1",1000,-20.,20.,1000,-0.63,-0.27);
      TH2D* Prof3=new TH2D("E_Vs_S","E_Vs_S",1000,-20.,20.,1000,0.,7.);
      cout<<"out of E   Z1"<<endl;
      for(size_t j=0;j<eigval.size();j++)
      {

          const double val=pow(eigvec.at(BASISSIZE-1,j),2);
          double val2=0;
          double val3=0;
          for(size_t k=0;k<BASISSIZE;k++)
          {
              val2+=(__builtin_popcount(BasisMat.at(k))- (L*1.-__builtin_popcount(BasisMat.at(k))))/L*norm(eigvec.at(k,j));
          }
           // cout<<real(eigval(j))<<" "<<val2<<endl;
          mat fullBasis(1,pow(2,L));
          size_t cont=0;
          for(int f=0;f<pow(2,L);f++)
          {
              //cout<<f<<" "<<cont<<endl;
            if(f!=BasisMat.at(cont))
            {
                fullBasis.at(f)=0.;
                continue;
            }
            fullBasis.at(f)=eigvec.at(cont,j);
            if(cont<BASISSIZE-1)cont++;
          }
          //cout<<fullBasis<<endl;
          mat B = reshape(fullBasis, pow(2,L/2),pow(2,L/2));
          mat Reduced=B*B.t();


          vec eigval2;
          mat eigvec2;


          eig_sym(eigval2, eigvec2, Reduced);


          for (size_t h=0;h<eigval2.size();h++)
          {
              const double elem=eigval2.at(h);
              if(elem>0.000000000000000001)val3+=elem*log(elem);
          }
          //cout<<eigval(j)<<" "<<-val3<<endl;
          Prof->Fill(real(eigval(j)),log10(val));
          Prof2->Fill(real(eigval(j)),val2);
          //Prof3->Fill(eigval(j),-val3);


      }

    Prof->Write();
    Prof2->Write();
    Prof3->Write();
    RootFile->Close();



}




