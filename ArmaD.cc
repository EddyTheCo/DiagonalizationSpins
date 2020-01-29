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
int Z2_index;
int NumberOfstatesLessHalfChain=2;
int createZ2(const int L)
{
 int var=0;
    for(int i=0;i<L-1;i++)
    {
        var+=((i%2)?0:1)*static_cast<int>(pow(2,i));
    }

    return var;
}

void CreateBasis(const int L);

int main()
{

    int  L;
    cout<<"Enter L"<<endl;
    cin>>L;



    auto start0 = chrono::high_resolution_clock::now();
    CreateBasis(L);

    auto start1 = chrono::high_resolution_clock::now();
    cout<<" Elapsed time Creating the Basis " << chrono::duration<double>(start1 - start0).count()<<"s"<<endl;
    cout<<"Creating Hamiltonian matrix"<<endl;

    const size_t BASISSIZE=BasisMat.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    mat Hamil(BASISSIZE,BASISSIZE);
    for(size_t i=0;i<BASISSIZE;i++)
    {
        for(size_t j=0;j<i;j++)
        {

            int var=(BasisMat.at(i))^(BasisMat.at(j));
            if(__builtin_popcount(var)==1)
            {
                Hamil.at(i,j)=1.;
                Hamil.at(j,i)=1.;
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
      for(size_t j=0;j<eigval.size();j++)
      {

          const double val=pow(eigvec.at(Z2_index,j),2);
          double val2=0;
          double val3=0;
          for(size_t k=0;k<BASISSIZE;k++)
          {
              val2+=(__builtin_popcount(BasisMat.at(k))- (L*1.-__builtin_popcount(BasisMat.at(k))))/L*pow(eigvec.at(k,j),2);
          }
          mat B = reshape(eigvec.col(j), NumberOfstatesLessHalfChain,BASISSIZE/NumberOfstatesLessHalfChain);
          mat Reduced=B*B.t();

          vec eigval2;
          mat eigvec2;

          eig_sym(eigval2, eigvec2, Reduced);

          for (size_t h=0;h<eigval2.size();h++)
          {
              const double elem=eigval2.at(h);
              if(elem>0.0000000000001)val3+=elem*log(elem);
          }

          Prof->Fill(eigval(j),log10(val));
          Prof2->Fill(eigval(j),val2);
          Prof3->Fill(eigval(j),-val3);


      }

    Prof->Write();
    Prof2->Write();
    Prof3->Write();
    RootFile->Close();



}

void recur(const int &number, const int &L ,const int &k,const int &Z_2)
{
    static bool var=false;
    if(k>=L-1)return;
    for(size_t m=(2+k);m<L-1;m++)
    {
        const int number2=number|(1<<(m));

        if(number2>pow(2,L/2))
        {
            BasisMat.push_back(number2);
            if(number2==Z_2)
            {
                var=true;
                Z2_index=BasisMat.size()-1;

            }
        }
        else
        {
            NumberOfstatesLessHalfChain++;
            BasisMat.insert(BasisMat.begin(),number2);
            if(var)
            {
                Z2_index++;
            }
        }




        recur(number2,L,m,Z_2);
    }




}
void CreateBasis(const int L)
{
    int Z_2=createZ2(L);
    BasisMat.push_back(0);
    BasisMat.push_back(1);

        recur(1,L,0,Z_2);


}


