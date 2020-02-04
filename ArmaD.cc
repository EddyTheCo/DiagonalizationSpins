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
vector<int>Mperiodicity;
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
    size_t R;
    int m=-1;
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

                R=j+1;
                break;
            }
        }

    }

    int reflected=reflect(number);
    for (size_t j=0;j<R;j++)
    {
        if(reflected<number)
        {
            return;
        }
        else
        {
            if(reflected==number)
            {
                m=j;
                break;
            }
        }
        reflected=leftRotate(reflected,1);
    }


    BasisMat.push_back(number);
    Rperiodicity.push_back(R);
    Mperiodicity.push_back(m);
    //cout<<"pushing "<<number<<" "<<R<<" "<<m<<endl;

    return;
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

            checkstate(i);
        }

    }


}

void representative(const int &number,int &rep,size_t &l)
{
    rep=number;
    int t=number;
    for(size_t i=1;i<L;i++)
    {
        t=leftRotate(t,1);
        if(t<rep)
        {
            rep=t;
            l=i;
        }

    }
     t=reflect(number);
    for(size_t i=0;i<L;i++)
    {

        if(t<rep)
        {
            rep=t;
            l=i;
        }
        t=leftRotate(t,1);
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
        cout<<BasisMat.at(i)<<" "<<Rperiodicity.at(i)<<" "<<Mperiodicity.at(i)<<endl;
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
        int five=5;

        for(size_t pos=0;pos<L;pos++)
        {
            five=leftRotate(five,1);
            int number=BasisMat.at(i);
            //cout<<"number=   "<<number<<endl;
            if(!(number&(five)))
            {
                number ^= 1UL << ((pos+2)%L);
                //cout<<"number2=   "<<number<<endl;
                int rep;
                size_t l=0;
                representative(number,rep,l);

                //cout<<"rep=     "<<rep<<endl;

                std::vector<int>::iterator it = std::find(BasisMat.begin(), BasisMat.end(), rep);
                if (it != BasisMat.end())
                {
                    int index = std::distance(BasisMat.begin(), it);
                    //cout<<"index=  "<<index<<endl;
                    double val=((Mperiodicity.at(i)==-1)?2.:1.)/((Mperiodicity.at(index)==-1)?2.:1.);
                    Hamil(i,index)+=sqrt(Rperiodicity.at(i)*val/Rperiodicity.at(index));

                }

            }

        }
    }
    //cout<<Hamil<<endl;

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
      //cout<<"out of E   Z1"<<endl;
      for(size_t j=0;j<eigval.size();j++)
      {

          const double val=eigvec.at(BASISSIZE-1,j);
          double val2=0;
          double val3=0;
          for(size_t k=0;k<BASISSIZE;k++)
          {
              val2+=(2.*__builtin_popcount(BasisMat.at(k))-L)*1.*pow(eigvec.at(k,j),2)/L;
          }
           // cout<<real(eigval(j))<<" "<<val2<<endl;
         mat fullBasis(1,pow(2,L),fill::zeros);

          for(int f=0;f<pow(2,L);f++)
          {
              int theRep;
              size_t thel;
              representative(f,theRep,thel);
              std::vector<int>::iterator it = std::find(BasisMat.begin(), BasisMat.end(), theRep);
              if (it != BasisMat.end())
              {
                  int index = std::distance(BasisMat.begin(), it);
                  fullBasis.at(f)+=eigvec.at(index,j)*sqrt(((Mperiodicity.at(index)==-1)?0.5:1.)/Rperiodicity.at(index));
              }


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
           cout<<eigval(j)<<" "<<-val3<<endl;

          Prof->Fill(real(eigval(j)),log10(val));
          Prof2->Fill(real(eigval(j)),val2);
          Prof3->Fill(eigval(j),-val3);


      }

    Prof->Write();
    Prof2->Write();
    Prof3->Write();
    RootFile->Close();



}




