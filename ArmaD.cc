#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<TH2D.h>
#include<TFile.h>
#define ARMA_DONT_USE_WRAPPER
#include <armadillo>
#include <complex>
#include <cmath>
#include <chrono>
#include<vector>

using namespace std;
using namespace arma;


vector<int> BasisMat;
int createZ2(const size_t L)
{
 int var=0;
    for(size_t i=0;i<L;i++)
    {
        var+=((i%2)?1:0)*static_cast<int>(pow(2,i));
    }

    return var;
}


void CreateBasis(const size_t L);

int main()
{

    size_t  L;
    cout<<"Enter L"<<endl;
    cin>>L;

    TFile *RootFile = new TFile("RootFile.root","RECREATE");

    cout<<"Creating Hamiltonian matrix"<<endl;
    auto start1 = chrono::high_resolution_clock::now();

    CreateBasis(L);


    const size_t BASISSIZE=BasisMat.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    mat Hamil(BASISSIZE,BASISSIZE);


    for(size_t i=0;i<BASISSIZE;i++)
    {
        for(size_t j=0;j<BASISSIZE;j++)
        {
            int var=(BasisMat.at(i))^(BasisMat.at(j));
            Hamil.at(i,j)=__builtin_popcount(var)/(L*1.0);

        }
    }

    //cout<<"HAMILTONIAN"<<endl;
    //cout<<Hamil<<endl;
    auto start2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = start2 - start1;
    cout<<" Elapsed time Creating Hamiltonian Matrix: " << elapsed.count()<<"s"<<endl;
    cout<<"Start the diagonalization"<<endl;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, Hamil);

    auto start3 = chrono::high_resolution_clock::now();
    elapsed = start3 - start2;
    cout<<" Elapsed time Diagonalizing Hamiltonian Matrix: " << elapsed.count()<<"s"<<endl;
    cout<<"writing eigenvectors"<<endl;
    eigvec.save("Eigenvectors.bin");
    cout<<"writing eigenvalues"<<endl;
    eigval.save("Eigenvalues.bin");

    //cout<<"eigenvalues="<<endl;
    //cout<<eigval<<endl;
    //cout<<"eigenvectors="<<endl;
    //cout<<eigvec<<endl;
    rowvec Z2(BASISSIZE);
    int Z_2=createZ2(L);

    for(size_t j=0;j<BASISSIZE;j++)
    {
        Z2.at(j)=__builtin_popcount((BasisMat.at(j))&Z_2)/(L*1.0);
    }

//cout<<"Z2="<<endl;
//cout<<Z2<<endl;
    TH2D* Prof=new TH2D("E_Vs_Log10(overlap)","",100,0.,0.,100,0.,0.);
    Prof->SetCanExtend(TH1::kYaxis);
    Prof->SetCanExtend(TH1::kXaxis);
    for(size_t j=0;j<BASISSIZE;j++)
    {
        const double val=pow(as_scalar(Z2*eigvec.col(j)),2);
        //cout<<"Z2*eigvect #"<<j<<" = "<<Z2*eigvec.col(j)<<endl;
        //cout<<"eige="<<eigvec.col(j);
        cout<<"log10(val)="<<log10(val)<<" eigenvalue="<<eigval(j)<<endl;
        Prof->Fill(eigval(j),log10(val));

    }

    Prof->Write();
    RootFile->Close();



}

void CreateBasis(const size_t L)
{
    BasisMat.push_back(0);
    size_t BasisMatSIZE=BasisMat.size();
cout<<"Creating Basis"<<endl;
vector<int> rep;
    for (size_t basisIndex=0;basisIndex<BasisMatSIZE;basisIndex++)
    {

        int number=BasisMat.at(basisIndex);


        for(size_t site=0;site<L;site++)
        {

            if(!(number&(1<<site)))
            {
                int number2=number|(1<<site);

                bool put=true;
                for(size_t h=0;h<rep.size();h++)
                {
                    if(number2==rep.at(h))
                    {
                        put=false;
                        break;
                    }

                }
                if(put)
                {
                    if(__builtin_popcount((number2&((1<<(L-1))+1)))>1)
                    {

                        put=false;
                    }
                    else
                    {
                        for(size_t k=0;k<L-1;k++)
                        {

                            if(__builtin_popcount(number2&(3<<k))>1)
                            {
                                put=false;
                                break;
                            }
                        }
                    }
                    if(put)
                    {

                        BasisMat.push_back(number2);
                        rep.push_back(number2);
                    }
                }
            }
        }
        if(basisIndex==BasisMatSIZE-1&&BasisMatSIZE!=BasisMat.size())
        {
            basisIndex=BasisMatSIZE-1;
            BasisMatSIZE=BasisMat.size();
            rep.clear();

        }

    }
}


