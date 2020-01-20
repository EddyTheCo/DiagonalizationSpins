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
#include<bits/stdc++.h>

using namespace std;
using namespace arma;


vector<int> BasisMat;
int createZ2(const int L)
{
 int var=0;
    for(int i=0;i<L;i++)
    {
        var+=((i%2)?1:0)*static_cast<int>(pow(2,i));
    }

    return var;
}

int findDisContinuous1(int n,const int L)
{
    //cout<<"Number="<<n<<endl;
    n=~n;
    int count=0;



    if(((n&((1<<(L-2))+1))==(1<<(L-2))+1)&&((n&(1<<(L-1)))!=(1<<(L-1))))
    {
        count++;
    }
    //cout<<"Count1="<<count<<endl;
    if(((n&((1<<(L-1))+2))==(1<<(L-1))+2)&&((n&(1))!=1))
    {
        count++;
    }
    //cout<<"Count2="<<count<<endl;
        for(int k=0;k<L-2;k++)
        {

            if(((n&(5<<k))==(5<<k))&&((n&(2<<k))!=(2<<k)))
            {
                count++;
            }
        }
   // cout<<"Count3="<<count<<endl;
//cout<<"Count="<<count<<endl;

    return count;
}

void CreateBasis(const int L);

int main()
{

    int  L;
    cout<<"Enter L"<<endl;
    cin>>L;

    TFile *RootFile = new TFile("RootFile.root","RECREATE");

    cout<<"Creating Hamiltonian matrix"<<endl;
    auto start1 = chrono::high_resolution_clock::now();

    CreateBasis(L);


    const int BASISSIZE=BasisMat.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    mat Hamil(BASISSIZE,BASISSIZE);

    for(int i=0;i<BASISSIZE;i++)
    {
        for(int j=0;j<=i;j++)
        {
            int var=(BasisMat.at(i))^(BasisMat.at(j));
            //cout<<"Calculating "<<BasisMat.at(i)<<" xor "<<BasisMat.at(j)<<" = "<<((BasisMat.at(i))^(BasisMat.at(j)))<<endl;
            if(var!=0)var=findDisContinuous1(var,L);
           // cout<<"Hamil="<<var<<endl;
            Hamil.at(i,j)=var/(L*1.);
            Hamil.at(j,i)=var/(L*1.);

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


    rowvec Z2(BASISSIZE);
    int Z_2=createZ2(L);
    int NOTZ2=~Z_2;
    for(int j=0;j<BASISSIZE;j++)
    {
        Z2.at(j)=__builtin_popcount((BasisMat.at(j))^NOTZ2)/(L*1.);
    }


    TH2D* Prof=new TH2D("E_Vs_Log10(overlap)","",1000,-500,500,100,-9,0);


    for(int j=0;j<BASISSIZE;j++)
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

void CreateBasis(const int L)
{
    BasisMat.push_back(0);
    int BasisMatSIZE=BasisMat.size();
cout<<"Creating Basis"<<endl;
vector<int> rep;
    for (int basisIndex=0;basisIndex<BasisMatSIZE;basisIndex++)
    {

        int number=BasisMat.at(basisIndex);


        for(int site=0;site<L;site++)
        {

            if(!(number&(1<<site)))
            {
                int number2=number|(1<<site);

                bool put=true;
                for(int h=0;h<rep.size();h++)
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
                        for(int k=0;k<L-1;k++)
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


