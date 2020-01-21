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
int Z2_index;
int createZ2(const int L)
{
 int var=0;
    for(int i=0;i<L;i++)
    {
        var+=((i%2)?1:0)*static_cast<int>(pow(2,i));
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
    mat Hamil(BASISSIZE,BASISSIZE,fill::zeros);
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



   /* mat Hamil(BASISSIZE,BASISSIZE);

    vector<int>::iterator itini = BasisMat.begin()+1;
    const vector<int>::iterator last = BasisMat.end();
    for(size_t i=1;i<BASISSIZE;i++)
    {

            const auto itHmil=Hamil.submat( i-1, i, i-1, BASISSIZE-1 ).begin();
           const vector<int>::iterator varIndex = itini-1;
            transform(itini, last,itHmil, [varIndex] (const int& in1) {

                return ((__builtin_popcount(in1^(*varIndex))==1)?1.0:0.0); });
            itini++;

    }

    Hamil=symmatu(Hamil);



    //cout<<Hamil<<endl;
    auto start3 = chrono::high_resolution_clock::now();
   chrono::duration<double> elapsed = start3 - start2;
    cout<<" Elapsed time Creating Hamiltonian Matrix2: " << elapsed.count()<<"s"<<endl;
    */

   cout<<"Start the diagonalization"<<endl;

    vec eigval;
    mat eigvec;

    eig_sym(eigval, eigvec, Hamil);


    auto start3 = chrono::high_resolution_clock::now();

    cout<<" Elapsed time Diagonalizing Hamiltonian Matrix: " << chrono::duration<double>(start3 - start2).count()<<"s"<<endl;
    cout<<"writing eigenvectors"<<endl;
    eigvec.save("Eigenvectors.bin");
    cout<<"writing eigenvalues"<<endl;
    eigval.save("Eigenvalues.bin");


    TFile *RootFile = new TFile("RootFile.root","RECREATE");

    TH2D* Prof=new TH2D("E_Vs_Log10(overlap)","",1000,-20,20,100,-10,0);


    for(size_t j=0;j<BASISSIZE;j++)
    {

        const double val=pow(as_scalar(eigvec.col(j).row(Z2_index)),2);
        //cout<<"Z2*eigvect #"<<j<<" = "<<Z2*eigvec.col(j)<<endl;
        //cout<<"eige="<<eigvec.col(j);
        //cout<<"log10(val)="<<log10(val)<<" eigenvalue="<<eigval(j)<<endl;
        Prof->Fill(eigval(j),log10(val));

    }

    Prof->Write();
    RootFile->Close();



}
void recur(const int &number, const int &L ,const int &k,const int &Z_2)
{

    int number2=number|(1<<k);

    if(!(k==L-1&&(number&1)))
    {
        BasisMat.push_back(number2);

        if(number2==Z_2)
        {
            Z2_index=BasisMat.size()-1;
        }
    }
    if(k>L-1)return;
    for(int i=k+2;i<L;i++)
    {
        recur(number2,L,i,Z_2);
    }
}
void CreateBasis(const int L)
{
    int Z_2=createZ2(L);
    BasisMat.push_back(0);
    int number=0;


    for(int site=0;site<L;site++)
    {
        recur(number,L,site,Z_2);

    }
}


