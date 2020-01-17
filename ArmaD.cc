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


using namespace std;
using namespace arma;


int main()
{

    size_t  L;
    cout<<"Enter L"<<endl;
    cin>>L;

    TFile *RootFile = new TFile("RootFile.root","RECREATE");

    cout<<"Creating Hamiltonian matrix"<<endl;
    auto start1 = chrono::high_resolution_clock::now();

    sp_mat P = zeros<sp_mat>(2,2);
    P(1,1)=1.;
    sp_mat X = zeros<sp_mat>(2,2);
    X(0,1)=1.;
    X(1,0)=1.;

    sp_mat ident = speye<sp_mat>(pow(2,L-3),pow(2,L-3));

    sp_mat Hamil=kron(kron(kron(X,P),ident),P) + kron(P,kron(ident,kron(P,X))) ;
    for(size_t site=0;site<L-2;site++)
    {
        sp_mat identL = speye<sp_mat>(pow(2,site),pow(2,site));
        sp_mat identR = speye<sp_mat>(pow(2,L-site-3),pow(2,L-site-3));
        Hamil += kron(identL,kron(kron(kron(P,X),P),identR));
    }

    auto start2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = start2 - start1;
    cout<<" Elapsed time Creating Hamiltonian Matrix: " << elapsed.count()<<"s"<<endl;
    cout<<"Start the diagonalization"<<endl;

    vec eigval;
    mat eigvec;

    eigs_sym(eigval, eigvec, Hamil,pow(2,L)-1);

    auto start3 = chrono::high_resolution_clock::now();
     elapsed = start3 - start2;
    cout<<" Elapsed time Diagonalizing Hamiltonian Matrix: " << elapsed.count()<<"s"<<endl;
    cout<<"writing eigenvectors"<<endl;
    eigvec.save("Eigenvectors.bin");
    cout<<"writing eigenvalues"<<endl;
    eigval.save("Eigenvalues.bin");
    //cout<<eigval<<endl;

    rowvec up={1.,0.};
    rowvec down={0.,1.};

    rowvec Z2=up;

    for(size_t j=1;j<L;j++)
    {
        Z2=kron(Z2,(j%2)?down:up);
    }
    TH2D* Prof=new TH2D("E_Vs_Log10(overlap)","",100.,-20.,20.,100,0.,-9.);
    for(size_t j=0;j<eigvec.n_cols;j++)
    {
        const double val=pow(as_scalar(Z2*eigvec.col(j)),2);
        if(log10(val)>-10.)Prof->Fill(eigval(j),log10(val));

    }

    Prof->Write();
    RootFile->Close();



}
