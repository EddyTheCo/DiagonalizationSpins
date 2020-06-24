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
#include<random>

using namespace std;
using namespace arma;

#define PXP
#define LENGHT 12
#define INTTYPE int
#define OBC
#define W 1

const string ap="FULL";



#include<bits/stdc++.h>
//#include <bitset>


std::random_device rd;
std::mt19937::result_type seed = rd() ^ (
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::system_clock::now().time_since_epoch()
                ).count() +
            (std::mt19937::result_type)
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now().time_since_epoch()
                ).count() );
std::mt19937 mt(seed);
uniform_real_distribution<double> dist(-W,W);
void saveRandom(void)
{

        ofstream fout(".seed.dat");
        fout<<mt;
        fout.close();
}

     INTTYPE basi;

#ifdef OBC
#define LEFTROT(a,b) leftRotateOBC(a,b)
#else
#define LEFTROT(a,b) leftRotate(a,b)
#endif

#ifdef PXP
INTTYPE five=5;
#define NOFPS 1
#define CHECKUPS(i) (i&(LEFTROT(i,1)))
#endif
#ifdef PPXPP
INTTYPE five=27;
#define CHECKUPS(i) ((i&(LEFTROT(i,1)))||(i&(LEFTROT(i,2))))
#define NOFPS 2
#endif
#ifdef PPPXPPP
INTTYPE five=119;
#define CHECKUPS(i) ((i&(LEFTROT(i,1)))||(i&(LEFTROT(i,2)))||(i&(LEFTROT(i,3))))
#define NOFPS 3
#endif
#ifdef PPPPXPPPP
INTTYPE five=495;
#define CHECKUPS(i) ((i&(LEFTROT(i,1)))||(i&(LEFTROT(i,2)))||(i&(LEFTROT(i,3)))||(i&(LEFTROT(i,4))))
#define NOFPS 4
#endif
#ifdef PPPPPXPPPPP
INTTYPE five=2015;
#define CHECKUPS(i) ((i&(LEFTROT(i,1)))||(i&(LEFTROT(i,2)))||(i&(LEFTROT(i,3)))||(i&(LEFTROT(i,4)))||(i&(LEFTROT(i,5))))
#define NOFPS 5
#endif

vector<INTTYPE> BasisMat;
vector<int> BasisMatRestrictedLHalf;

array<double,LENGHT> wi;

void PRINTINFO(const size_t &Ba)
{
    ofstream INFO("INFO"+ap,std::ofstream::out);
saveRandom();
#ifdef PXP
INFO<<" pxp"<<endl;
#endif
#ifdef PPXPP
INFO<<" ppxpp"<<endl;
#endif
#ifdef PPPXPPP
INFO<<" pppxppp"<<endl;
#endif
#ifdef PPPPXPPPP
INFO<<" ppppxpppp"<<endl;
#endif
#ifdef PPPPPXPPPPP
INFO<<" pppppxppppp"<<endl;
#endif
INFO<<"Disordered W="<<W<<endl;
INFO<<"L="<<LENGHT<<endl;
INFO<<"basis="<<Ba<<endl;

#ifdef OBC
INFO<<"OBC"<<endl;
#else
INFO<<"PBC"<<endl;
#endif
INFO.close();
}



INTTYPE RRotate(void)
{
    INTTYPE val=0;
    for(size_t i=0;i<LENGHT;i++)
    {
        val=val|(1ULL<<i);
    }
    return val;

}
const INTTYPE refRotate=RRotate();

INTTYPE leftRotate(const INTTYPE &n, const int &d)
{
    if(d>=0)
    {
        return (((n << d)|(n >> (LENGHT - d)))&(refRotate));
    }
    else
    {
        return (((n >> (-1*d))|(n << (LENGHT + d)))&(refRotate));
    }

}
INTTYPE leftRotateOBC(const INTTYPE &n, const int &d)
{
    if(d>=0)
    {
        return (((n << d))&(refRotate));
    }
    else
    {
        return (((n >> (-1*d)))&(refRotate));
    }

}

void CreateBasis()
{

    BasisMat.push_back(0);


    for(size_t i=0;i<LENGHT;i++)
    {
        wi.at(i)=dist(mt);
        cout<<wi.at(i)<<endl;
        const size_t BAS=BasisMat.size();


        for(size_t g=0;g<BAS;g++)
        {

            const INTTYPE num=(BasisMat.at(g))|leftRotate(1,i);


            if(!CHECKUPS(num))
            {

                BasisMat.push_back(num);


            }


        }




    }





}


int main()
{


    auto start0 = chrono::high_resolution_clock::now();


    CreateBasis();


    auto start1 = chrono::high_resolution_clock::now();
    cout<<" Elapsed time Creating the Basis " << chrono::duration<double>(start1 - start0).count()<<"s"<<endl;
    cout<<"Creating Hamiltonian matrix"<<endl;

    const size_t BASISSIZE=BasisMat.size();


    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    PRINTINFO(BASISSIZE);
    mat Hamil(BASISSIZE,BASISSIZE,fill::zeros);

    for(vector<INTTYPE>::iterator  it=BasisMat.begin();it!= BasisMat.end();it++)
    {
        for(size_t pos=0;pos<LENGHT;pos++)
        {

            INTTYPE fiveVar=LEFTROT(five,pos-NOFPS);

            INTTYPE number=*it;

            int UpDown=((number) & (1ULL<<(pos)))?1:-1;
            if(!(number&(fiveVar)))
            {

                number ^= 1ULL << pos;
                std::vector<INTTYPE>::iterator it2 = std::lower_bound(BasisMat.begin(), BasisMat.end(), number);
                if (it2 != BasisMat.end())
                {
                    const size_t index = std::distance(BasisMat.begin(), it2);
                    const size_t i = std::distance(BasisMat.begin(), it);
                    Hamil(i,index)+=1.;

                }
                else
                {
                    cout<<"Error1"<<endl;
                }

            }
            Hamil(pos,pos)+=wi.at(pos)*UpDown;



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


    eigval.save("eigvalues.dat");



}




