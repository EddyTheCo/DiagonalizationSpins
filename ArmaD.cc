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

#include<bits/stdc++.h>

using namespace std;
using namespace arma;




Row<int>  BasisMat;
Row<size_t> BasisRPerio;
Row<int> BasisMPerio;
vector<int> BasisMatRestrictedLHalf;

size_t Read(const string str)
{
    size_t var;
    cout<<"input "<<str<<endl;
    cin>>var;
    return var;
}
const size_t  L=Read("L");
const size_t kmomentum=Read("k");
int RRotate(void)
{
    int val=0;
    for(size_t i=0;i<L;i++)
    {
        val+=pow(2,i);
    }
    return val;

}
const int refRotate=RRotate();
int leftRotate(const int &n, const size_t &d)
{
    return ((n << d)|(n >> (L - d)))&(refRotate);
}

int reflect(const int &number)
{
    int reflect=0;
    for(size_t j=0;j<L;j++)
    {
        if(number&(1<<j))reflect=reflect|(1<<(L-1-j));
    }

    return reflect;
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

    const size_t sz = BasisMat.size();
        BasisMat.resize(sz+1);
        BasisMat(sz) = number;
        BasisRPerio.resize(sz+1);
        BasisRPerio(sz) = R;
        BasisMPerio.resize(sz+1);
        BasisMPerio(sz) = m;



    return;
}
void CreateBasis()
{

    for (int i = 0; i < pow(2,L); i++)
    {

        if(!(i&(leftRotate(i,1))))
        {
            checkstate(i);
            if(i<pow(2,L/2))
            {
                BasisMatRestrictedLHalf.push_back(i);
            }

        }

    }


}

void representative(const int &number,int &rep)
{
    rep=number;
    int t=number;
    for(size_t i=1;i<L;i++)
    {
        t=leftRotate(t,1);
        if(t<rep)
        {
            rep=t;

        }

    }
     t=reflect(number);
    for(size_t i=0;i<L;i++)
    {

        if(t<rep)
        {
            rep=t;

        }
        t=leftRotate(t,1);
    }

}
int main()
{

    ofstream outData("outData",std::ofstream::out);

    auto start0 = chrono::high_resolution_clock::now();
    CreateBasis();




    auto start1 = chrono::high_resolution_clock::now();
    cout<<" Elapsed time Creating the Basis " << chrono::duration<double>(start1 - start0).count()<<"s"<<endl;
    cout<<"Creating Hamiltonian matrix"<<endl;

    const size_t BASISSIZE=BasisMat.size();
    const size_t RESTBASISLHALFSIZE=BasisMatRestrictedLHalf.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    mat Hamil(BASISSIZE,BASISSIZE,fill::zeros);
    int five=5;

    for(size_t pos=0;pos<L;pos++)
    {
        five=leftRotate(five,1);


        BasisMat.for_each( [five,pos,&Hamil](const int number) {

            if(!(number&(five)))
            {

                const int varnumber =number^ 1UL << ((pos+2)%L);
                int rep;
                representative(varnumber,rep);

                uvec q1 = find(BasisMat == rep,1);
                uvec q2 = find(BasisMat == number,1);

                Hamil.at(q2.at(0),q1.at(0))+=sqrt(1.*BasisRPerio.at(q2.at(0))*((BasisMPerio.at(q2.at(0))==-1)?2.:1.)
                                                  /(BasisRPerio.at(q1.at(0))*((BasisMPerio.at(q1.at(0))==-1)?2.:1.)));


            }


            } );





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


    eigval.save("eigvalues.dat");
    outData<< left << setw(12)<<"Energy "<< left << setw(12) << "Z1 " << left << setw(12) <<"Z4 " <<left << setw(12) << "X1 " << left << setw(12) <<"S " <<endl;



      for(size_t j=0;j<eigval.size();j++)
      {
//cout<<"WORKING ON EIGVECTOR "<<eigval.at(j)<<endl;

          //const double val=eigvec.at(BASISSIZE-1,j)/sqrt(2.);

          double Xval=0;

          double val3=0;
          const Col<double> vectNORM=pow(eigvec.col(j),2);
          Row<int> Zv=BasisMat;
                  Zv.for_each( []( int &number) {
                number= (2*__builtin_popcount(number)-L);
              } );


          const double Zval=dot(Zv,vectNORM);


          Row<int> Zv4=BasisMat;
                  Zv4.for_each( []( int &number) {
                      number=number^leftRotate(number,4);
                number= (L-2*__builtin_popcount(number));
              } );
          const double Z4val=dot(Zv4,vectNORM);


          for(size_t l=0;l<L;l++)
          {

              BasisMat.for_each( [l,vectNORM,&Xval](const int &number) {

                      const int varnumber =number^ 1UL << l;
                      if(!(varnumber&(leftRotate(varnumber,1))))
                      {
                          int rep;
                          representative(varnumber,rep);
                          uvec q1 = find(BasisMat == rep,1);
                          Xval+=vectNORM.at(q1.at(0));

                      }
                  } );

          }


          /*
          mat fullBasis(RESTBASISLHALFSIZE,RESTBASISLHALFSIZE,fill::zeros);
          size_t irow=0;

          for(vector<int>::iterator  itR=BasisMatRestrictedLHalf.begin();itR!= BasisMatRestrictedLHalf.end();itR++)
          {
              size_t icolumn=0;
              for(vector<int>::iterator  itL=BasisMatRestrictedLHalf.begin();itL!= BasisMatRestrictedLHalf.end();itL++)
              {

                    const int num=((*itR)|leftRotate(*itL,L/2));

                    if(!(num&(leftRotate(num,1))))
                    {

                        int theRep;
                        size_t thel;
                        representative(num,theRep,thel);

                        std::vector<BasNumber>::iterator it = std::lower_bound(BasisMat.begin(), BasisMat.end(), theRep,ComparStruct);

                        int index = std::distance(BasisMat.begin(), it);
                        fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                   sqrt(((it->Mperio==-1)?0.5:1.)/it->Rperio);


                    }


                    icolumn++;

              }

              irow++;

          }


           mat Reduced=fullBasis*fullBasis.t();


           vec eigval2;
           mat eigvec2;


           eig_sym(eigval2, eigvec2, Reduced);
           val3=trace(eigval2.t()*trunc_log(eigval2));

*/
           outData<< left << setw(12)<<eigval(j)<< left << setw(12) << Zval<< left << setw(12) <<Z4val<< left << setw(12) << Xval<< left << setw(12) <<-val3<<endl;







      }

      auto start4 = chrono::high_resolution_clock::now();

      cout<<" Elapsed time working on eigenvectors: " << chrono::duration<double>(start4 - start3).count()<<"s"<<endl;

    outData.close();



}




