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



struct BasNumber
{
     int basi;
    size_t Rperio;
    int Mperio;
};

bool ComparStruct(const BasNumber &a,const  int &b) {
   return(a.basi<b);
}
vector<BasNumber> BasisMat;

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

    const BasNumber a = { .basi = number, .Rperio = R, .Mperio = m };

    BasisMat.push_back(a);

    return;
}
void CreateBasis()
{

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

    ofstream outData("outData",std::ofstream::out);

    auto start0 = chrono::high_resolution_clock::now();
    CreateBasis();

    auto start1 = chrono::high_resolution_clock::now();
    cout<<" Elapsed time Creating the Basis " << chrono::duration<double>(start1 - start0).count()<<"s"<<endl;
    cout<<"Creating Hamiltonian matrix"<<endl;

    const size_t BASISSIZE=BasisMat.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    mat Hamil(BASISSIZE,BASISSIZE,fill::zeros);
    for(vector<BasNumber>::iterator  it=BasisMat.begin();it!= BasisMat.end();it++)
    {
        int five=5;

        for(size_t pos=0;pos<L;pos++)
        {
            five=leftRotate(five,1);
            int number=it->basi;

            if(!(number&(five)))
            {
                number ^= 1UL << ((pos+2)%L);

                int rep;
                size_t l=0;
                representative(number,rep,l);



                std::vector<BasNumber>::iterator it2 = std::lower_bound(BasisMat.begin(), BasisMat.end(), rep,ComparStruct);
                if (it2 != BasisMat.end())
                {
                    const int index = std::distance(BasisMat.begin(), it2);
                    const int i = std::distance(BasisMat.begin(), it);

                    const double val=((it->Mperio==-1)?2.:1.)/((it2->Mperio==-1)?2.:1.);
                    Hamil(i,index)+=sqrt(it->Rperio*val/it2->Rperio);

                }
                else
                {
                    cout<<"Error1"<<endl;
                }

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


    eigval.save("eigvalues.dat");
    outData<< left << setw(12)<<"Energy "<< left << setw(12) << "Z1 " << left << setw(12) <<"Z4 " <<left << setw(12) << "X1 " << left << setw(12) <<"S " <<endl;



      for(size_t j=0;j<eigval.size();j++)
      {

          //const double val=eigvec.at(BASISSIZE-1,j)/sqrt(2.);
          double Zval=0;
          double Xval=0;
          double Z4val=0;
          double val3=0;
          for(size_t k=0;k<BASISSIZE;k++)
          {
              const double vectornorm=pow(eigvec.at(k,j),2);
              int numero=BasisMat.at(k).basi;
              Zval+=(2.*__builtin_popcount(numero)-L)*1.*vectornorm;
              numero=(numero^leftRotate(numero,4));
              const int nup=  __builtin_popcount(numero);
              Z4val+=(L-2.*nup)*vectornorm;

                    for(size_t l=0;l<L;l++)
                    {
                        int nu=BasisMat.at(k).basi;
                        nu ^= 1UL << l;
                        if(!(nu&(leftRotate(nu,1))))
                        {
                            int rep;
                            size_t l1=0;
                            representative(nu,rep,l1);
                            std::vector<BasNumber>::iterator it = std::lower_bound(BasisMat.begin(), BasisMat.end(), rep,ComparStruct);
                            if (it != BasisMat.end())
                            {
                                const int index = std::distance(BasisMat.begin(), it);
                                Xval+=eigvec.at(index,j)*eigvec.at(k,j);
                            }
                            else
                            {
                                cout<<"Error1"<<endl;
                            }
                        }

                    }

          }

         const size_t maxRow=BasisMat.back().basi/pow(2,L/2);
         const size_t maxcolumn=(BasisMat.back().basi)% static_cast<size_t>(pow(2,L/2));
         mat fullBasis(maxRow+1,maxcolumn+1,fill::zeros);

          for(size_t f=0;f<BASISSIZE;f++)
          {
              const size_t Row=BasisMat.at(f).basi/pow(2,L/2);
              const size_t Column=BasisMat.at(f).basi% static_cast<size_t>(pow(2,L/2));
              fullBasis.at(Row,Column)=eigvec.at(f,j);
          }

          mat Reduced=fullBasis*fullBasis.t();
          vec eigval2;
          mat eigvec2;


          eig_sym(eigval2, eigvec2, Reduced);


          for (size_t h=0;h<eigval2.size();h++)
          {
              const double elem=eigval2.at(h);
              if(elem>0.000000000000000001)val3+=elem*log(elem);
          }
            outData<< left << setw(12)<<eigval(j)<< left << setw(12) << Zval<< left << setw(12) <<Z4val<< left << setw(12) << Xval<< left << setw(12) <<-val3<<endl;


          //cout<<eigval(j)<<" "<<-val3<<endl;




      }


    outData.close();



}




