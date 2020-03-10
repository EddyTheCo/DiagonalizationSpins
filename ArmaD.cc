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


using namespace std;
using namespace arma;

#define PXP
#define LENGHT 20
#define INTTYPE int
#define USEPARITY


#ifdef USEPARITY
#define PARITY1
#ifdef PARITY1
const string ap="P1";
#else
const string ap="P2";
#endif
#else
const string ap="FULL";
#endif


#include<bits/stdc++.h>
//#include <bitset>




struct BasNumber
{
     INTTYPE basi;
    size_t Rperio;
#ifdef PARITY1
    int Mperio;
#endif
};

#ifdef PXP
INTTYPE five=5;
#define NOFPS 1
#define CHECKUPS(i) (i&(leftRotate(i,1)))
#endif
#ifdef PPXPP
INTTYPE five=27;
#define CHECKUPS(i) ((i&(leftRotate(i,1)))||(i&(leftRotate(i,2))))
#define NOFPS 2
#endif
#ifdef PPPXPPP
INTTYPE five=119;
#define CHECKUPS(i) ((i&(leftRotate(i,1)))||(i&(leftRotate(i,2)))||(i&(leftRotate(i,3))))
#define NOFPS 3
#endif
#ifdef PPPPXPPPP
INTTYPE five=495;
#define CHECKUPS(i) ((i&(leftRotate(i,1)))||(i&(leftRotate(i,2)))||(i&(leftRotate(i,3)))||(i&(leftRotate(i,4))))
#define NOFPS 4
#endif
#ifdef PPPPPXPPPPP
INTTYPE five=2015;
#define CHECKUPS(i) ((i&(leftRotate(i,1)))||(i&(leftRotate(i,2)))||(i&(leftRotate(i,3)))||(i&(leftRotate(i,4)))||(i&(leftRotate(i,5))))
#define NOFPS 5
#endif
bool ComparStruct(const BasNumber &a,const  INTTYPE &b) {
   return(a.basi<b);
}
vector<BasNumber> BasisMat;
vector<int> BasisMatRestrictedLHalf;



void PRINTINFO(const size_t &Ba)
{
    ofstream INFO("INFO"+ap,std::ofstream::out);

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
#ifdef USEPARITY
#ifdef PARITY1
INFO<<"parity 1"<<endl;
#else
INFO<<"parity 2"<<endl;
#endif
#else
INFO<<"Without parity"<<endl;
#endif
INFO<<"L="<<LENGHT<<endl;
INFO<<"basis="<<Ba<<" elements"<<endl;
INFO.close();
}


//const size_t kmomentum=Read("k 0 or L/2");
#ifdef USEPARITY
INTTYPE reflect(const INTTYPE &number)
{
    INTTYPE reflect=0;
    for(size_t j=0;j<LENGHT;j++)
    {
        if(number&(1ULL<<j))reflect=reflect|(1ULL<<(LENGHT-1-j));

    }

    return reflect;
}
#endif
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

INTTYPE leftRotate(const INTTYPE &n, const INTTYPE &d)
{

    return (((n << d)|(n >> (LENGHT - d)))&(refRotate));

}


void checkstate(const INTTYPE &number)
{

    INTTYPE shift=number;
    size_t R;

    for(size_t j=0;j<LENGHT;j++)
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
                //if((kmomentum%L/(j+1))!=0)return;
                    R=j+1;
                    break;
            }
        }

    }
#ifdef USEPARITY
    int m=-1;
    INTTYPE reflected=reflect(number);

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

#ifdef PARITY2
        if(m!=-1)return;
        const BasNumber a = { .basi = number, .Rperio = R};
#else
    const BasNumber a = { .basi = number, .Rperio = R, .Mperio = m };
#endif


#else
const BasNumber a = { .basi = number, .Rperio = R};
#endif

    BasisMat.push_back(a);

    return;
}
void CreateBasis()
{

    INTTYPE first=0;
    checkstate(first);
    BasisMatRestrictedLHalf.push_back(0);
    vector<INTTYPE> VAR;
    VAR.push_back(0);
    for(size_t i=0;i<LENGHT;i++)
    {

        const size_t BAS=VAR.size();


        for(size_t g=0;g<BAS;g++)
        {

            const INTTYPE num=(VAR.at(g))|leftRotate(1,i);


            if(!CHECKUPS(num))
            {
                VAR.push_back(num);

                checkstate(num);
                if(num<pow(2,LENGHT/2))
                {
                    BasisMatRestrictedLHalf.push_back(num);
                }

            }


        }




    }





}

void representative(const INTTYPE &number,INTTYPE &rep
                    #ifdef PARITY2
                    ,double &signo
                    #endif
                    )
{
    rep=number;
    INTTYPE t=number;
    for(size_t i=1;i<LENGHT;i++)
    {
        t=leftRotate(t,1);
        if(t<rep)
        {
            rep=t;
        }

    }

#ifdef USEPARITY
     t=reflect(number);
    for(size_t i=0;i<LENGHT;i++)
    {

        if(t<rep)
        {
            rep=t;
#ifdef PARITY2
            signo=-1.;
#endif

        }
#ifdef PARITY2
        if(t==number)
        {
            rep=-1;
            return;
        }
#endif

        t=leftRotate(t,1);

    }
#endif
}
int main()
{

    ofstream Zvalout("Zval"+ap,std::ofstream::out);
    ofstream outData("outData"+ap,std::ofstream::out);
    ofstream ZvalMultout("ZvalMUlt"+ap,std::ofstream::out);
    ofstream Xvalout("Xval"+ap,std::ofstream::out);
    auto start0 = chrono::high_resolution_clock::now();
    CreateBasis();


    auto start1 = chrono::high_resolution_clock::now();
    cout<<" Elapsed time Creating the Basis " << chrono::duration<double>(start1 - start0).count()<<"s"<<endl;
    cout<<"Creating Hamiltonian matrix"<<endl;

    const size_t BASISSIZE=BasisMat.size();

    const size_t RESTBASISLHALFSIZE=BasisMatRestrictedLHalf.size();
    cout<<"SIZE OF THE BASIS:"<<BASISSIZE<<endl;
    PRINTINFO(BASISSIZE);
    mat Hamil(BASISSIZE,BASISSIZE,fill::zeros);
    for(vector<BasNumber>::iterator  it=BasisMat.begin();it!= BasisMat.end();it++)
    {
        for(size_t pos=0;pos<LENGHT;pos++)
        {
            five=leftRotate(five,1);
            INTTYPE number=it->basi;

            if(!(number&(five)))
            {

                number ^= 1ULL << ((pos+1+NOFPS)%LENGHT);

                INTTYPE rep;

#ifdef PARITY2
                double signo=1.;
                representative(number,rep,signo);
#else
                representative(number,rep);
#endif
#ifdef PARITY2
                if(rep!=-1)
#endif
                {
                    std::vector<BasNumber>::iterator it2 = std::lower_bound(BasisMat.begin(), BasisMat.end(), rep,ComparStruct);
                    if (it2 != BasisMat.end())
                    {
                        const size_t index = std::distance(BasisMat.begin(), it2);
                        const size_t i = std::distance(BasisMat.begin(), it);

#ifdef USEPARITY
#ifdef PARITY1
                        const double val=((it->Mperio==-1)?2.:1.)/((it2->Mperio==-1)?2.:1.);
                        Hamil(i,index)+=sqrt(it->Rperio*val/it2->Rperio)/*((kmomentum==0)?1.:(pow(-1.,l)))*/;
#endif
#ifdef PARITY2
                        Hamil(i,index)+=sqrt(it->Rperio*1./it2->Rperio)/*((kmomentum==0)?1.:(pow(-1.,l)))*/*signo;
#endif

#else
                        Hamil(i,index)+=sqrt(it->Rperio*1./it2->Rperio)/*((kmomentum==0)?1.:(pow(-1.,l)))*/;
#endif

                    }
                    else
                    {
                        cout<<"Error1"<<endl;
                    }

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

#ifdef PARITY2
    outData<< left << setw(14)<<"Energy "<<left << setw(14)  <<"S " <<endl;
#else
    outData<< left << setw(14)<<"Energy "<< left << setw(14) << "NEELOverlap "<<left << setw(14)  <<"S " <<endl;
#endif


      for(size_t j=0;j<eigval.size();j++)
      {
          cout<<"WORKING ON EIGVECTOR "<<j<<endl;
#ifndef PARITY2
          const double overlapNeel=log10(eigvec.at(BASISSIZE-1,j)*eigvec.at(BASISSIZE-1,j));
#endif
          std::array<double,30> Zval={0.};
            std::array<double,30> ZvalMult={0.};
            std::array<double,30> Xval={0.};

          double Sval=0;
          for(size_t k=0;k<BASISSIZE;k++)
          {
              const double vectornorm=pow(eigvec.at(k,j),2);
              INTTYPE numero=BasisMat.at(k).basi;
              const INTTYPE constnum=numero;
              Zval.at(0)+=(2.*__builtin_popcount(numero)-LENGHT)*vectornorm;
              INTTYPE mult=constnum;
              for(size_t m=1;m<=LENGHT/2;m++)
              {
                  numero=leftRotate(numero,1);
                  mult^=(numero);
                  Zval.at(m)+=(LENGHT-2.*__builtin_popcount(constnum^numero))*vectornorm;
                  ZvalMult.at(m-1)+=(LENGHT-2.*__builtin_popcount(mult))*vectornorm;
                  mult^=(refRotate);

                  for(size_t pos=0;pos<LENGHT;pos++)
                  {

                          INTTYPE number=BasisMat.at(k).basi;
                          number ^= 1ULL << ((pos)%LENGHT);
                          number ^= 1ULL << ((pos+m)%LENGHT);
                          if(!CHECKUPS(number))
                          {
                              INTTYPE rep;

              #ifdef PARITY2
                              double signo=1.;
                              representative(number,rep,signo);
              #else
                              representative(number,rep);
              #endif
              #ifdef PARITY2
                              if(rep!=-1)
              #endif
                              {
                                  std::vector<BasNumber>::iterator it2 = std::lower_bound(BasisMat.begin(), BasisMat.end(), rep,ComparStruct);
                                  if (it2 != BasisMat.end())
                                  {
                                      const size_t index = std::distance(BasisMat.begin(), it2);


              #ifdef USEPARITY
              #ifdef PARITY1
                                      const double val=((BasisMat.at(k).Mperio==-1)?2.:1.)/((it2->Mperio==-1)?2.:1.);
                                      Xval.at(m-1)+=sqrt(BasisMat.at(k).Rperio*val/it2->Rperio)*eigvec.at(k,j)*eigvec.at(index,j)/*((kmomentum==0)?1.:(pow(-1.,l)))*/;
              #endif
              #ifdef PARITY2
                                      Xval.at(m-1)+=sqrt(BasisMat.at(k).Rperio*1./it2->Rperio)*eigvec.at(k,j)*eigvec.at(index,j)/*((kmomentum==0)?1.:(pow(-1.,l)))*/*signo;
              #endif

              #else
                                      Xval.at(m-1)+=sqrt(BasisMat.at(k).Rperio*1./it2->Rperio)*eigvec.at(k,j)*eigvec.at(index,j)/*((kmomentum==0)?1.:(pow(-1.,l)))*/;
              #endif

                                  }

                          }

                    }


        }


              }


          }




          mat fullBasis(RESTBASISLHALFSIZE,RESTBASISLHALFSIZE,fill::zeros);
          size_t irow=0;

          for(vector<int>::iterator  itR=BasisMatRestrictedLHalf.begin();itR!= BasisMatRestrictedLHalf.end();itR++)
          {
              size_t icolumn=0;
              for(vector<int>::iterator  itL=BasisMatRestrictedLHalf.begin();itL!= BasisMatRestrictedLHalf.end();itL++)
              {

                    const INTTYPE num=((*itR)|leftRotate(*itL,LENGHT/2));

                    if(!CHECKUPS(num))
                    {
                        INTTYPE theRep;

#ifdef PARITY2
                        double signo=1.;
                        representative(num,theRep,signo);
#else
                        representative(num,theRep);
#endif
#ifdef PARITY2
                        if(theRep!=-1)
#endif
                        {
                            std::vector<BasNumber>::iterator it = std::lower_bound(BasisMat.begin(), BasisMat.end(), theRep,ComparStruct);

                            const size_t index = std::distance(BasisMat.begin(), it);
                            #ifdef USEPARITY
#ifdef PARITY1
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)/*((kmomentum==0)?1.:(pow(-1.,thel)))*/*
                                                       sqrt(((it->Mperio==-1)?0.5:1.)/it->Rperio);
#endif
#ifdef PARITY2
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)/*((kmomentum==0)?1.:(pow(-1.,thel)))*/*
                                                       sqrt(0.5/it->Rperio)*signo;
#endif

                            #else
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)/*((kmomentum==0)?1.:(pow(-1.,thel)))*/*
                                                       sqrt(1./it->Rperio);
                            #endif

                        }


                    }


                    icolumn++;

              }

              irow++;

          }


           mat Reduced=fullBasis*fullBasis.t();


           vec eigval2;
           mat eigvec2;


           eig_sym(eigval2, eigvec2, Reduced);


           Sval=trace(eigval2.t()*trunc_log(eigval2));
#ifndef PARITY2
            outData<< left << setw(14)<<eigval(j)<< left << setw(14)<<overlapNeel<< left << setw(14) <<-Sval<<endl;
#else
                    outData<< left << setw(14)<<eigval(j)<< left << setw(14) <<-Sval<<endl;
#endif
    Zvalout<< left << setw(14);
    ZvalMultout<< left << setw(14);
    Xvalout<< left << setw(14);

    for(size_t p=0;p<LENGHT/2;p++)
    {
        Zvalout<<Zval.at(p)<< left << setw(14);
        ZvalMultout<<ZvalMult.at(p)<< left << setw(14);
        Xvalout<<Xval.at(p)<< left << setw(14);
    }
    Zvalout<<endl;
    ZvalMultout<<endl;
    Xvalout<<endl;


      }

      auto start4 = chrono::high_resolution_clock::now();

      cout<<" Elapsed time working on eigenvectors: " << chrono::duration<double>(start4 - start3).count()<<"s"<<endl;

    outData.close();
    Zvalout.close();
    ZvalMultout.close();
    Xvalout.close();



}




