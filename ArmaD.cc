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

#define PPPXPPP
#define LENGHT 28
#define INTTYPE int
#define USEPARITY
#define PBC


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
#ifdef PBC
    size_t Rperio;
#ifdef PARITY1
    int Mperio;
#endif
#endif
};
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
INFO<<"basis="<<Ba<<endl;
#ifdef OBC
INFO<<"OBC"<<endl;
#else
INFO<<"PBC"<<endl;
#endif
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

void checkstate(const INTTYPE &number)
{
#ifdef PBC
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
#endif

#ifdef USEPARITY

    INTTYPE reflected=reflect(number);
#ifdef PBC
    int m=-1;
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
    if(reflected<number)
    {
        return;
    }
#ifdef PARITY2
    if(reflected==number)return;
#endif
    const BasNumber a = { .basi = number};
#endif

#else
#ifdef PBC
const BasNumber a = { .basi = number, .Rperio = R};
#else
    const BasNumber a = { .basi = number};
#endif
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
#ifdef PBC

    for(size_t i=1;i<LENGHT;i++)
    {
        t=leftRotate(t,1);
        if(t<rep)
        {
            rep=t;
        }

    }
#endif
#ifdef USEPARITY
    t=reflect(number);
#ifdef PBC
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
#else
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

#endif
#endif
}
int main()
{

    ofstream Zvalout("Zval"+ap,std::ofstream::out);
    ofstream Sparcityout("Sparcity"+ap,std::ofstream::out);
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

            INTTYPE fiveVar=LEFTROT(five,pos-NOFPS);

            INTTYPE number=it->basi;

            if(!(number&(fiveVar)))
            {


                number ^= 1ULL << pos;

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
#ifdef PBC
                        const double val=((it->Mperio==-1)?2.:1.)/((it2->Mperio==-1)?2.:1.);
                        Hamil(i,index)+=sqrt(it->Rperio*val/it2->Rperio)/*((kmomentum==0)?1.:(pow(-1.,l)))*/;
#else
                        Hamil(i,index)+=sqrt(((it->basi!=reflect(it->basi))?2.:1.)/((rep!=reflect(rep))?2.:1.));
#endif

#endif
#ifdef PARITY2
#ifdef PBC
                        Hamil(i,index)+=sqrt(it->Rperio*1./it2->Rperio)/*((kmomentum==0)?1.:(pow(-1.,l)))*/*signo;
#else
                        Hamil(i,index)+=signo;
#endif
#endif

#else
#ifdef PBC
                        Hamil(i,index)+=sqrt(it->Rperio*1./it2->Rperio)/*((kmomentum==0)?1.:(pow(-1.,l)))*/;
#else
                        Hamil(i,index)+=1;
#endif
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
    eigvec.save("eigvectors.dat");

    for(size_t j=0;j<eigval.size();j++)
    {
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
#ifdef PBC
                      fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                 sqrt(((it->Mperio==-1)?0.5:1.)/it->Rperio);
#else
                      fullBasis.at(irow,icolumn)=eigvec.at(index,j)*(sqrt((theRep!=reflect(theRep)?0.5:1.)));
#endif
#endif
#ifdef PARITY2
#ifdef PBC
                      fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                 sqrt(0.5/it->Rperio)*signo;
#else
                      fullBasis.at(irow,icolumn)=eigvec.at(index,j)*signo*sqrt(0.5);
#endif
#endif

                      #else
#ifdef PBC
                      fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                 sqrt(1./it->Rperio);
#else
                      fullBasis.at(irow,icolumn)=eigvec.at(index,j);
#endif
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
     double summa=0;
     for(size_t h=0;h<eigval2.size();h++ )
     {
        if(eigval2.at(h)>0.000000000001)
        {
           summa++;
        }

     }

     Sparcityout<<eigval.at(j)<<" "<<summa/eigval2.size()<<endl;

     if(summa/eigval2.size()<0.2)
     {
         vec va=eigvec.col(j);
         va.save("eigvector"+to_string(j)+".dat");

     }


}

Sparcityout.close();

/*
    outData<< left << setw(14)<<"Energy "<<left << setw(14)  <<"S " <<endl;
    Zvalout<<" Z1   Z1Z2"<<endl;
    ZvalMultout<<" Z1Z2   Z1Z23"<<endl;
    Xvalout<<" X1X2   X1X3"<<endl;

      for(size_t j=0;j<eigval.size();j++)
      {
          cout<<"WORKING ON EIGVECTOR "<<j<<endl;
          std::array<double,LENGHT/2+1> Zval={0.};
            std::array<double,LENGHT/2+1> ZvalMult={0.};
            std::array<double,LENGHT/2+1> Xval={0.};

          double Sval=0;
          for(size_t k=0;k<BASISSIZE;k++)
          {
              const double vectornorm=pow(eigvec.at(k,j),2);
              INTTYPE numero=BasisMat.at(k).basi;
              const INTTYPE constnum=numero;
              Zval.at(0)+=(2.*__builtin_popcount(numero)-LENGHT)*vectornorm;
              INTTYPE mult=constnum;
              for(size_t m=1;m<LENGHT/2;m++)
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
#ifdef PBC
                                      const double val=((BasisMat.at(k).Mperio==-1)?2.:1.)/((it2->Mperio==-1)?2.:1.);
                                      Xval.at(m-1)+=sqrt(BasisMat.at(k).Rperio*val/it2->Rperio)*eigvec.at(k,j)*eigvec.at(index,j);
#else
                                      Xval.at(m-1)+=sqrt(((BasisMat.at(k).basi!=reflect(BasisMat.at(k).basi))?2.:1.)/((rep!=reflect(rep))?2.:1.))*eigvec.at(k,j)*eigvec.at(index,j);
#endif
              #endif
              #ifdef PARITY2
#ifdef PBC
                                      Xval.at(m-1)+=sqrt(BasisMat.at(k).Rperio*1./it2->Rperio)*eigvec.at(k,j)*eigvec.at(index,j)*signo;
#else
                                      Xval.at(m-1)+=eigvec.at(k,j)*eigvec.at(index,j)*signo;
#endif
              #endif

              #else
#ifdef PBC
                                      Xval.at(m-1)+=sqrt(BasisMat.at(k).Rperio*1./it2->Rperio)*eigvec.at(k,j)*eigvec.at(index,j);
#else
                                      Xval.at(m-1)+=eigvec.at(k,j)*eigvec.at(index,j);
#endif
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
#ifdef PBC
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                       sqrt(((it->Mperio==-1)?0.5:1.)/it->Rperio);
#else
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)*(sqrt((theRep!=reflect(theRep)?0.5:1.)));
#endif
#endif
#ifdef PARITY2
#ifdef PBC
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                       sqrt(0.5/it->Rperio)*signo;
#else
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)*signo*sqrt(0.5);
#endif
#endif

                            #else
#ifdef PBC
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j)*
                                                       sqrt(1./it->Rperio);
#else
                            fullBasis.at(irow,icolumn)=eigvec.at(index,j);
#endif
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
           outData<< left << setw(14)<<eigval(j)<< left << setw(14) <<-Sval<<endl;

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

*/

}




