#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define ARMA_DONT_USE_WRAPPER
#define ARMA_NO_DEBUG
#include <armadillo>




using namespace std;
using namespace arma;

void Thouless(void)
{
    mat eigV;
    eigV.load("eigvalues.dat");

    vec x = conv_to< vec >::from(eigV);
    vec y= vec(x.n_elem);



    for(size_t i=0;i<x.n_elem;i++)
    {
        y.at(i)=i+1;
    }

    vec p = polyfit(x,y,10);

    vec y2 = polyval(p,x);


    vec GaussFilte=exp(-2*square(y2-mean(y2))/var(y2));


    vec Zvec=square(GaussFilte);


    double Z=sum(Zvec);


    ofstream THOUX("THOUX");
    ofstream THOUY("THOUY");
    ofstream THOUZ("THOUZ");

    for (double tau=0; tau<1.2;tau+=0.001)
    {
        cx_vec fase(x.n_elem);

        for (size_t i = 0; i < x.n_elem; i++) {
            fase.at(i)=polar(1.0,-y2.at(i)*tau);

        }


        double suma=norm(dot(GaussFilte,fase));
        THOUX<<tau<<endl;
        THOUY<<tau <<" "<<suma<<endl;
    }
    THOUZ<<Z<<endl;
    THOUX.close();
    THOUY.close();
    THOUZ.close();




}

int main()
{
    Thouless();
}
