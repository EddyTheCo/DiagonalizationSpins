#include<TGraphErrors.h>
#include <TFrame.h>
#include <TCanvas.h>
#include<TAxis.h>
#include<TObjArray.h>





using namespace std;

#ifndef MyTGraph_h
#define MyTGraph_h 1


class MyTGraph:public TGraphErrors
{
	public:
    MyTGraph(Int_t n, const Double_t *x, const Double_t *y, const Double_t *ex=0, const Double_t *ey=0, const char *title="",double T=0);
~MyTGraph();




};

#endif
