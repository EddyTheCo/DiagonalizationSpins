#include "myTGraph.hh"
#include<TGraphErrors.h>
#include <TFrame.h>
#include <TCanvas.h>
#include<TAxis.h>
#include<TObjArray.h>
#include<myTGraph.hh>
#include<TF1.h>


using namespace std;



MyTGraph::MyTGraph(Int_t n, const Double_t *x, const Double_t *y, const Double_t *ex, const Double_t *ey, const char *title, double T):TGraphErrors(n,x,y, ex, ey)
{
    TCanvas *c1 = new TCanvas("c1","Root",200,100,900,600);
    //c1->SetLogy();
    //c1->SetLogx();
   // if(T==2.269)
    //{
      //  TF1 *f = new TF1("f", "x^[0]");
        //this->Fit(f);
    //}
    //this->SetLineColor();
    //this->SetLineWidth(m+1);
    this->SetMarkerColor(223);
    this->SetMarkerStyle(22);
    this->GetXaxis()->SetTitle("l");
    this->GetYaxis()->SetTitle("tao");
    this->SetTitle(title);
    this->Draw("ap");
     // TCanvas::Update() draws the frame, after which one can change it
     c1->Update();
     c1->GetFrame()->SetFillColor(21);
     c1->GetFrame()->SetBorderSize(12);
     c1->Modified();
     c1->SaveAs("picture.jpg");



}
MyTGraph::~MyTGraph()
{
	
	
}
