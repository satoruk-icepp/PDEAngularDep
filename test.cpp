#include "Fresnel.h"

Double_t wavelength = 175*nm;
Double_t dLayer=60*nm;
void test(){
   TCanvas* cr = new TCanvas("cr","cr",600,600);
   // TF1* frp=new TF1("frp",FresnelReflection_P,0,90*degree,8);
   // TF1* frs=new TF1("frs",FresnelReflection_S,0,90*degree,8);
   std::complex<Double_t> nLXe={1.64,0};
   std::complex<Double_t> nLayer={1.61,0};
   std::complex<Double_t> nSi={0.692533,2.4685};
   TGraphErrors* gr= new TGraphErrors();
   Double_t par[]={real(nLXe),imag(nLXe),real(nLayer),imag(nLayer),real(nSi),imag(nSi),wavelength,dLayer};
   for (Int_t iang = 0; iang < 180; iang++) {
      // Double_t refangle=RefractionAngle(nLXe,nLayer,iang*degree);
      Double_t x[]={iang*degree/2.};
      
      // gr->SetPoint(iang,iang,norm(FresnelCoefficientR_P(nLayer,nSi,refangle)));
      gr->SetPoint(iang,iang,FresnelReflection_S(x,par));
   }
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(2);
   gr->Draw("ap");
   
}

