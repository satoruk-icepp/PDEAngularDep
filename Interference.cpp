#include "Fresnel.h"


// std::complex<Double_t> nLXe={1,0};
// std::complex<Double_t> nLayer={1.2,0};
// std::complex<Double_t> nSi={1.44,0};
// Double_t wavelength = 633*nm;
// Double_t wavelength = 175*nm;
// Double_t dLayer=500*microm;
// Double_t dLayer=500*nm;
// Double_t dLayer=1*microm;
Double_t dLayer=60*nm;
void Interference(){
   gStyle->SetOptStat(0);
   // TFile* fin = new TFile("./Projects/MPPCPDE.root","read");
   // TFile* flp = new TFile("./Projects/MPPCAngDep/pde_angle_attcorr.root","read");
   // TH2D* htemp = (TH2D*)fin->Get("hangledep");
   // TGraph* grtemp = (TGraph*)flp->Get("Graph");
   // TGraph* grnorm = new TGraph();
   // for (Int_t ipoint = 0; ipoint < grtemp->GetN(); ipoint++) {
   //    grnorm->SetPoint(ipoint,grtemp->GetX()[ipoint],grtemp->GetY()[ipoint]/0.18);
   // }
   
   Int_t Nconf =3;
   TF1* frt;
   TF1* frp;
   TF1* frs;
   TF1* fps;
   TF1* fcomp;
   // std::vector<Double_t> dLayer={60*nm,80*nm,200*nm};
   // for (Int_t iconf = 0; iconf < 1; iconf++) {
   frt=new TF1("frt",FresnelReflection_T,0,90*degree,8);
   frp=new TF1("frp",FresnelReflection_P,0,90*degree,8);
   frs=new TF1("frs",FresnelReflection_S,0,90*degree,8);
   fps=new TF1("fps",FresnelPhaseShift_P,0,90*degree,8);
   // }
   std::vector<std::complex<Double_t>> nLayers={
      // {1.64,0},
      {1.,0},
      // {1.2,0},
      // {1.4,0}
      {1.61,0},
      // {1.,0}
      // {1.61,0}
      {0.692533,2.4685}
      // {1.64,0}
   };
   // Double_t WLCenter = 633*nm;
   Double_t WLCenter = 175.9*nm;
   Double_t WLSigma = 4.3*nm;
   TF1* fwl = new TF1("fwl","TMath::Gaus(x,[0],[1])",WLCenter-3*WLSigma,WLCenter+3*WLSigma);
   // fwl->SetParameters(175.9,10.2);
   
   fwl->SetParameters(WLCenter,WLSigma);
   Int_t Nwls = 1000;
   TH1D* hrt=new TH1D("hrt","hrt",90,0,90);
   TH1D* hrp=new TH1D("hrp","hrp",90,0,90);
   TH1D* hrs=new TH1D("hrs","hrs",90,0,90);
   TH1D* hwl = new TH1D("hwl","hwl",100,(WLCenter-3*WLSigma)/nm,(WLCenter+3*WLSigma)/nm);
   TH1D* hps = new TH1D("hps","hps",90,0,90);
   // TGraphErrors* gps = new TGraphErrors();
   for (Int_t iwl = 0; iwl < Nwls; iwl++) {
      Double_t wavelength = fwl->GetRandom();
      hwl->Fill(wavelength/nm);
      frt->SetParameters(real(nLayers[0]),imag(nLayers[0]),real(nLayers[1]),imag(nLayers[1]),real(nLayers[2]),imag(nLayers[2]),wavelength,dLayer);
      frp->SetParameters(real(nLayers[0]),imag(nLayers[0]),real(nLayers[1]),imag(nLayers[1]),real(nLayers[2]),imag(nLayers[2]),wavelength,dLayer);
      frs->SetParameters(real(nLayers[0]),imag(nLayers[0]),real(nLayers[1]),imag(nLayers[1]),real(nLayers[2]),imag(nLayers[2]),wavelength,dLayer);
      fps->SetParameters(real(nLayers[0]),imag(nLayers[0]),real(nLayers[1]),imag(nLayers[1]),real(nLayers[2]),imag(nLayers[2]),wavelength,dLayer);
      for (Int_t iangle = 0; iangle < 90; iangle++) {
         hrt->Fill(iangle,frt->Eval(iangle*TMath::DegToRad())/(Double_t)Nwls);
         hrp->Fill(iangle,frp->Eval(iangle*TMath::DegToRad())/(Double_t)Nwls);
         hrs->Fill(iangle,frs->Eval(iangle*TMath::DegToRad())/(Double_t)Nwls);
         hps->Fill(iangle,fps->Eval(iangle*TMath::DegToRad())/(Double_t)Nwls);
      }
   }
   // TCanvas* cweighted = new TCanvas("cwe","cwe",1200,600);
   TCanvas* cweighted = new TCanvas("cwe","cwe",600,600);
   // cweighted->Divide(2,1);
   // cweighted->cd(1);
   hrt->SetMarkerColor(1);
   hrt->SetMarkerStyle(20);
   hrt->SetMaximum(1);
   hrt->SetMinimum(0);
   
   hrt->Draw("hist p");
   hrp->SetMarkerColor(2);
   hrp->SetMarkerStyle(20);
   hrp->Draw("same hist p");
   hrs->SetMarkerColor(3);
   hrs->SetMarkerStyle(20);
   hrs->Draw("same hist p");
   hrt->SetTitle("R_{total}");
   hrs->SetTitle("R_{s}");
   hrp->SetTitle("R_{p}");
   // gPad->BuildLegend(0.1,0.6,0.4,0.9);//left top
   gPad->BuildLegend(0.6,0.13,0.9,0.4);//left top
   gPad->SetBottomMargin(0.13);
   hrt->SetTitle("R vs Incident angle #theta;#theta[deg];R");
   hrt->GetXaxis()->SetTitleSize(0.08);
   hrt->GetXaxis()->SetTitleOffset(0.7);
   hrt->GetYaxis()->SetTitleSize(0.08);
   hrt->GetYaxis()->SetTitleOffset(0.5);
   TCanvas* cot = new TCanvas("cot","cot");
   hwl->SetTitle("Generated wavelength #lambda;#lambda[nm];# of #lambda");
   gPad->SetBottomMargin(0.13);
   hwl->GetXaxis()->SetTitleSize(0.08);
   hwl->GetXaxis()->SetTitleOffset(0.7);
   hwl->GetYaxis()->SetTitleSize(0.08);
   hwl->GetYaxis()->SetTitleOffset(0.5);
gStyle->SetTitleAlign(33);
gStyle->SetTitleX(0.99);
   hwl->Draw();
   // cweighted->cd(2);
   // // fwl->Draw();
   // // hwl->Draw();
   // hps->Draw();
   // TCanvas* cr = new TCanvas("cr","cr",600,600);
   // TF1* frp=new TF1("frp",FresnelReflection_P,0,90*degree,8);
   // TF1* frs=new TF1("frs",FresnelReflection_S,0,90*degree,8);
   // std::complex<Double_t> nQuartz={1.64,0};
   // std::complex<Double_t> nLXe={1.64,0};
   // std::complex<Double_t> nLayer={1.61,0};
   // std::complex<Double_t> nSi={0.692533,2.4685};
   
   
   
   // fcomp=new TF1("trans",Trans,0,90,8);      
   // fcomp->SetParameters(real(nLayers[0]),imag(nLayers[0]),real(nLayers[1]),imag(nLayers[1]),real(nLayers[2]),imag(nLayers[2]),wavelength,dLayer);
   // frt[iconf]->SetName("total reflection");
   // frt->SetLineColor(2+iconf);
   // frt->SetLineWidth(1);
   // frt->SetNpx(100000);
   // frt->SetMaximum(1);
   // frt->SetMinimum(0);
   // frs->SetLineColor(3+iconf);
   // frp->SetLineColor(4+iconf);
   // frt->GetXaxis()->SetNdivisions(-503);
   // frt->GetXaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"0");
   // frt->GetXaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"30");
   // frt->GetXaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"60");
   // frt->GetXaxis()->ChangeLabel(-1,-1,-1,-1,-1,-1,"90");
   // // if (iconf==0) {
   // frt->Draw();
   // // }else{
   // // frt[iconf]->Draw("same");
   // //    frp[iconf]->Draw("same");
   // // frs[iconf]->Draw("same");
   // frt[iconf]->SetTitle("R_{total}");
   // frs[iconf]->SetTitle("R_{s}");
   // frp[iconf]->SetTitle("R_{p}");
   // gPad->BuildLegend(0.1,0.8,0.4,0.9);
   // // frt[iconf]->SetTitle("R_{total}");
   // frt[iconf]->SetTitle("Reflectivity@gXe;Incident angle[deg];R");
   // }
   
   // std::cout<<"at 0: "<<1-frt[0]->Eval(0)<<std::endl;
   
   // TCanvas* cps = new TCanvas("cps","cps");
   // cps->cd();
   // TF1* fPS= new TF1("phaseshift",PhaseShift,0,90*degree,8);
   // fPS->SetParameters(real(nLXe),imag(nLXe),real(nLayer),imag(nLayer),real(nSi),imag(nSi),wavelength,dLayer);
   // fPS->Draw();
   // 
   // gStyle->SetOptStat(0);
   // TCanvas* ccomp=new TCanvas("ccomp","ccomp");
   // htemp->SetTitle("Angular Dependence of PDE;Incident Angle[deg];Normalized PDE");
   // htemp->GetXaxis()->SetTitleSize(0.08);
   // htemp->GetXaxis()->SetTitleOffset(0.5);
   // htemp->GetYaxis()->SetTitleSize(0.08);
   // htemp->GetYaxis()->SetTitleOffset(0.5);
   // htemp->Draw("colz");
   // grnorm->SetMarkerColor(2);
   // grnorm->SetMarkerStyle(4);
   // grnorm->Draw("p same");
   // fcomp->Draw("same");
   // auto legend = new TLegend(0.1,0.7,0.4,0.9);
   // legend->AddEntry(htemp,"Measured in actual machine","f");
   // legend->AddEntry(grnorm,"Measured in large prototype","p");
   // legend->AddEntry(fcomp,"Expected","l");
   // legend->Draw();
   // frp->SetParameters(real(nLXe),imag(nLXe),real(nLayer),imag(nLayer),real(nSi),imag(nSi),wavelength,dLayer);
   // frs->SetParameters(real(nLXe),imag(nLXe),real(nLayer),imag(nLayer),real(nSi),imag(nSi),wavelength,dLayer);
   
   // frp->SetTitle("Reflectivity;Incident angle[rad];R");
   
   // frp->SetMaximum(1);
   // frp->SetMinimum(0);
   // frp->SetLineColor(kRed);
   // frp->SetNpx(10000);
   // frp->SetLineWidth(5);
   // frp->Draw();
   // 
   // frs->SetLineColor(kBlue);
   // frs->SetLineWidth(5);
   // frs->SetNpx(10000);
   // frs->Draw("same");
   
}
