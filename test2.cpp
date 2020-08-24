Double_t PoisGaus(Double_t* x, Double_t* par);

void test2(){
    TF1* fPois= new TF1("fPois","[0]*TMath::Poisson(x,[1])",0,100);
  TF1* fGaus= new TF1("fGaus","TMath::Gaus(x,0,[0])",-20,20);
  TF1Convolution *f_conv = new TF1Convolution(fPois,fGaus,-5,100,true);
   f_conv->SetRange(-5.,100.);
   f_conv->SetNofPointsFFT(1000);
   TF1   *f = new TF1("f",*f_conv, -5, 100, f_conv->GetNpar());
   f->SetParameters(10,5,0.5);
  TF1* ftest = new TF1("ftest",PoisGaus,-5,100,5);
  ftest->SetParameters(100,5,3,-5,45);
  f->Draw();
}

Double_t PoisGaus(Double_t* x, Double_t* par) {
   Double_t val = 0.;
   int nbin = 100;
   Double_t lowl = par[3];
   Double_t higl = par[4];
   for (int i = 0; i < nbin; i++) {
	 double y = lowl+(higl-lowl)*i/(double)(nbin+1.);
      val += par[0]*TMath::Poisson(y,par[1]) * TMath::Gaus(x[0] - y,0, par[2]) / nbin;
   }
   // cout<<x[0]<<" "<<val<<endl;
   return val;
}
// Double_t PoisGaus(Double_t* x,Double_t* par){
//   TF1* fPois= new TF1("fPois","[0]*TMath::Poisson(x,[1])",0,100);
//   TF1* fGaus= new TF1("fGaus","TMath::Gaus(x,0,[0])",-20,20);
//   TF1Convolution *f_conv = new TF1Convolution(fPois,fGaus,-5,100,true);
//    f_conv->SetRange(-5.,100.);
//    f_conv->SetNofPointsFFT(1000);
//    TF1   *f = new TF1("f",*f_conv, 0., 5., f_conv->GetNpar());
//    f->SetParameters(par[0],par[1],par[2]);
//   return 
// }
