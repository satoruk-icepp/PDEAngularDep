void PlotRI(){
   ifstream fin;
   fin.open("SiliconRefractiveIndices.dat",std::ios::in);
   Double_t wavelength;
   Double_t n;
   Double_t k;
   TCanvas* cri = new TCanvas("cri","cri");
   TGraph* grn=new TGraph();
   TGraph* grk=new TGraph();
   while (!fin.eof()) {
      fin>>wavelength>>n>>k;
      grn->SetPoint(grn->GetN(),wavelength*TMath::Power(10,6),n);
      grk->SetPoint(grk->GetN(),wavelength*TMath::Power(10,6),k);
   }
   grn->SetTitle("Silicon Refractive Indices;#lambda[nm];n or k");
   gStyle->SetTitleFontSize(.08);
   gStyle->SetTitleOffset(.6, "XY");
   gStyle->SetLabelSize(.05, "XY");
   gPad->SetBottomMargin(0.15);
   gPad->SetLeftMargin(0.15);
   grn->GetXaxis()->SetTitleSize(0.1);
   grn->GetYaxis()->SetTitleSize(0.1);
   grn->SetMaximum(6);
   grn->SetMinimum(0);
   grn->SetMarkerColor(kRed);
   grn->SetMarkerStyle(20);
   grn->Draw("ap");
   std::cout<<"n at"
   <<" 170nm: "<<grn->Eval(170)
   <<" 175nm: "<<grn->Eval(175)
   <<" 180nm: "<<grn->Eval(180)<<std::endl;
   std::cout<<"k at"
   <<" 170nm: "<<grk->Eval(170)
   <<" 175nm: "<<grk->Eval(175)
   <<" 180nm: "<<grk->Eval(180)<<std::endl;
   grk->SetMarkerColor(kBlue);
   grk->SetMarkerStyle(20);
   grk->Draw("p same");
   auto legend = new TLegend(0.15,0.7,0.3,0.9);
   legend->AddEntry(grn,"n","lep");
   legend->AddEntry(grk,"k","lep");
   legend->Draw();
   // cri->Update();
}