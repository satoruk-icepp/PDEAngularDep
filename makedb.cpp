void makedb(){
   ifstream infile;
   ofstream outfile;
   infile.open("MPPCPDE.dat",std::ios::in);
   outfile.open("XECQE.dat",std::ios::out);
   std::map<Int_t,Double_t> pdemap;
   std::map<Int_t,Double_t> pdeerrmap;
   Int_t ch;
   Double_t PDE;
   Double_t PDEerr;
   Int_t DBid=1;
TH2D* hPDE = new TH2D("hpde","hpde",44,0,44,93,0,93);
TGraphErrors* grpde = new TGraphErrors();
   while (!infile.eof()) {
      infile>>ch>>PDE>>PDEerr;
      pdemap[ch]=PDE;
      pdeerrmap[ch]=PDEerr;
      std::cout<<"ch: "<<ch<<" PDE: "<<PDE<<std::endl;
if (ch<4092) {
   hPDE->Fill(ch%44,ch/44,PDE);
}

   }
   
   for (Int_t ipm = 0; ipm < 4760; ipm++) {
      if (pdemap.count(ipm)<1) {
         outfile<<DBid<<","<<ipm<<","<<0<<","<<0<<std::endl;
      }else{
         outfile<<DBid<<","<<ipm<<","<<pdemap[ipm]<<","<<pdeerrmap[ipm]<<std::endl;
      }
   }
hPDE->Draw("colz");
}