Double_t degree    = TMath::Pi()/180.;
Double_t radtodeg  = 1./degree;
Double_t MPPCSize  = 12;
std::complex<Double_t> nxe(1.64,0);//LXe,quartz
std::complex<Double_t> nsi(0.7,2.5);//Silicon

Double_t CalculateReflectivity_P (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t theta);
Double_t CalculateReflectivity_S (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t theta);
Double_t PMIncidentAngle(TVector3& view, TVector3& center,TVector3& normal);
Double_t TransmittanceOfSilicon (std::complex<Double_t> n1, std::complex<Double_t> n2,Double_t theta);
Double_t PMSolidAngleMPPC(TVector3& view, TVector3& center,TVector3& normal);
Double_t PMAngleDependenceOfPDE(Double_t theta);
Bool_t JudgeDetection(TVector3& Direction, TVector3& MPPCCenter);

Double_t RS(Double_t *x,Double_t *par){
return CalculateReflectivity_S(nxe,nsi,x[0]);
}

Double_t RP(Double_t *x,Double_t *par){
return CalculateReflectivity_P(nxe,nsi,x[0]);
}

void AngConv(){
TCanvas* cr  = new TCanvas("cr","cr",1000,600);
cr->Divide(2,1);
TF1* frs = new TF1("frs",RS,0,90,0);
TF1* frp = new TF1("frp",RP,0,90,0);
cr->cd(1);
frs->Draw();
cr->cd(2);
frp->Draw();

	// gRandom->SetSeed(0);
	// TVector3 Origin(0,0,0);
	// TCanvas* canvas1= new TCanvas("canvas1","canvas1",600,600);
	// TCanvas* canvas2= new TCanvas("canvas2","canvas2",600,600);
	// Int_t Npho=500000;
   // 
   // 
	// Double_t MPPCtheta = 0* degree;
	// Double_t MPPCphi   = 0* degree;
	// Double_t MPPCR     = 20;
	// TH2D* hXY= new TH2D("hXY","hXY",100,-100,100,100,-100,100);
	// TH1D* htheta= new TH1D("htheta","htheta",100,0,TMath::Pi());
	// TGraph* grSolidAngle= new TGraph();
	// TGraph* grMCSolid= new TGraph();
	// TGraph* grADPDE= new TGraph();
	// TGraph* grNphe = new TGraph();
	// TGraph* grTS= new TGraph();
	// for (int iangle = 0; iangle < 30; iangle++) {
	// 	Double_t count  = 0;
	// 	Double_t ndirect = 0;
	// 	MPPCtheta=iangle * 3 * degree;
	// 	// grADPDE->SetPoint(grADPDE->GetN(),MPPCtheta*radtodeg,PMAngleDependenceOfPDE(MPPCtheta*radtodeg));
	// 	// TVector3 Direction;
	// 	TVector3 MPPCNorm(0,0,1);
	// 	TVector3 MPPCCenter(
	// 		MPPCR*TMath::Sin(MPPCtheta) * TMath::Cos(MPPCphi),
	// 		MPPCR*TMath::Sin(MPPCtheta) * TMath::Sin(MPPCphi),
	// 		MPPCR*TMath::Cos(MPPCtheta)
	// 	);
	// 	Double_t MPPCXMin =MPPCCenter.X()-MPPCSize/2.;
	// 	Double_t MPPCXMax =MPPCCenter.X()+MPPCSize/2.;
	// 	Double_t MPPCYMin =MPPCCenter.Y()-MPPCSize/2.;
	// 	Double_t MPPCYMax =MPPCCenter.Y()+MPPCSize/2.;
	// 	Double_t MPPCZ    =MPPCCenter.Z();
	// 	std::cout<<"MPPCXmin:"<<MPPCXMin<<std::endl;
	// 	std::cout<<"MPPCXmax:"<<MPPCXMax<<std::endl;
	// 	std::cout<<"MPPCYmin:"<<MPPCYMin<<std::endl;
	// 	std::cout<<"MPPCYmax:"<<MPPCYMax<<std::endl;
	// 	std::cout<<"MPPCZ:   "<<MPPCZ<<std::endl;
	// 	// TVector3 MPPCDirection(
	// 	//
	// 	// )
	// 	for (int ipho = 0; ipho < Npho; ipho++) {
	// 		Double_t costheta = gRandom->Uniform(0,1);
	// 		Double_t phi      = gRandom->Uniform(-TMath::Pi(),TMath::Pi());
	// 		// Double_t costheta = 0;
	// 		// Double_t phi      = 0;
	// 		htheta->Fill(phi);
	// 		TVector3 Direction(
	// 			TMath::Sqrt(1-costheta*costheta) * TMath::Cos(phi),
	// 			TMath::Sqrt(1-costheta*costheta) * TMath::Sin(phi),
	// 			costheta
	// 		);
	// 		// Double_t Scale=0;
	// 		// if(Direction.Z()<0) continue;
	// 		// Scale= MPPCCenter.Z()/Direction.Z();
	// 		//
	// 		//
	// 		// // std::cout<<"scale: "<<Scale<<std::endl;
	// 		// // std::cout<<"Direction "<<Direction.X()<<" "<<Direction.Y()<<" "<<Direction.Z()<<std::endl;
	// 		// Double_t XatMPPC  = Direction.X()*Scale;
	// 		// Double_t YatMPPC  = Direction.Y()*Scale;
	// 		// if (iangle==0) {
	// 		// 	hXY->Fill(XatMPPC,YatMPPC);
	// 		// }
	// 		if (JudgeDetection(Direction,MPPCCenter)) {
	// 			Double_t incangle=PMIncidentAngle(Origin,MPPCCenter,MPPCNorm);
	// 			count=count+TransmittanceOfSilicon(incangle*radtodeg)/TransmittanceOfSilicon(0);
	// 			ndirect=ndirect+1;
	// 		}
   // 
	// 		// if (XatMPPC>MPPCXMin&&XatMPPC<MPPCXMax) {
	// 		// 	if (YatMPPC>MPPCYMin&&YatMPPC<MPPCYMax){
	// 		// 		count=count+PMAngleDependenceOfPDE(incangle*radtodeg);
	// 		// 	}
	// 		// }
	// 	}
	// 	std::cout<<"count: "<<2*count<<std::endl;
	// 	std::cout<<"direct: "<<2*ndirect<<std::endl;
	// 	Double_t Solidangle=PMSolidAngleMPPC(Origin,MPPCCenter,MPPCNorm);
	// 	std::cout<<"solid angle: "<<Solidangle<<std::endl;
	// 	grMCSolid    -> SetPoint(grMCSolid->GetN(),MPPCtheta*radtodeg,ndirect/(double)Npho/2.);
	// 	grSolidAngle -> SetPoint(grSolidAngle->GetN(),MPPCtheta*radtodeg,Solidangle);
	// 	grNphe ->SetPoint(grNphe->GetN(),MPPCtheta*radtodeg,count/(double)Npho/2./Solidangle);
	// 	grTS   ->SetPoint(grTS->GetN()  ,MPPCtheta*radtodeg,TransmittanceOfSilicon(MPPCtheta*radtodeg)/TransmittanceOfSilicon(0));
	// }
	// // hXY->Draw("colz");
	// canvas1->cd();
	// grSolidAngle->SetMinimum(0);
	// grSolidAngle->SetMarkerColor(2);
	// grSolidAngle->SetMarkerStyle(20);
	// grSolidAngle->Draw("ap");
	// grMCSolid->SetMarkerColor(4);
	// grMCSolid->SetMarkerStyle(20);
	// grMCSolid->Draw("same p");
	// canvas2->cd();
	// grNphe->SetMarkerColor(2);
	// grNphe->SetMarkerStyle(20);
	// grNphe->Draw("ap");
	// grTS->SetMarkerColor(4);
	// grTS->SetMarkerStyle(20);
	// grTS->Draw("same p");
	// grADPDE->SetTitle("Fresnel Equation;Incident Angle[deg];Transmission Probability");
	// grADPDE->SetMarkerColor(4);
	// grADPDE->SetMarkerStyle(20);
	// grADPDE->Draw("ap");

	// grTS->SetTitle("Fresnel Equation;Incident Angle[deg];Transmission Probability");
	// grTS->SetMarkerColor(4);
	// grTS->SetMarkerStyle(20);
	// grTS->SetMinimum(0);
	// grTS->Draw("ap");
	// htheta->Draw();
}

Double_t PMSolidAngleMPPC(TVector3& view, TVector3& center,TVector3& normal)
{
	// This function returns solid angle of MPPC.

	// check direction.
	TVector3 center_view = center - view;
	if ((center_view) * normal <= 0)
	return 0;

	TVector3 center_chip;
	// Double_t rtnval=0;

	//Set U,V direction. Temporary hard corded.
	TVector3 unit[3]; // unit vectors.
	unit[0].SetXYZ(0,1,0); //U direction
	unit[1] = (unit[0].Cross(normal)).Unit();//V direction
	unit[2] = normal.Unit();//W direction

	Double_t ChipSize=12;
	Double_t sin_a1;
	Double_t sin_a2;
	Double_t sin_b1;
	Double_t sin_b2;
	Double_t solid_total2=0;

	TVector3 vcorner1 = center + ChipSize/2. * unit[0] + ChipSize/2. * unit[1];
	TVector3 vcorner2 = center -ChipSize/2. * unit[0] - ChipSize/2. * unit[1];
	TVector3 v1       = view-vcorner1;
	TVector3 v2       = view-vcorner2;
	sin_a1 = v1.Dot(unit[0]) / sqrt(pow(v1.Dot(unit[0]),2) + pow(v1.Dot(unit[2]),2));
	sin_a2 = v2.Dot(unit[0]) / sqrt(pow(v2.Dot(unit[0]),2) + pow(v2.Dot(unit[2]),2));
	sin_b1 = v1.Dot(unit[1]) / sqrt(pow(v1.Dot(unit[1]),2) + pow(v1.Dot(unit[2]),2));
	sin_b2 = v2.Dot(unit[1]) / sqrt(pow(v2.Dot(unit[1]),2) + pow(v2.Dot(unit[2]),2));
	solid_total2 += TMath::Abs((asin(sin_a1*sin_b1) + asin(sin_a2*sin_b2) - asin(sin_a1*sin_b2) - asin(sin_b1*sin_a2)) / (4*TMath::Pi()));
	return solid_total2;
}

Double_t CalculateReflectivity_P (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t theta) {
	//Calculate reflectivity for P wave based on Fresnel equation
	theta *= degree;
	return
	pow(std::abs( n2*n2*cos(theta) - n1*std::sqrt(n2*n2-n1*n1*sin(theta)*sin(theta)) ), 2) /
	pow(std::abs( n2*n2*cos(theta) + n1*std::sqrt(n2*n2-n1*n1*sin(theta)*sin(theta)) ), 2);
}

Double_t CalculateReflectivity_S (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t theta) {
	//Calculate reflectivity for S wave based on Fresnel equation
	theta *= degree;
	return
	pow(std::abs( n1*cos(theta) -std::sqrt(n2*n2 - n1*n1*sin(theta)*sin(theta)) ), 2) /
	pow(std::abs( n1*cos(theta) +std::sqrt(n2*n2 - n1*n1*sin(theta)*sin(theta)) ), 2);
}

Double_t TransmittanceOfSilicon (std::complex<Double_t> n1, std::complex<Double_t> n2,Double_t theta){

	return 1.- CalculateReflectivity_P(n1,n2,theta)/2.0 - CalculateReflectivity_S(n1,n2,theta)/2.0;
}





Double_t PMIncidentAngle(TVector3& view, TVector3& center,TVector3& normal)
{
	TVector3 center_view = center - view;
	Double_t incidentAngle = center_view.Angle(normal);
	return incidentAngle;
}

Double_t PMAngleDependenceOfPDE(Double_t theta){
	/*
	Angle dependence of PDE
	Measured angle dependence of PDE in MPPC is icorrected with
	optical transmittance from LXe to silicon
	which is calculated from complex refractive index
	*/

	Double_t par[3]={0.15,21.45,-6.1e-5};//LP result
	Double_t rtnval=0;
	if(theta<0 ||theta>90)return 0;

	Double_t trans=TransmittanceOfSilicon(nxe,nsi,theta*TMath::DegToRad());
	if(trans<1e-4)return 0.;

	if(theta<par[1]){
		rtnval=1.;
	}else{
		rtnval=(par[0]+(theta-par[1])*(theta-par[1])*par[2])/par[0];

	}
	rtnval*=TransmittanceOfSilicon(nxe,nsi,0)/trans;
	// std::cout<<"angle: "<<theta<<" trans: "<<trans<<" rtnval: "<<rtnval<<std::endl;

	if(rtnval>1){
		return 1;
	}else if(rtnval<1e-4){
		return 0;
	}else{
		return rtnval;
	}
}

Bool_t JudgeDetection(TVector3& Direction, TVector3& MPPCCenter){
	Bool_t rtnbool = false;
	Double_t Scale=0;
	Double_t MPPCXMin =MPPCCenter.X()-MPPCSize/2.;
	Double_t MPPCXMax =MPPCCenter.X()+MPPCSize/2.;
	Double_t MPPCYMin =MPPCCenter.Y()-MPPCSize/2.;
	Double_t MPPCYMax =MPPCCenter.Y()+MPPCSize/2.;
	if(Direction.Z()>0){
		Scale= MPPCCenter.Z()/Direction.Z();
		Double_t XatMPPC  = Direction.X()*Scale;
		Double_t YatMPPC  = Direction.Y()*Scale;
		// Double_t incangle = PMIncidentAngle(Origin,MPPCCenter,MPPCNorm);
		if (XatMPPC>MPPCXMin&&XatMPPC<MPPCXMax){
			if (YatMPPC>MPPCYMin&&YatMPPC<MPPCYMax){
				rtnbool= true;
			}
		}
	}
	return rtnbool;

	// count=count+PMAngleDependenceOfPDE(incangle*radtodeg);
}