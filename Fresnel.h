const Double_t degree = TMath::Pi()/180.;
const Double_t nm = 1./TMath::Power(10,9);
const Double_t microm = 1./TMath::Power(10,6);
Double_t FresnelReflection_P(Double_t* x,Double_t* par);
Double_t FresnelReflection_S(Double_t* x,Double_t* par);
Double_t FresnelReflection_T(Double_t* x,Double_t* par);
Double_t ConventionalR_P(Double_t* x,Double_t* par);
Double_t ConventionalR_S(Double_t* x,Double_t* par);
Double_t ConventionalR_T(Double_t* x,Double_t* par);
std::complex<Double_t> FresnelCoefficientT_P (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle);
std::complex<Double_t> FresnelCoefficientT_S (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle);
std::complex<Double_t> FresnelCoefficientR_P (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle);
std::complex<Double_t> FresnelCoefficientR_S (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle);
std::complex<Double_t> FresnelRefAmplitude_P(Double_t* x,Double_t* par);
std::complex<Double_t> FresnelRefAmplitude_S(Double_t* x,Double_t* par);
std::complex<Double_t> RelativePermittivity(std::complex<Double_t> n1);
Double_t RefractionAngle(std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t theta1);
Double_t ZeroCenter(Double_t value);
Double_t PhaseShift(Double_t *x,Double_t* par);

Double_t Trans(Double_t* x, Double_t* par){
   Double_t z[1]={0};
   Double_t y[1]={x[0]*TMath::DegToRad()};
   return (2-(FresnelReflection_P(y,par)+FresnelReflection_S(y,par)))/
   (2-(FresnelReflection_P(z,par)+FresnelReflection_S(z,par)));
}

Double_t ZeroCenter(Double_t value){
   if (value>TMath::Pi()/2.)
   {
      return value-TMath::Pi();
   }else{
      return value;
   }
}

Double_t ConventionalR_P(Double_t* x,Double_t* par){
   Double_t theta=x[0];
   std::complex<Double_t> nLXe  ={par[0],par[1]};
   std::complex<Double_t> nLayer={par[2],par[3]};
   std::complex<Double_t> nSi   ={par[4],par[5]};
   // Double_t dLayer = par[6];
   std::complex<Double_t> rP_01 = FresnelCoefficientR_P(nLXe  ,nLayer,theta);
   Double_t refangle = RefractionAngle(nLXe,nLayer,theta);
   if (refangle<0) return 1.;
   std::complex<Double_t> rP_12 = FresnelCoefficientR_P(nLayer,nSi   ,refangle);
   std::complex<Double_t> rP_10 = FresnelCoefficientR_P(nLayer,nLXe   ,refangle);
   Double_t RP_01=norm(rP_01);
   Double_t RP_12=norm(rP_12);
   Double_t RP_10=norm(rP_10);
   // Double_t delta = 4*TMath::Pi()*real(nLayer)*dLayer*TMath::Cos(refangle)/wavelength;
   // std::complex<Double_t> midelta = {0,-delta};
   // // std::complex<Double_t> midelta = {0,0};
   // std::complex<Double_t> One = {1,0};
   Double_t Rtotal=RP_01+(1-RP_01)*RP_12*(1-RP_10)/(1-RP_10*RP_12);
   // std::complex<Double_t> Amplitude_P = (rP_01+rP_12*exp(midelta))/(One+rP_01*rP_12*exp(midelta));
   // std::complex<Double_t> Amplitude_P = rP_01+rP_12*(One-rP_01*rP_01)+rP_01*rP_12*rP_12*(One-rP_01*rP_01);
   // std::cout<<"angle: "<<theta/degree<<std::endl;
   // std::cout<<" RP01: "<<real(rP_01)<<"+i"<<imag(rP_01)<<std::endl;
   // std::cout<<" RP12: "<<real(rP_12)<<"+i"<<imag(rP_12)<<std::endl;
   return Rtotal;
}

Double_t ConventionalR_S(Double_t* x,Double_t* par){
   Double_t theta=x[0];
   std::complex<Double_t> nLXe  ={par[0],par[1]};
   std::complex<Double_t> nLayer={par[2],par[3]};
   std::complex<Double_t> nSi   ={par[4],par[5]};
   std::complex<Double_t> rS_01 = FresnelCoefficientR_S(nLXe  ,nLayer,theta);
   Double_t refangle = RefractionAngle(nLXe,nLayer,theta);
   if (refangle<0) return 1.;
   std::complex<Double_t> rS_12 = FresnelCoefficientR_S(nLayer,nSi   ,refangle);
   std::complex<Double_t> rS_10 = FresnelCoefficientR_S(nLayer,nLXe   ,refangle);
   Double_t RS_01=norm(rS_01);
   Double_t RS_12=norm(rS_12);
   Double_t RS_10=norm(rS_10);
   Double_t Rtotal=RS_01+(1-RS_01)*RS_12*(1-RS_10)/(1-RS_10*RS_12);
   return Rtotal;
}

Double_t ConventionalR_T(Double_t* x,Double_t* par){
   return (ConventionalR_P(x,par)+ConventionalR_S(x,par))/2.;
}

Double_t FresnelReflection_T(Double_t* x,Double_t* par){
   return (FresnelReflection_P(x,par)+FresnelReflection_S(x,par))/2.;
}

Double_t FresnelReflection_P(Double_t* x,Double_t* par){
   return norm(FresnelRefAmplitude_P(x,par));
}

Double_t FresnelReflection_S(Double_t* x,Double_t* par){
   return norm(FresnelRefAmplitude_S(x,par));
}

Double_t FresnelPhaseShift_P(Double_t* x,Double_t* par){
return arg(FresnelRefAmplitude_P(x,par));
}

Double_t FresnelPhaseShift_S(Double_t* x,Double_t* par){
return arg(FresnelRefAmplitude_S(x,par));
}

std::complex<Double_t> FresnelRefAmplitude_P(Double_t* x,Double_t* par){
   Double_t theta=x[0];
   std::complex<Double_t> nLXe  ={par[0],par[1]};
   std::complex<Double_t> nLayer={par[2],par[3]};
   std::complex<Double_t> nSi   ={par[4],par[5]};
   Double_t wavelength = par[6];
   Double_t dLayer     = par[7];
   // std::cout<<"fp: parameters"<<std::endl;
   std::complex<Double_t> rP_01 = FresnelCoefficientR_P(nLXe  ,nLayer,theta);
   std::complex<Double_t> tP_01 = FresnelCoefficientT_P(nLXe  ,nLayer,theta);
   
   Double_t refangle = RefractionAngle(nLXe,nLayer,theta);
   if (refangle<0) return 1.;
   std::complex<Double_t> rP_12 = FresnelCoefficientR_P(nLayer,nSi   ,refangle);
   std::complex<Double_t> tP_10 = FresnelCoefficientT_P(nLayer,nLXe  ,refangle);
   std::complex<Double_t> tP_12 = FresnelCoefficientT_P(nLayer,nSi   ,refangle);
   Double_t delta = PhaseShift(x,par);
   // Double_t delta = 0;
   std::complex<Double_t> midelta = {0,-delta};
   std::complex<Double_t> One = {1.,0};
   std::complex<Double_t> Amplitude_P = (rP_01+rP_12*exp(midelta))/(One+rP_01*rP_12*exp(midelta));
   return Amplitude_P;
}

std::complex<Double_t> FresnelRefAmplitude_S(Double_t* x,Double_t* par){
   Double_t theta=x[0];
   std::complex<Double_t> nLXe  ={par[0],par[1]};
   std::complex<Double_t> nLayer={par[2],par[3]};
   std::complex<Double_t> nSi   ={par[4],par[5]};
   Double_t wavelength = par[6];
   Double_t dLayer     = par[7];
   // std::cout<<"fp: parameters"<<std::endl;
   std::complex<Double_t> rS_01 = FresnelCoefficientR_S(nLXe  ,nLayer,theta);
   Double_t refangle = RefractionAngle(nLXe,nLayer,theta);
   std::complex<Double_t> rS_12 = FresnelCoefficientR_S(nLayer,nSi   ,refangle);
   Double_t delta = PhaseShift(x,par);
   // Double_t delta = 0;
   std::complex<Double_t> midelta = {0,-delta};
   std::complex<Double_t> One = {1.,0};
   std::complex<Double_t> Amplitude_S = (rS_01+rS_12*exp(midelta))/(One+rS_01*rS_12*exp(midelta));
   // std::complex<Double_t> Amplitude_S = rS_12;
   // Int_t Nreflection=200;
   // std::complex<Double_t> Amplitude_S=rS_01;
   // std::complex<Double_t> BaseTerm=rS_12*(One-rS_01*rS_01)*exp(midelta);
   // for (Int_t iref = 0; iref < Nreflection; iref++) {
   //    Amplitude_S+=BaseTerm;
   //    BaseTerm *= rS_01*rS_12;
   // }
   // std::complex<Double_t> Amplitude_S = rS_01+rS_12*(One-rS_01*rS_01)+rS_01*rS_12*rS_12*(One-rS_01*rS_01);
   return Amplitude_S;
   
}


Double_t PhaseShift(Double_t *x,Double_t* par){
   Double_t incidentAngle=x[0];
   std::complex<Double_t> nLXe  ={par[0],par[1]};
   std::complex<Double_t> nLayer={par[2],par[3]};
   std::complex<Double_t> nSi   ={par[4],par[5]};
   Double_t wavelength = par[6];
   Double_t dLayer     = par[7];
   
   Double_t ra = RefractionAngle(nLXe,nLayer,incidentAngle);
   return 4.*TMath::Pi()*real(nLayer)*dLayer*TMath::Cos(ra)/wavelength;
}
// std::complex<Double_t> delta(std::complex<Double_t> n,Double_t thick){
// return 4*TMath::Pi()*n*thick/wavelength;
// }


std::complex<Double_t> FresnelCoefficientR_P (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle) {
   //Calculate reflectivity for P wave based on Fresnel equation
   std::complex<Double_t> RP1=RelativePermittivity(n1);
   std::complex<Double_t> RP2=RelativePermittivity(n2);
   return 
   (RP2*cos(incangle)-sqrt(RP1*(RP2-RP1*sin(incangle)*sin(incangle))))/
   (RP2*cos(incangle)+sqrt(RP1*(RP2-RP1*sin(incangle)*sin(incangle))));
}

std::complex<Double_t> FresnelCoefficientR_S (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle) {
   //Calculate reflectivity for S wave based on Fresnel equation
   std::complex<Double_t> RP1=RelativePermittivity(n1);
   std::complex<Double_t> RP2=RelativePermittivity(n2);
   return 
   (sqrt(RP1)*cos(incangle)-sqrt(RP2-RP1*sin(incangle)*sin(incangle)))/
   (sqrt(RP1)*cos(incangle)+sqrt(RP2-RP1*sin(incangle)*sin(incangle)));
}

std::complex<Double_t> FresnelCoefficientT_P (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle) {
   //Calculate reflectivity for P wave based on Fresnel equation
   std::complex<Double_t> RP1=RelativePermittivity(n1);
   std::complex<Double_t> RP2=RelativePermittivity(n2);
   return 
   2.*sqrt(RP1)*sqrt(RP2)*cos(incangle)/
   (RP2*cos(incangle)+sqrt(RP1*(RP2-RP1*sin(incangle)*sin(incangle))));
}

std::complex<Double_t> FresnelCoefficientT_S (std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t incangle) {
   //Calculate reflectivity for S wave based on Fresnel equation
   std::complex<Double_t> RP1=RelativePermittivity(n1);
   std::complex<Double_t> RP2=RelativePermittivity(n2);
   return 
   2.*sqrt(RP1)*cos(incangle)/
   (sqrt(RP1)*cos(incangle)+sqrt(RP2-RP1*sin(incangle)*sin(incangle)));
}

Double_t RefractionAngle(std::complex<Double_t> n1, std::complex<Double_t> n2, Double_t theta1){
   Double_t sinref=real(n1)*TMath::Sin(theta1)/real(n2);
   if (sinref>1.)
   {
      return TMath::Pi()/2.;
   }
   return TMath::ASin(sinref);
}

std::complex<Double_t> RelativePermittivity(std::complex<Double_t> n1){
   return n1*n1;
}