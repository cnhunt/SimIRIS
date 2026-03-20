#include "eloss.h"
#include "catima/catima.h"

Double_t emptyArray[100] = {};

StringCode hashString(const TString& str)
{
	if (str == "Foil") return StringCode::Foil;
	if (str == "Al") return StringCode::Al;
	if (str == "B") return StringCode::B;
	if (str == "C4H10") return StringCode::C4H10;
	if (str == "CsI") return StringCode::Foil;
	if (str == "Target") return StringCode::Target;
	if (str == "Mylar") return StringCode::Mylar;
	if (str == "P") return StringCode::P;
	if (str == "Si") return StringCode::Si;
	if (str == "Si3N4") return StringCode::Si3N4;
	if (str == "SiO2") return StringCode::SiO2;
	return StringCode::unknown;
}

Double_t* getElossE(nucleus &P, TString material)
{
	switch(hashString(material)) {
		case StringCode::Foil:
			return P.EL.eFoil;
		case StringCode::Al:
			return P.EL.eAl;
		case StringCode::B:
			return P.EL.eB;
		case StringCode::C4H10:
			return P.EL.eC4H10;
		case StringCode::CsI:
			return P.EL.eCsI;
		case StringCode::Target:
			return P.EL.eTgt;
		case StringCode::Mylar:
			return P.EL.eMy;
		case StringCode::P:
			return P.EL.eP;
		case StringCode::Si:
			return P.EL.eSi;
		case StringCode::Si3N4:
			return P.EL.eSi3N4;
		case StringCode::SiO2:
			return P.EL.eSiO2;
		default:
			return emptyArray;
	}
}

Double_t* getElossDeDx(nucleus &P, TString material)
{
	switch(hashString(material)) {
		case StringCode::Foil:
			return P.EL.dedxFoil;
		case StringCode::Al:
			return P.EL.dedxAl;
		case StringCode::B:
			return P.EL.dedxB;
		case StringCode::C4H10:
			return P.EL.dedxC4H10;
		case StringCode::CsI:
			return P.EL.dedxCsI;
		case StringCode::Target:
			return P.EL.dedxTgt;
		case StringCode::Mylar:
			return P.EL.dedxMy;
		case StringCode::P:
			return P.EL.dedxP;
		case StringCode::Si:
			return P.EL.dedxSi;
		case StringCode::Si3N4:
			return P.EL.dedxSi3N4;
		case StringCode::SiO2:
			return P.EL.dedxSiO2;
		default:
			return emptyArray;
	}
}


Double_t eval(Double_t in, Double_t x[100], Double_t y[100])
{
	Double_t dxin=0., dx=0., dy=0., e=0.;
	if(in<=0.){
		e = 0.;
	}
	else if(in<x[0]){
		e = y[0]*in/x[0];
	}
	else if(in>x[99]){
		dxin = in-x[99];
		dx = x[99]-x[98];
		dy = y[99]-y[98];
		e = y[99]+dy*dxin/dx;
	}
	else{
		for(int i=1; i<100;i++){   
			if(in>x[i-1]&&in<x[i]){
				dxin = in-x[i-1];
				dx = x[i]-x[i-1];
				dy = y[i]-y[i-1];
				e = y[i-1]+dy*dxin/dx;
				break;
			}
		}
	}
	return e;
}

//Make it a method for a particle class.
Double_t eloss(nucleus P, Double_t TZoverA, Double_t ein, Double_t th , Double_t x[100], Double_t y[100])//initial energy and thickness are given as arguments 
{

	return eloss_Lise(P, TZoverA, ein, th, x, y);
}

Double_t eloss(nucleus P, Double_t TZoverA, Double_t ein, Double_t th , TString material)//initial energy and thickness are given as arguments 
{

	return eloss_Lise(P, TZoverA, ein, th, getElossE(P, material), getElossDeDx(P, material));
}

Double_t eloss_Lise(nucleus P, Double_t TZoverA, Double_t ein, Double_t th , Double_t x[100], Double_t y[100])//initial energy and thickness are given as arguments 
{

	Double_t k;
	Double_t Bohr;
	Double_t sgm;
	Double_t strg;

	if(ein==0.) return 0;
	if(th==0.) return 0.;
	
	TRandom3 *rndm = new TRandom3(0);
	//Energy loss calculation including energy Straggling
	Double_t dx =th/100.; //in mg/cm2
	Bohr = TMath::Sqrt(0.157 * dx * P.Z*P.Z *TZoverA)/1000.;
	Double_t de = 0; //energy loss
	Double_t en= ein; //the energy variable
	for (int i=0; i<100; i++){
	  	de = (dx * eval(en,x,y));//energy loss in dx
		if(de>en){
		   	en=0.;	
			break;
		}
		k= 1.1+0.47*TMath::Log10(en/Double_t(P.A));
		sgm = k*Bohr;
		strg = rndm->Gaus(de,sgm);
		de = (strg>0.) ? strg : 0.;
	  	if(de>en){
		   	en=0.;	
			break;
		}
		en = en - de; // energy remaining after dx
	}
	rndm->Delete();
	return ein-en;
}

Double_t elossFi(Double_t efi, Double_t th, Double_t x[100], Double_t y[100])//final energy and thickness are given as arguments 
{
	if (th==0) return efi;
	//Energy loss calculation
	Double_t dx =th/100.; //in 
	Double_t de = 0; //energy loss
	Double_t en= efi; //the energy variable
	Double_t pos = 0.;// the position variable
	
	while (pos<= th){
		de =  (dx* eval(en,x,y))/2.; //
		en = en + de;
		de =(dx* eval(en,x,y))/2.;//energy loss in dx
		en = en +de; // energy remaining after dx
		pos = pos +dx;
	}
 	return en-efi;
}
