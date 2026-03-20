// Eloss.h
#ifndef EneLoss_H
#define EneLoss_H

#include <TMath.h>
#include <TRandom3.h>
#include "nucleus.h"

enum class StringCode
{
		Foil,	
		Al,	
		B,	
		C4H10,	
		CsI,
		Target,	
		Mylar,	
		P,	
		Si,	
		Si3N4,	
		SiO2,
        unknown	
};


StringCode hashString(const TString& str);

Double_t* getElossE(nucleus&, TString);
Double_t* getElossDeDx(nucleus&, TString);

Double_t eval(Double_t, Double_t[100], Double_t[100]);
Double_t eloss(nucleus, Double_t, Double_t, Double_t, Double_t[100], Double_t[100]);
Double_t eloss(nucleus, Double_t, Double_t, Double_t, TString);
Double_t eloss_Lise(nucleus, Double_t, Double_t, Double_t, Double_t[100], Double_t[100]);
//Double_t eloss_CATIMA(nucleus, Double_t, Double_t, Double_t, Double_t[100], Double_t[100]);
Double_t elossFi(Double_t, Double_t, Double_t[100], Double_t[100]);
#endif
// end
