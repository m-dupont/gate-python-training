#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TF2.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom2.h>
#include <TFile.h>

void AnalyseRootFile() 
{	
  
  //NaI detector resolution (experimental)
  Double_t a = -1e-05;
  Double_t b = 0.08;
  Double_t c = 6.;
  TRandom2 *ran       = new TRandom2();

  
  //GATE
  //Get actor results
  TFile *f1 = new TFile("../output/Energy.root","read");
  TH1D *hEdep      = (TH1D *) f1->Get("edepHisto");

  //Create new histograms
  TH1F *hEdepCorr = new TH1F("hEdepCorr","hEdepCorr",2000,0,2);
  TH1D *hEdepCanal = new TH1D("hEdepCanal","hEdepCanal",1024,0,1024);
  TH1D *hCanal = new TH1D("hCanal","hCanal",1024,0,1024);

  //Add detector resolution and convert into canal
  for(int i=1;i < hEdep->GetNbinsX(); i++)
    {
      Double_t nE = hEdep->GetBinContent(i);
      Double_t Edep = i*0.5;

      //Loop on histograms bins
      for(int j = 0;j < nE; j++)
	{
	  
	  //Add detector resolution
	  Double_t sigma_edep = (a*(Edep*Edep)+b*Edep+c)/2.37; //FWHM = 2.37xsigma = axE^2+bxE+c
	  Double_t Edepcorr =  ran->Gaus(Edep,sigma_edep) ;

	  //Convert into canal
	  Double_t Edepcanal = (Edepcorr+12.)/3.18; // calibration linéaire
	  	  
	  //Fill histograms
	  if(i>10)
	    {
	      hEdepCanal->Fill(Edepcanal); 
	      hEdepCorr->Fill(Edepcorr/1000.);
	    }
	}
    }

  //Genie2000 (experimental measurements)  
  TH1F *hGenie2K = new TH1F("hGenie2K","hGenie2K",1024,0,1024);
  string parametersline;
  std::ifstream fin("137Cs_2cm_5min_bdfsoustr.txt");
  Int_t count = 0.;
  int nbin = 0;
  while(getline(fin,parametersline))
    {
      sscanf(parametersline.c_str(),"%d",&count);
      if(nbin > 10) hGenie2K->Fill(nbin,count);
      nbin+=1;
    }

  //Normalize to maximum
  hCanal->Scale(hGenie2K->GetMaximum()/hCanal->GetMaximum());
  hEdepCanal->Scale(hGenie2K->GetMaximum()/hEdepCanal->GetMaximum());
		
  //Draw histograms	
  gStyle->SetOptStat(0000);
  TCanvas *c2 = new TCanvas("c2","c2",900,400);
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLogy(1);
  hEdep->SetTitle("");
  hEdep->SetXTitle("Energy [MeV]");
  hEdep->SetMarkerStyle(20);
  hEdep->SetMarkerColor(2);
  hEdep->SetLineColor(2);
  hEdep->Draw("HIST");
  hEdepCorr->SetMarkerStyle(20);
  hEdepCorr->SetMarkerColor(4);
  hEdepCorr->SetLineColor(4);
  hEdepCorr->Draw("HIST SAME");
  TLegend* leg = new TLegend(0.55,0.75,0.9,0.9);
  leg->SetLineStyle(1);
  leg->SetFillStyle(0);
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetTextSize(0.035);
  leg->AddEntry(hEdep,"w/o energy resolution","pl");
  leg->AddEntry(hEdepCorr,"with energy resolution","pl");
  leg->Draw("same");
  c2->cd(2);
  gPad->SetLogy(1);
  hGenie2K->SetTitle("");
  hGenie2K->SetXTitle("Canal");
  hGenie2K->SetMarkerStyle(20);
  hGenie2K->SetMarkerColor(1);
  hGenie2K->SetLineColor(1);
  hGenie2K->SetLineWidth(2);
  hGenie2K->Draw("HIST");
  hEdepCanal->SetMarkerStyle(20);
  hEdepCanal->SetMarkerColor(4);
  hEdepCanal->SetLineColor(4);
  hEdepCanal->SetLineWidth(2);
  hEdepCanal->Draw("HIST SAME");
  TLegend* leg2 = new TLegend(0.65,0.75,0.9,0.9);
  leg2->SetLineStyle(1);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(1);
  leg2->SetLineColor(1);
  leg2->SetTextSize(0.035);
  leg2->AddEntry(hGenie2K,"Data","pl");
  leg2->AddEntry(hEdepCanal,"GATE","pl");
  leg2->Draw("same");
		
}
