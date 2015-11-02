#define HH4bAna_v2_cxx
#include "HH4bAna_v2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <iostream>

void HH4bAna_v2::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L HH4bAna_v2.C
//      Root > HH4bAna_v2 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  TFile* Output = new TFile("Outfile.root", "RECREATE");
  Output->cd();

  TH1D* hnjets = new TH1D("hnjets","N jets", 6, -0.5, 5.5);
  TH1D* hmassGR = new TH1D("hmassGR","Graviton Mass", 60, 500, 3500);
  TH1D* hdetaGR = new TH1D("hdetaGR","Graviton cos theta*", 20, 0, 2);

  TH1D* hmassJet0 = new TH1D("hmassJet1","Jet1 Pruned Mass", 40, 0, 200);
  TH1D* hmassJet1 = new TH1D("hmassJet2","Jet2 Pruned Mass", 40, 0, 200);


   // Jet Id
   TH1D* h_jetAK8_nhf    = new TH1D("h_jetAK8_nhf"   , "Neutral hadron fraction", 40, 0.   , 2.) ;
   TH1D* h_jetAK8_chf    = new TH1D("h_jetAK8_chf"   , "Charged hadron fraction", 40, 0.   , 2.) ;
   TH1D* h_jetAK8_phf    = new TH1D("h_jetAK8_phf"   , "Photon fraction"    , 40, 0.   , 2.) ;
   TH1D* h_jetAK8_emf    = new TH1D("h_jetAK8_cef"   , "Charged em fraction"    , 40, 0.   , 2.) ;
   TH1D* h_jetAK8_muf    = new TH1D("h_jetAK8_muf"   , "Muon fraction"          , 40, 0.  , 2.) ;
   TH1D* h_jetAK8_nconst = new TH1D("h_jetAK8_nconst", "No. of constituents"    , 40, -0.5, 40.5) ;

   TH1D* h_nsj_GR = new TH1D("h_nsj_GR", "No. of subjets in 2 leading jets"    , 7, -0.5, 6.5) ;
   TH1D* h_nsj_csvl_GR = new TH1D("h_nsj_csvl_GR", "No. of subjets in 2 leading jets that are b-tagged"    , 7, -0.5, 6.5) ;
   
   TH1D* h_jetAK8_tau21_0   = new TH1D("h_jetAK8_tau21_0"   , "tau21 for jet 0", 10, 0., 1.) ;
   TH1D* h_jetAK8_tau21_1   = new TH1D("h_jetAK8_tau21_1"   , "tau21 for jet 1", 10, 0., 1.) ;

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      double  njetsAK8 = jetAK8_N;
      if (jentry%1000 == 0) std::cout << "njets = " << njetsAK8 << std::endl;
	
      if (njetsAK8 < 2) continue;
     


     TLorentzVector jetp4_1 ;
     jetp4_1.SetPtEtaPhiM(jetAK8_pt->at(0), jetAK8_eta->at(0), jetAK8_phi->at(0), jetAK8_mass->at(0));
     TLorentzVector jetp4_2 ;
     jetp4_2.SetPtEtaPhiM(jetAK8_pt->at(1), jetAK8_eta->at(1), jetAK8_phi->at(1), jetAK8_mass->at(1));

     bool b_jetId_0 = jetAK8_IDTight->at(0);
     bool b_jetId_1 = jetAK8_IDTight->at(1);
     bool b_jetId = b_jetId_0 && b_jetId_1;
     if (jentry%1000 == 0) cout << b_jetId_0 << endl;

     bool b_phase_space_jets = ( fabs(jetp4_1.Eta()) < 2.4 ) && ( fabs(jetp4_2.Eta()) < 2.4 ) &&  ( jetp4_1.Pt() > 200 ) && ( jetp4_2.Pt() > 200);

     TLorentzVector jetp4_GR = jetp4_1 + jetp4_2;

     double dEta = jetp4_1.Eta() - jetp4_2.Eta();
     double mass = jetp4_GR.M();

     bool b_CosThetaStar = fabs(dEta) < 1.3; 

     if ( !b_jetId || !b_phase_space_jets || !b_CosThetaStar ) continue;

     double m0 = jetAK8_pruned_mass->at(0);
     double m1 = jetAK8_pruned_mass->at(1);

     bool massSelect0 = m0 < 145 && m0 > 105;
     bool massSelect1 = m1 < 145 && m1 > 105;

     int nsj_csvl = 0;
     int nsj = 0;

     for (int ijet = 0; ijet < 2; ijet++){
       if (subjetAK8_softdrop_N->at(ijet) != 2) continue;
       for (int isj = 0; isj < subjetAK8_softdrop_N->at(ijet); isj++){
	 vector<float> vcsv =  subjetAK8_softdrop_csv->at(ijet);
	 float csv = vcsv[isj];
	 if (csv > 0.605) nsj_csvl++;
	 nsj++;
       }
     }

     double tau21_0 = jetAK8_tau2->at(0)/jetAK8_tau1->at(0);
     double tau21_1 = jetAK8_tau2->at(1)/jetAK8_tau1->at(1);

     h_nsj_GR->Fill(nsj);
     h_nsj_csvl_GR->Fill(nsj_csvl);
     h_jetAK8_tau21_0->Fill(tau21_0);
     h_jetAK8_tau21_1->Fill(tau21_1);

     // fill Histos
     //   if (!massSelect0 || !massSelect1) continue;

     hnjets->Fill(njetsAK8);
     hdetaGR->Fill(fabs(dEta));
     hmassGR->Fill(jetp4_GR.M());
       
     for (int ijet = 0; ijet < 2; ++ijet) {
       h_jetAK8_nhf -> Fill(jetAK8_nhf->at(ijet)) ;
       h_jetAK8_chf -> Fill(jetAK8_chf->at(ijet)) ;
       h_jetAK8_emf -> Fill(jetAK8_emf->at(ijet)) ;
       h_jetAK8_phf -> Fill(jetAK8_phf->at(ijet)) ;
       h_jetAK8_muf -> Fill (jetAK8_muf->at(ijet)) ;
       h_jetAK8_nconst -> Fill (jetAK8_cm->at(ijet) + jetAK8_nm->at(ijet)) ;
       hmassJet0 -> Fill(m0);
       hmassJet1 -> Fill(m1);
     }


   }

   TCanvas* Canv = new TCanvas();
   Canv->cd();
   hnjets->Draw();
   Canv->SaveAs("njets.png");

   Canv->Clear();
   hmassGR->Draw();
   Canv->SaveAs("massGR.png");

   Canv->Clear();
   hdetaGR->Draw();
   Canv->SaveAs("detaGR.png");
   
   Output->Write();
   Output->Close();


}


