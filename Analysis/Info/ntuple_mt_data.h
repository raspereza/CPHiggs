//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 23 08:49:21 2026 by ROOT version 6.30/07
// from TTree ntuple/
// found on file: root://eoscms.cern.ch//eos/cms/store/group/phys_tau/lrussell/forAliaksei/ForFakeFactors/Run3_2022EE/mt/Muon_Run2022E/nominal/merged.root
//////////////////////////////////////////////////////////

#ifndef ntuple_h
#define ntuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ntuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Long64_t        event;
   Long64_t        run;
   Long64_t        lumi;
   Double_t        charge_1;
   Short_t         charge_2;
   Double_t        pt_1;
   Double_t        eta_1;
   Double_t        phi_1;
   Double_t        mass_1;
   Double_t        pt_2;
   Double_t        eta_2;
   Double_t        phi_2;
   Double_t        mass_2;
   Bool_t          os;
   Double_t        dR;
   Double_t        dphi;
   Double_t        pt_tt;
   Double_t        pt_vis;
   Double_t        phi_vis;
   Double_t        eta_vis;
   Double_t        mt_1;
   Double_t        mt_2;
   Double_t        mt_lep;
   Double_t        mt_tot;
   Double_t        m_vis;
   Double_t        met_pt;
   Double_t        met_phi;
   Double_t        met_covXX;
   Double_t        met_covXY;
   Double_t        met_covYY;
   Double_t        met_dphi_1;
   Double_t        met_dphi_2;
   Bool_t          trg_singlemuon;
   Bool_t          trg_mt_cross;
   Long64_t        idDeepTau2018v2p5VSjet_2;
   Long64_t        idDeepTau2018v2p5VSmu_2;
   Long64_t        idDeepTau2018v2p5VSe_2;
   Long64_t        decayMode_2;
   Short_t         decayModePNet_2;
   Double_t        n_jets;
   Double_t        n_bjets;
   Double_t        mjj;
   Double_t        jdeta;
   Double_t        sjdphi;
   Double_t        dijetpt;
   Double_t        jpt_1;
   Double_t        jeta_1;
   Double_t        jphi_1;
   Double_t        jpt_2;
   Double_t        jeta_2;
   Double_t        jphi_2;
   Double_t        aco_mu_pi;
   Double_t        aco_mu_rho;
   Double_t        aco_mu_a1_FASTMTT_MassConstraint;
   Double_t        alphaAngle_mu_pi_1;
   Double_t        alphaAngle_mu_pi_2;
   Double_t        alphaAngle_mu_rho_1;
   Double_t        alphaAngle_mu_rho_2;
   Double_t        alphaAngle_mu_a1_FASTMTT_MassConstraint_1;
   Double_t        alphaAngle_mu_a1_FASTMTT_MassConstraint_2;
   Double_t        Ostar_aco_mu_pi;
   Double_t        Ostar_aco_mu_rho;
   Double_t        Ostar_aco_mu_a1_FASTMTT_MassConstraint;
   Double_t        PV_x;
   Double_t        PV_y;
   Double_t        PV_z;
   Double_t        PVBS_x;
   Double_t        PVBS_y;
   Double_t        PVBS_z;
   Double_t        ip_x_1;
   Double_t        ip_y_1;
   Double_t        ip_z_1;
   Double_t        ip_x_2;
   Double_t        ip_y_2;
   Double_t        ip_z_2;
   Double_t        ip_LengthSig_1;
   Double_t        ip_LengthSig_2;
   Double_t        hasRefitSV_1;
   Bool_t          hasRefitSV_2;
   Double_t        sv_x_1;
   Double_t        sv_y_1;
   Double_t        sv_z_1;
   Double_t        sv_x_2;
   Double_t        sv_y_2;
   Double_t        sv_z_2;
   Double_t        FastMTT_mass;
   Double_t        FastMTT_pt;
   Double_t        FastMTT_pt_1;
   Double_t        FastMTT_pt_2;
   Double_t        FastMTT_mass_constraint;
   Double_t        FastMTT_pt_constraint;
   Double_t        FastMTT_pt_1_constraint;
   Double_t        FastMTT_pt_2_constraint;
   Double_t        pi_pt_1;
   Double_t        pi_eta_1;
   Double_t        pi_phi_1;
   Double_t        pi_mass_1;
   Double_t        pi_charge_1;
   Double_t        pi_pdgId_1;
   Double_t        pi_Energy_1;
   Double_t        pi2_pt_1;
   Double_t        pi2_eta_1;
   Double_t        pi2_phi_1;
   Double_t        pi2_mass_1;
   Double_t        pi2_charge_1;
   Double_t        pi2_pdgId_1;
   Double_t        pi2_Energy_1;
   Double_t        pi3_pt_1;
   Double_t        pi3_eta_1;
   Double_t        pi3_phi_1;
   Double_t        pi3_mass_1;
   Double_t        pi3_charge_1;
   Double_t        pi3_pdgId_1;
   Double_t        pi3_Energy_1;
   Double_t        pi0_pt_1;
   Double_t        pi0_eta_1;
   Double_t        pi0_phi_1;
   Double_t        pi0_mass_1;
   Double_t        pi0_charge_1;
   Double_t        pi0_pdgId_1;
   Double_t        pi0_Energy_1;
   Double_t        pi_pt_2;
   Double_t        pi_eta_2;
   Double_t        pi_phi_2;
   Double_t        pi_mass_2;
   Double_t        pi_charge_2;
   Double_t        pi_pdgId_2;
   Double_t        pi_Energy_2;
   Double_t        pi2_pt_2;
   Double_t        pi2_eta_2;
   Double_t        pi2_phi_2;
   Double_t        pi2_mass_2;
   Double_t        pi2_charge_2;
   Double_t        pi2_pdgId_2;
   Double_t        pi2_Energy_2;
   Double_t        pi3_pt_2;
   Double_t        pi3_eta_2;
   Double_t        pi3_phi_2;
   Double_t        pi3_mass_2;
   Double_t        pi3_charge_2;
   Double_t        pi3_pdgId_2;
   Double_t        pi3_Energy_2;
   Double_t        pi0_pt_2;
   Double_t        pi0_eta_2;
   Double_t        pi0_phi_2;
   Double_t        pi0_mass_2;
   Double_t        pi0_charge_2;
   Double_t        pi0_pdgId_2;
   Double_t        pi0_Energy_2;
   Double_t        weight;
   Double_t        is_data;
   Double_t        pion_E_split_1;
   Double_t        pion_E_split_2;
   Double_t        iso_1;
   Double_t        abs_eta_1;
   Double_t        BDT_pred_score;
   Double_t        BDT_pred_class;
   Double_t        BDT_raw_score_tau;
   Double_t        BDT_raw_score_higgs;
   Double_t        BDT_raw_score_fake;
   Double_t        BDT_W_score;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_charge_1;   //!
   TBranch        *b_charge_2;   //!
   TBranch        *b_pt_1;   //!
   TBranch        *b_eta_1;   //!
   TBranch        *b_phi_1;   //!
   TBranch        *b_mass_1;   //!
   TBranch        *b_pt_2;   //!
   TBranch        *b_eta_2;   //!
   TBranch        *b_phi_2;   //!
   TBranch        *b_mass_2;   //!
   TBranch        *b_os;   //!
   TBranch        *b_dR;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_pt_tt;   //!
   TBranch        *b_pt_vis;   //!
   TBranch        *b_phi_vis;   //!
   TBranch        *b_eta_vis;   //!
   TBranch        *b_mt_1;   //!
   TBranch        *b_mt_2;   //!
   TBranch        *b_mt_lep;   //!
   TBranch        *b_mt_tot;   //!
   TBranch        *b_m_vis;   //!
   TBranch        *b_met_pt;   //!
   TBranch        *b_met_phi;   //!
   TBranch        *b_met_covXX;   //!
   TBranch        *b_met_covXY;   //!
   TBranch        *b_met_covYY;   //!
   TBranch        *b_met_dphi_1;   //!
   TBranch        *b_met_dphi_2;   //!
   TBranch        *b_trg_singlemuon;   //!
   TBranch        *b_trg_mt_cross;   //!
   TBranch        *b_idDeepTau2018v2p5VSjet_2;   //!
   TBranch        *b_idDeepTau2018v2p5VSmu_2;   //!
   TBranch        *b_idDeepTau2018v2p5VSe_2;   //!
   TBranch        *b_decayMode_2;   //!
   TBranch        *b_decayModePNet_2;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_n_bjets;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_jdeta;   //!
   TBranch        *b_sjdphi;   //!
   TBranch        *b_dijetpt;   //!
   TBranch        *b_jpt_1;   //!
   TBranch        *b_jeta_1;   //!
   TBranch        *b_jphi_1;   //!
   TBranch        *b_jpt_2;   //!
   TBranch        *b_jeta_2;   //!
   TBranch        *b_jphi_2;   //!
   TBranch        *b_aco_mu_pi;   //!
   TBranch        *b_aco_mu_rho;   //!
   TBranch        *b_aco_mu_a1_FASTMTT_MassConstraint;   //!
   TBranch        *b_alphaAngle_mu_pi_1;   //!
   TBranch        *b_alphaAngle_mu_pi_2;   //!
   TBranch        *b_alphaAngle_mu_rho_1;   //!
   TBranch        *b_alphaAngle_mu_rho_2;   //!
   TBranch        *b_alphaAngle_mu_a1_FASTMTT_MassConstraint_1;   //!
   TBranch        *b_alphaAngle_mu_a1_FASTMTT_MassConstraint_2;   //!
   TBranch        *b_Ostar_aco_mu_pi;   //!
   TBranch        *b_Ostar_aco_mu_rho;   //!
   TBranch        *b_Ostar_aco_mu_a1_FASTMTT_MassConstraint;   //!
   TBranch        *b_PV_x;   //!
   TBranch        *b_PV_y;   //!
   TBranch        *b_PV_z;   //!
   TBranch        *b_PVBS_x;   //!
   TBranch        *b_PVBS_y;   //!
   TBranch        *b_PVBS_z;   //!
   TBranch        *b_ip_x_1;   //!
   TBranch        *b_ip_y_1;   //!
   TBranch        *b_ip_z_1;   //!
   TBranch        *b_ip_x_2;   //!
   TBranch        *b_ip_y_2;   //!
   TBranch        *b_ip_z_2;   //!
   TBranch        *b_ip_LengthSig_1;   //!
   TBranch        *b_ip_LengthSig_2;   //!
   TBranch        *b_hasRefitSV_1;   //!
   TBranch        *b_hasRefitSV_2;   //!
   TBranch        *b_sv_x_1;   //!
   TBranch        *b_sv_y_1;   //!
   TBranch        *b_sv_z_1;   //!
   TBranch        *b_sv_x_2;   //!
   TBranch        *b_sv_y_2;   //!
   TBranch        *b_sv_z_2;   //!
   TBranch        *b_FastMTT_mass;   //!
   TBranch        *b_FastMTT_pt;   //!
   TBranch        *b_FastMTT_pt_1;   //!
   TBranch        *b_FastMTT_pt_2;   //!
   TBranch        *b_FastMTT_mass_constraint;   //!
   TBranch        *b_FastMTT_pt_constraint;   //!
   TBranch        *b_FastMTT_pt_1_constraint;   //!
   TBranch        *b_FastMTT_pt_2_constraint;   //!
   TBranch        *b_pi_pt_1;   //!
   TBranch        *b_pi_eta_1;   //!
   TBranch        *b_pi_phi_1;   //!
   TBranch        *b_pi_mass_1;   //!
   TBranch        *b_pi_charge_1;   //!
   TBranch        *b_pi_pdgId_1;   //!
   TBranch        *b_pi_Energy_1;   //!
   TBranch        *b_pi2_pt_1;   //!
   TBranch        *b_pi2_eta_1;   //!
   TBranch        *b_pi2_phi_1;   //!
   TBranch        *b_pi2_mass_1;   //!
   TBranch        *b_pi2_charge_1;   //!
   TBranch        *b_pi2_pdgId_1;   //!
   TBranch        *b_pi2_Energy_1;   //!
   TBranch        *b_pi3_pt_1;   //!
   TBranch        *b_pi3_eta_1;   //!
   TBranch        *b_pi3_phi_1;   //!
   TBranch        *b_pi3_mass_1;   //!
   TBranch        *b_pi3_charge_1;   //!
   TBranch        *b_pi3_pdgId_1;   //!
   TBranch        *b_pi3_Energy_1;   //!
   TBranch        *b_pi0_pt_1;   //!
   TBranch        *b_pi0_eta_1;   //!
   TBranch        *b_pi0_phi_1;   //!
   TBranch        *b_pi0_mass_1;   //!
   TBranch        *b_pi0_charge_1;   //!
   TBranch        *b_pi0_pdgId_1;   //!
   TBranch        *b_pi0_Energy_1;   //!
   TBranch        *b_pi_pt_2;   //!
   TBranch        *b_pi_eta_2;   //!
   TBranch        *b_pi_phi_2;   //!
   TBranch        *b_pi_mass_2;   //!
   TBranch        *b_pi_charge_2;   //!
   TBranch        *b_pi_pdgId_2;   //!
   TBranch        *b_pi_Energy_2;   //!
   TBranch        *b_pi2_pt_2;   //!
   TBranch        *b_pi2_eta_2;   //!
   TBranch        *b_pi2_phi_2;   //!
   TBranch        *b_pi2_mass_2;   //!
   TBranch        *b_pi2_charge_2;   //!
   TBranch        *b_pi2_pdgId_2;   //!
   TBranch        *b_pi2_Energy_2;   //!
   TBranch        *b_pi3_pt_2;   //!
   TBranch        *b_pi3_eta_2;   //!
   TBranch        *b_pi3_phi_2;   //!
   TBranch        *b_pi3_mass_2;   //!
   TBranch        *b_pi3_charge_2;   //!
   TBranch        *b_pi3_pdgId_2;   //!
   TBranch        *b_pi3_Energy_2;   //!
   TBranch        *b_pi0_pt_2;   //!
   TBranch        *b_pi0_eta_2;   //!
   TBranch        *b_pi0_phi_2;   //!
   TBranch        *b_pi0_mass_2;   //!
   TBranch        *b_pi0_charge_2;   //!
   TBranch        *b_pi0_pdgId_2;   //!
   TBranch        *b_pi0_Energy_2;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_is_data;   //!
   TBranch        *b_pion_E_split_1;   //!
   TBranch        *b_pion_E_split_2;   //!
   TBranch        *b_iso_1;   //!
   TBranch        *b_abs_eta_1;   //!
   TBranch        *b_BDT_pred_score;   //!
   TBranch        *b_BDT_pred_class;   //!
   TBranch        *b_BDT_raw_score_tau;   //!
   TBranch        *b_BDT_raw_score_higgs;   //!
   TBranch        *b_BDT_raw_score_fake;   //!
   TBranch        *b_BDT_W_score;   //!

   ntuple(TTree *tree=0);
   virtual ~ntuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntuple_cxx
ntuple::ntuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms.cern.ch//eos/cms/store/group/phys_tau/lrussell/forAliaksei/ForFakeFactors/Run3_2022EE/mt/Muon_Run2022E/nominal/merged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms.cern.ch//eos/cms/store/group/phys_tau/lrussell/forAliaksei/ForFakeFactors/Run3_2022EE/mt/Muon_Run2022E/nominal/merged.root");
      }
      f->GetObject("ntuple",tree);

   }
   Init(tree);
}

ntuple::~ntuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("charge_1", &charge_1, &b_charge_1);
   fChain->SetBranchAddress("charge_2", &charge_2, &b_charge_2);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("mass_1", &mass_1, &b_mass_1);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("mass_2", &mass_2, &b_mass_2);
   fChain->SetBranchAddress("os", &os, &b_os);
   fChain->SetBranchAddress("dR", &dR, &b_dR);
   fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("pt_vis", &pt_vis, &b_pt_vis);
   fChain->SetBranchAddress("phi_vis", &phi_vis, &b_phi_vis);
   fChain->SetBranchAddress("eta_vis", &eta_vis, &b_eta_vis);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("mt_lep", &mt_lep, &b_mt_lep);
   fChain->SetBranchAddress("mt_tot", &mt_tot, &b_mt_tot);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("met_pt", &met_pt, &b_met_pt);
   fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
   fChain->SetBranchAddress("met_covXX", &met_covXX, &b_met_covXX);
   fChain->SetBranchAddress("met_covXY", &met_covXY, &b_met_covXY);
   fChain->SetBranchAddress("met_covYY", &met_covYY, &b_met_covYY);
   fChain->SetBranchAddress("met_dphi_1", &met_dphi_1, &b_met_dphi_1);
   fChain->SetBranchAddress("met_dphi_2", &met_dphi_2, &b_met_dphi_2);
   fChain->SetBranchAddress("trg_singlemuon", &trg_singlemuon, &b_trg_singlemuon);
   fChain->SetBranchAddress("trg_mt_cross", &trg_mt_cross, &b_trg_mt_cross);
   fChain->SetBranchAddress("idDeepTau2018v2p5VSjet_2", &idDeepTau2018v2p5VSjet_2, &b_idDeepTau2018v2p5VSjet_2);
   fChain->SetBranchAddress("idDeepTau2018v2p5VSmu_2", &idDeepTau2018v2p5VSmu_2, &b_idDeepTau2018v2p5VSmu_2);
   fChain->SetBranchAddress("idDeepTau2018v2p5VSe_2", &idDeepTau2018v2p5VSe_2, &b_idDeepTau2018v2p5VSe_2);
   fChain->SetBranchAddress("decayMode_2", &decayMode_2, &b_decayMode_2);
   fChain->SetBranchAddress("decayModePNet_2", &decayModePNet_2, &b_decayModePNet_2);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("n_bjets", &n_bjets, &b_n_bjets);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("jdeta", &jdeta, &b_jdeta);
   fChain->SetBranchAddress("sjdphi", &sjdphi, &b_sjdphi);
   fChain->SetBranchAddress("dijetpt", &dijetpt, &b_dijetpt);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   fChain->SetBranchAddress("aco_mu_pi", &aco_mu_pi, &b_aco_mu_pi);
   fChain->SetBranchAddress("aco_mu_rho", &aco_mu_rho, &b_aco_mu_rho);
   fChain->SetBranchAddress("aco_mu_a1_FASTMTT_MassConstraint", &aco_mu_a1_FASTMTT_MassConstraint, &b_aco_mu_a1_FASTMTT_MassConstraint);
   fChain->SetBranchAddress("alphaAngle_mu_pi_1", &alphaAngle_mu_pi_1, &b_alphaAngle_mu_pi_1);
   fChain->SetBranchAddress("alphaAngle_mu_pi_2", &alphaAngle_mu_pi_2, &b_alphaAngle_mu_pi_2);
   fChain->SetBranchAddress("alphaAngle_mu_rho_1", &alphaAngle_mu_rho_1, &b_alphaAngle_mu_rho_1);
   fChain->SetBranchAddress("alphaAngle_mu_rho_2", &alphaAngle_mu_rho_2, &b_alphaAngle_mu_rho_2);
   fChain->SetBranchAddress("alphaAngle_mu_a1_FASTMTT_MassConstraint_1", &alphaAngle_mu_a1_FASTMTT_MassConstraint_1, &b_alphaAngle_mu_a1_FASTMTT_MassConstraint_1);
   fChain->SetBranchAddress("alphaAngle_mu_a1_FASTMTT_MassConstraint_2", &alphaAngle_mu_a1_FASTMTT_MassConstraint_2, &b_alphaAngle_mu_a1_FASTMTT_MassConstraint_2);
   fChain->SetBranchAddress("Ostar_aco_mu_pi", &Ostar_aco_mu_pi, &b_Ostar_aco_mu_pi);
   fChain->SetBranchAddress("Ostar_aco_mu_rho", &Ostar_aco_mu_rho, &b_Ostar_aco_mu_rho);
   fChain->SetBranchAddress("Ostar_aco_mu_a1_FASTMTT_MassConstraint", &Ostar_aco_mu_a1_FASTMTT_MassConstraint, &b_Ostar_aco_mu_a1_FASTMTT_MassConstraint);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PVBS_x", &PVBS_x, &b_PVBS_x);
   fChain->SetBranchAddress("PVBS_y", &PVBS_y, &b_PVBS_y);
   fChain->SetBranchAddress("PVBS_z", &PVBS_z, &b_PVBS_z);
   fChain->SetBranchAddress("ip_x_1", &ip_x_1, &b_ip_x_1);
   fChain->SetBranchAddress("ip_y_1", &ip_y_1, &b_ip_y_1);
   fChain->SetBranchAddress("ip_z_1", &ip_z_1, &b_ip_z_1);
   fChain->SetBranchAddress("ip_x_2", &ip_x_2, &b_ip_x_2);
   fChain->SetBranchAddress("ip_y_2", &ip_y_2, &b_ip_y_2);
   fChain->SetBranchAddress("ip_z_2", &ip_z_2, &b_ip_z_2);
   fChain->SetBranchAddress("ip_LengthSig_1", &ip_LengthSig_1, &b_ip_LengthSig_1);
   fChain->SetBranchAddress("ip_LengthSig_2", &ip_LengthSig_2, &b_ip_LengthSig_2);
   fChain->SetBranchAddress("hasRefitSV_1", &hasRefitSV_1, &b_hasRefitSV_1);
   fChain->SetBranchAddress("hasRefitSV_2", &hasRefitSV_2, &b_hasRefitSV_2);
   fChain->SetBranchAddress("sv_x_1", &sv_x_1, &b_sv_x_1);
   fChain->SetBranchAddress("sv_y_1", &sv_y_1, &b_sv_y_1);
   fChain->SetBranchAddress("sv_z_1", &sv_z_1, &b_sv_z_1);
   fChain->SetBranchAddress("sv_x_2", &sv_x_2, &b_sv_x_2);
   fChain->SetBranchAddress("sv_y_2", &sv_y_2, &b_sv_y_2);
   fChain->SetBranchAddress("sv_z_2", &sv_z_2, &b_sv_z_2);
   fChain->SetBranchAddress("FastMTT_mass", &FastMTT_mass, &b_FastMTT_mass);
   fChain->SetBranchAddress("FastMTT_pt", &FastMTT_pt, &b_FastMTT_pt);
   fChain->SetBranchAddress("FastMTT_pt_1", &FastMTT_pt_1, &b_FastMTT_pt_1);
   fChain->SetBranchAddress("FastMTT_pt_2", &FastMTT_pt_2, &b_FastMTT_pt_2);
   fChain->SetBranchAddress("FastMTT_mass_constraint", &FastMTT_mass_constraint, &b_FastMTT_mass_constraint);
   fChain->SetBranchAddress("FastMTT_pt_constraint", &FastMTT_pt_constraint, &b_FastMTT_pt_constraint);
   fChain->SetBranchAddress("FastMTT_pt_1_constraint", &FastMTT_pt_1_constraint, &b_FastMTT_pt_1_constraint);
   fChain->SetBranchAddress("FastMTT_pt_2_constraint", &FastMTT_pt_2_constraint, &b_FastMTT_pt_2_constraint);
   fChain->SetBranchAddress("pi_pt_1", &pi_pt_1, &b_pi_pt_1);
   fChain->SetBranchAddress("pi_eta_1", &pi_eta_1, &b_pi_eta_1);
   fChain->SetBranchAddress("pi_phi_1", &pi_phi_1, &b_pi_phi_1);
   fChain->SetBranchAddress("pi_mass_1", &pi_mass_1, &b_pi_mass_1);
   fChain->SetBranchAddress("pi_charge_1", &pi_charge_1, &b_pi_charge_1);
   fChain->SetBranchAddress("pi_pdgId_1", &pi_pdgId_1, &b_pi_pdgId_1);
   fChain->SetBranchAddress("pi_Energy_1", &pi_Energy_1, &b_pi_Energy_1);
   fChain->SetBranchAddress("pi2_pt_1", &pi2_pt_1, &b_pi2_pt_1);
   fChain->SetBranchAddress("pi2_eta_1", &pi2_eta_1, &b_pi2_eta_1);
   fChain->SetBranchAddress("pi2_phi_1", &pi2_phi_1, &b_pi2_phi_1);
   fChain->SetBranchAddress("pi2_mass_1", &pi2_mass_1, &b_pi2_mass_1);
   fChain->SetBranchAddress("pi2_charge_1", &pi2_charge_1, &b_pi2_charge_1);
   fChain->SetBranchAddress("pi2_pdgId_1", &pi2_pdgId_1, &b_pi2_pdgId_1);
   fChain->SetBranchAddress("pi2_Energy_1", &pi2_Energy_1, &b_pi2_Energy_1);
   fChain->SetBranchAddress("pi3_pt_1", &pi3_pt_1, &b_pi3_pt_1);
   fChain->SetBranchAddress("pi3_eta_1", &pi3_eta_1, &b_pi3_eta_1);
   fChain->SetBranchAddress("pi3_phi_1", &pi3_phi_1, &b_pi3_phi_1);
   fChain->SetBranchAddress("pi3_mass_1", &pi3_mass_1, &b_pi3_mass_1);
   fChain->SetBranchAddress("pi3_charge_1", &pi3_charge_1, &b_pi3_charge_1);
   fChain->SetBranchAddress("pi3_pdgId_1", &pi3_pdgId_1, &b_pi3_pdgId_1);
   fChain->SetBranchAddress("pi3_Energy_1", &pi3_Energy_1, &b_pi3_Energy_1);
   fChain->SetBranchAddress("pi0_pt_1", &pi0_pt_1, &b_pi0_pt_1);
   fChain->SetBranchAddress("pi0_eta_1", &pi0_eta_1, &b_pi0_eta_1);
   fChain->SetBranchAddress("pi0_phi_1", &pi0_phi_1, &b_pi0_phi_1);
   fChain->SetBranchAddress("pi0_mass_1", &pi0_mass_1, &b_pi0_mass_1);
   fChain->SetBranchAddress("pi0_charge_1", &pi0_charge_1, &b_pi0_charge_1);
   fChain->SetBranchAddress("pi0_pdgId_1", &pi0_pdgId_1, &b_pi0_pdgId_1);
   fChain->SetBranchAddress("pi0_Energy_1", &pi0_Energy_1, &b_pi0_Energy_1);
   fChain->SetBranchAddress("pi_pt_2", &pi_pt_2, &b_pi_pt_2);
   fChain->SetBranchAddress("pi_eta_2", &pi_eta_2, &b_pi_eta_2);
   fChain->SetBranchAddress("pi_phi_2", &pi_phi_2, &b_pi_phi_2);
   fChain->SetBranchAddress("pi_mass_2", &pi_mass_2, &b_pi_mass_2);
   fChain->SetBranchAddress("pi_charge_2", &pi_charge_2, &b_pi_charge_2);
   fChain->SetBranchAddress("pi_pdgId_2", &pi_pdgId_2, &b_pi_pdgId_2);
   fChain->SetBranchAddress("pi_Energy_2", &pi_Energy_2, &b_pi_Energy_2);
   fChain->SetBranchAddress("pi2_pt_2", &pi2_pt_2, &b_pi2_pt_2);
   fChain->SetBranchAddress("pi2_eta_2", &pi2_eta_2, &b_pi2_eta_2);
   fChain->SetBranchAddress("pi2_phi_2", &pi2_phi_2, &b_pi2_phi_2);
   fChain->SetBranchAddress("pi2_mass_2", &pi2_mass_2, &b_pi2_mass_2);
   fChain->SetBranchAddress("pi2_charge_2", &pi2_charge_2, &b_pi2_charge_2);
   fChain->SetBranchAddress("pi2_pdgId_2", &pi2_pdgId_2, &b_pi2_pdgId_2);
   fChain->SetBranchAddress("pi2_Energy_2", &pi2_Energy_2, &b_pi2_Energy_2);
   fChain->SetBranchAddress("pi3_pt_2", &pi3_pt_2, &b_pi3_pt_2);
   fChain->SetBranchAddress("pi3_eta_2", &pi3_eta_2, &b_pi3_eta_2);
   fChain->SetBranchAddress("pi3_phi_2", &pi3_phi_2, &b_pi3_phi_2);
   fChain->SetBranchAddress("pi3_mass_2", &pi3_mass_2, &b_pi3_mass_2);
   fChain->SetBranchAddress("pi3_charge_2", &pi3_charge_2, &b_pi3_charge_2);
   fChain->SetBranchAddress("pi3_pdgId_2", &pi3_pdgId_2, &b_pi3_pdgId_2);
   fChain->SetBranchAddress("pi3_Energy_2", &pi3_Energy_2, &b_pi3_Energy_2);
   fChain->SetBranchAddress("pi0_pt_2", &pi0_pt_2, &b_pi0_pt_2);
   fChain->SetBranchAddress("pi0_eta_2", &pi0_eta_2, &b_pi0_eta_2);
   fChain->SetBranchAddress("pi0_phi_2", &pi0_phi_2, &b_pi0_phi_2);
   fChain->SetBranchAddress("pi0_mass_2", &pi0_mass_2, &b_pi0_mass_2);
   fChain->SetBranchAddress("pi0_charge_2", &pi0_charge_2, &b_pi0_charge_2);
   fChain->SetBranchAddress("pi0_pdgId_2", &pi0_pdgId_2, &b_pi0_pdgId_2);
   fChain->SetBranchAddress("pi0_Energy_2", &pi0_Energy_2, &b_pi0_Energy_2);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("is_data", &is_data, &b_is_data);
   fChain->SetBranchAddress("pion_E_split_1", &pion_E_split_1, &b_pion_E_split_1);
   fChain->SetBranchAddress("pion_E_split_2", &pion_E_split_2, &b_pion_E_split_2);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("abs_eta_1", &abs_eta_1, &b_abs_eta_1);
   fChain->SetBranchAddress("BDT_pred_score", &BDT_pred_score, &b_BDT_pred_score);
   fChain->SetBranchAddress("BDT_pred_class", &BDT_pred_class, &b_BDT_pred_class);
   fChain->SetBranchAddress("BDT_raw_score_tau", &BDT_raw_score_tau, &b_BDT_raw_score_tau);
   fChain->SetBranchAddress("BDT_raw_score_higgs", &BDT_raw_score_higgs, &b_BDT_raw_score_higgs);
   fChain->SetBranchAddress("BDT_raw_score_fake", &BDT_raw_score_fake, &b_BDT_raw_score_fake);
   fChain->SetBranchAddress("BDT_W_score", &BDT_W_score, &b_BDT_W_score);
   Notify();
}

Bool_t ntuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntuple_cxx
