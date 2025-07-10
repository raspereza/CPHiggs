//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  7 10:59:37 2025 by ROOT version 6.30/07
// from TTree ntuple/
// found on file: root://eoscms.cern.ch//eos/cms/store/group/phys_tau/lrussell/forAliaksei/Run3_2022/mt/GluGluHTo2Tau_UncorrelatedDecay_Filtered/nominal/merged.root
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
   Double_t        original_index_1;
   Double_t        original_index_2;
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
   Long64_t        idDeepTau2018v2p5noDAVSjet_2;
   Long64_t        idDeepTau2018v2p5noDAVSmu_2;
   Long64_t        idDeepTau2018v2p5noDAVSe_2;
   Double_t        rawDeepTau2018v2p5VSjet_2;
   Double_t        rawDeepTau2018v2p5VSmu_2;
   Double_t        rawDeepTau2018v2p5VSe_2;
   Double_t        rawDeepTau2018v2p5noDAVSjet_2;
   Double_t        rawDeepTau2018v2p5noDAVSmu_2;
   Double_t        rawDeepTau2018v2p5noDAVSe_2;
   Double_t        rawPNetVSjet_2;
   Double_t        rawPNetVSmu_2;
   Double_t        rawPNetVSe_2;
   Long64_t        decayMode_2;
   Short_t         decayModePNet_2;
   Double_t        probDM0PNet_2;
   Double_t        probDM1PNet_2;
   Double_t        probDM2PNet_2;
   Double_t        probDM10PNet_2;
   Double_t        probDM11PNet_2;
   Double_t        n_jets;
   Double_t        n_prebjets;
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
   Double_t        seeding_n_jets;
   Double_t        seeding_mjj;
   Double_t        seeding_jdeta;
   Double_t        seeding_sjdphi;
   Double_t        seeding_dijetpt;
   Double_t        seeding_jpt_1;
   Double_t        seeding_jeta_1;
   Double_t        seeding_jphi_1;
   Double_t        seeding_jpt_2;
   Double_t        seeding_jeta_2;
   Double_t        seeding_jphi_2;
   Double_t        aco_mu_pi;
   Double_t        aco_mu_rho;
   Double_t        aco_mu_a1;
   Double_t        alphaAngle_mu_pi_1;
   Double_t        alphaAngle_mu_pi_2;
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
   Double_t        PVBS_cov00;
   Double_t        PVBS_cov10;
   Double_t        PVBS_cov11;
   Double_t        PVBS_cov20;
   Double_t        PVBS_cov21;
   Double_t        PVBS_cov22;
   Double_t        sv_cov00_1;
   Double_t        sv_cov10_1;
   Double_t        sv_cov11_1;
   Double_t        sv_cov20_1;
   Double_t        sv_cov21_1;
   Double_t        sv_cov22_1;
   Double_t        sv_cov00_2;
   Double_t        sv_cov10_2;
   Double_t        sv_cov11_2;
   Double_t        sv_cov20_2;
   Double_t        sv_cov21_2;
   Double_t        sv_cov22_2;
   Double_t        ip_cov00_1;
   Double_t        ip_cov10_1;
   Double_t        ip_cov11_1;
   Double_t        ip_cov20_1;
   Double_t        ip_cov21_1;
   Double_t        ip_cov22_1;
   Double_t        ip_cov00_2;
   Double_t        ip_cov10_2;
   Double_t        ip_cov11_2;
   Double_t        ip_cov20_2;
   Double_t        ip_cov21_2;
   Double_t        ip_cov22_2;
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
   Double_t        leadTkPtOverTauPt_2;
   Double_t        weight;
   Double_t        genWeight;
   Double_t        LHEReweightingWeight_SM;
   Double_t        LHEReweightingWeight_PS;
   Double_t        LHEReweightingWeight_MM;
   Long64_t        nLHEjets;
   Long64_t        npNLOjets;
   Double_t        LHE_Vpt;
   Double_t        wt_cp_sm;
   Double_t        wt_cp_ps;
   Double_t        wt_cp_mm;
   Double_t        pion_E_split_1;
   Double_t        pion_E_split_2;
   Double_t        gen_boson_pT;
   Double_t        gen_boson_mass;
   Double_t        gen_boson_eta;
   Double_t        gen_boson_phi;
   Double_t        gen_taunus_pT;
   Double_t        gen_taunus_phi;
   Long64_t        genPartFlav_1;
   Long64_t        genPartFlav_2;
   Double_t        genPart_pt_1;
   Double_t        genPart_eta_1;
   Double_t        genPart_phi_1;
   Double_t        genPart_pdgId_1;
   Double_t        genPart_pt_2;
   Double_t        genPart_eta_2;
   Double_t        genPart_phi_2;
   Double_t        genPart_pdgId_2;
   Double_t        genVisTau_pt_1;
   Double_t        genVisTau_eta_1;
   Double_t        genVisTau_phi_1;
   Double_t        genVisTau_mass_1;
   Double_t        genVisTau_pt_2;
   Double_t        genVisTau_eta_2;
   Double_t        genVisTau_phi_2;
   Double_t        genVisTau_mass_2;
   Double_t        gen_decayMode_1;
   Double_t        gen_decayMode_2;
   Double_t        genIP_1_x;
   Double_t        genIP_1_y;
   Double_t        genIP_1_z;
   Double_t        genIP_2_x;
   Double_t        genIP_2_y;
   Double_t        genIP_2_z;
   Double_t        GenVsReco_PVBS_dxy;
   Double_t        GenVsReco_PVBS_dz;
   Double_t        GenVsReco_PV_dxy;
   Double_t        GenVsReco_PV_dz;
   Double_t        w_DY_soup;
   Double_t        w_WJ_soup;
   Double_t        w_DY_NLO_soup;
   Double_t        w_Pileup;
   Double_t        w_Top_pt_Reweighting;
   Double_t        w_Zpt_Reweighting_Imperial;
   Double_t        w_ggH_QuarkMass_Effects;
   Double_t        w_Electron_ID;
   Double_t        w_Electron_Reco;
   Double_t        w_Muon_ID;
   Double_t        w_Muon_Isolation;
   Double_t        w_Tau_ID_PNet;
   Double_t        w_Trigger;
   Double_t        w_SingleMuon_leg1;
   Double_t        w_CrossTrigger_Tau_leg2;
   Double_t        w_CrossTrigger_Muon_leg1;
   Double_t        iso_1;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_original_index_1;   //!
   TBranch        *b_original_index_2;   //!
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
   TBranch        *b_idDeepTau2018v2p5noDAVSjet_2;   //!
   TBranch        *b_idDeepTau2018v2p5noDAVSmu_2;   //!
   TBranch        *b_idDeepTau2018v2p5noDAVSe_2;   //!
   TBranch        *b_rawDeepTau2018v2p5VSjet_2;   //!
   TBranch        *b_rawDeepTau2018v2p5VSmu_2;   //!
   TBranch        *b_rawDeepTau2018v2p5VSe_2;   //!
   TBranch        *b_rawDeepTau2018v2p5noDAVSjet_2;   //!
   TBranch        *b_rawDeepTau2018v2p5noDAVSmu_2;   //!
   TBranch        *b_rawDeepTau2018v2p5noDAVSe_2;   //!
   TBranch        *b_rawPNetVSjet_2;   //!
   TBranch        *b_rawPNetVSmu_2;   //!
   TBranch        *b_rawPNetVSe_2;   //!
   TBranch        *b_decayMode_2;   //!
   TBranch        *b_decayModePNet_2;   //!
   TBranch        *b_probDM0PNet_2;   //!
   TBranch        *b_probDM1PNet_2;   //!
   TBranch        *b_probDM2PNet_2;   //!
   TBranch        *b_probDM10PNet_2;   //!
   TBranch        *b_probDM11PNet_2;   //!
   TBranch        *b_n_jets;   //!
   TBranch        *b_n_prebjets;   //!
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
   TBranch        *b_seeding_n_jets;   //!
   TBranch        *b_seeding_mjj;   //!
   TBranch        *b_seeding_jdeta;   //!
   TBranch        *b_seeding_sjdphi;   //!
   TBranch        *b_seeding_dijetpt;   //!
   TBranch        *b_seeding_jpt_1;   //!
   TBranch        *b_seeding_jeta_1;   //!
   TBranch        *b_seeding_jphi_1;   //!
   TBranch        *b_seeding_jpt_2;   //!
   TBranch        *b_seeding_jeta_2;   //!
   TBranch        *b_seeding_jphi_2;   //!
   TBranch        *b_aco_mu_pi;   //!
   TBranch        *b_aco_mu_rho;   //!
   TBranch        *b_aco_mu_a1;   //!
   TBranch        *b_alphaAngle_mu_pi_1;   //!
   TBranch        *b_alphaAngle_mu_pi_2;   //!
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
   TBranch        *b_PVBS_cov00;   //!
   TBranch        *b_PVBS_cov10;   //!
   TBranch        *b_PVBS_cov11;   //!
   TBranch        *b_PVBS_cov20;   //!
   TBranch        *b_PVBS_cov21;   //!
   TBranch        *b_PVBS_cov22;   //!
   TBranch        *b_sv_cov00_1;   //!
   TBranch        *b_sv_cov10_1;   //!
   TBranch        *b_sv_cov11_1;   //!
   TBranch        *b_sv_cov20_1;   //!
   TBranch        *b_sv_cov21_1;   //!
   TBranch        *b_sv_cov22_1;   //!
   TBranch        *b_sv_cov00_2;   //!
   TBranch        *b_sv_cov10_2;   //!
   TBranch        *b_sv_cov11_2;   //!
   TBranch        *b_sv_cov20_2;   //!
   TBranch        *b_sv_cov21_2;   //!
   TBranch        *b_sv_cov22_2;   //!
   TBranch        *b_ip_cov00_1;   //!
   TBranch        *b_ip_cov10_1;   //!
   TBranch        *b_ip_cov11_1;   //!
   TBranch        *b_ip_cov20_1;   //!
   TBranch        *b_ip_cov21_1;   //!
   TBranch        *b_ip_cov22_1;   //!
   TBranch        *b_ip_cov00_2;   //!
   TBranch        *b_ip_cov10_2;   //!
   TBranch        *b_ip_cov11_2;   //!
   TBranch        *b_ip_cov20_2;   //!
   TBranch        *b_ip_cov21_2;   //!
   TBranch        *b_ip_cov22_2;   //!
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
   TBranch        *b_leadTkPtOverTauPt_2;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_LHEReweightingWeight_SM;   //!
   TBranch        *b_LHEReweightingWeight_PS;   //!
   TBranch        *b_LHEReweightingWeight_MM;   //!
   TBranch        *b_nLHEjets;   //!
   TBranch        *b_npNLOjets;   //!
   TBranch        *b_LHE_Vpt;   //!
   TBranch        *b_wt_cp_sm;   //!
   TBranch        *b_wt_cp_ps;   //!
   TBranch        *b_wt_cp_mm;   //!
   TBranch        *b_pion_E_split_1;   //!
   TBranch        *b_pion_E_split_2;   //!
   TBranch        *b_gen_boson_pT;   //!
   TBranch        *b_gen_boson_mass;   //!
   TBranch        *b_gen_boson_eta;   //!
   TBranch        *b_gen_boson_phi;   //!
   TBranch        *b_gen_taunus_pT;   //!
   TBranch        *b_gen_taunus_phi;   //!
   TBranch        *b_genPartFlav_1;   //!
   TBranch        *b_genPartFlav_2;   //!
   TBranch        *b_genPart_pt_1;   //!
   TBranch        *b_genPart_eta_1;   //!
   TBranch        *b_genPart_phi_1;   //!
   TBranch        *b_genPart_pdgId_1;   //!
   TBranch        *b_genPart_pt_2;   //!
   TBranch        *b_genPart_eta_2;   //!
   TBranch        *b_genPart_phi_2;   //!
   TBranch        *b_genPart_pdgId_2;   //!
   TBranch        *b_genVisTau_pt_1;   //!
   TBranch        *b_genVisTau_eta_1;   //!
   TBranch        *b_genVisTau_phi_1;   //!
   TBranch        *b_genVisTau_mass_1;   //!
   TBranch        *b_genVisTau_pt_2;   //!
   TBranch        *b_genVisTau_eta_2;   //!
   TBranch        *b_genVisTau_phi_2;   //!
   TBranch        *b_genVisTau_mass_2;   //!
   TBranch        *b_gen_decayMode_1;   //!
   TBranch        *b_gen_decayMode_2;   //!
   TBranch        *b_genIP_1_x;   //!
   TBranch        *b_genIP_1_y;   //!
   TBranch        *b_genIP_1_z;   //!
   TBranch        *b_genIP_2_x;   //!
   TBranch        *b_genIP_2_y;   //!
   TBranch        *b_genIP_2_z;   //!
   TBranch        *b_GenVsReco_PVBS_dxy;   //!
   TBranch        *b_GenVsReco_PVBS_dz;   //!
   TBranch        *b_GenVsReco_PV_dxy;   //!
   TBranch        *b_GenVsReco_PV_dz;   //!
   TBranch        *b_w_DY_soup;   //!
   TBranch        *b_w_WJ_soup;   //!
   TBranch        *b_w_DY_NLO_soup;   //!
   TBranch        *b_w_Pileup;   //!
   TBranch        *b_w_Top_pt_Reweighting;   //!
   TBranch        *b_w_Zpt_Reweighting_Imperial;   //!
   TBranch        *b_w_ggH_QuarkMass_Effects;   //!
   TBranch        *b_w_Electron_ID;   //!
   TBranch        *b_w_Electron_Reco;   //!
   TBranch        *b_w_Muon_ID;   //!
   TBranch        *b_w_Muon_Isolation;   //!
   TBranch        *b_w_Tau_ID_PNet;   //!
   TBranch        *b_w_Trigger;   //!
   TBranch        *b_w_SingleMuon_leg1;   //!
   TBranch        *b_w_CrossTrigger_Tau_leg2;   //!
   TBranch        *b_w_CrossTrigger_Muon_leg1;   //!
   TBranch        *b_iso_1;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoscms.cern.ch//eos/cms/store/group/phys_tau/lrussell/forAliaksei/Run3_2022/mt/GluGluHTo2Tau_UncorrelatedDecay_Filtered/nominal/merged.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoscms.cern.ch//eos/cms/store/group/phys_tau/lrussell/forAliaksei/Run3_2022/mt/GluGluHTo2Tau_UncorrelatedDecay_Filtered/nominal/merged.root");
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
   fChain->SetBranchAddress("original_index_1", &original_index_1, &b_original_index_1);
   fChain->SetBranchAddress("original_index_2", &original_index_2, &b_original_index_2);
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
   fChain->SetBranchAddress("idDeepTau2018v2p5noDAVSjet_2", &idDeepTau2018v2p5noDAVSjet_2, &b_idDeepTau2018v2p5noDAVSjet_2);
   fChain->SetBranchAddress("idDeepTau2018v2p5noDAVSmu_2", &idDeepTau2018v2p5noDAVSmu_2, &b_idDeepTau2018v2p5noDAVSmu_2);
   fChain->SetBranchAddress("idDeepTau2018v2p5noDAVSe_2", &idDeepTau2018v2p5noDAVSe_2, &b_idDeepTau2018v2p5noDAVSe_2);
   fChain->SetBranchAddress("rawDeepTau2018v2p5VSjet_2", &rawDeepTau2018v2p5VSjet_2, &b_rawDeepTau2018v2p5VSjet_2);
   fChain->SetBranchAddress("rawDeepTau2018v2p5VSmu_2", &rawDeepTau2018v2p5VSmu_2, &b_rawDeepTau2018v2p5VSmu_2);
   fChain->SetBranchAddress("rawDeepTau2018v2p5VSe_2", &rawDeepTau2018v2p5VSe_2, &b_rawDeepTau2018v2p5VSe_2);
   fChain->SetBranchAddress("rawDeepTau2018v2p5noDAVSjet_2", &rawDeepTau2018v2p5noDAVSjet_2, &b_rawDeepTau2018v2p5noDAVSjet_2);
   fChain->SetBranchAddress("rawDeepTau2018v2p5noDAVSmu_2", &rawDeepTau2018v2p5noDAVSmu_2, &b_rawDeepTau2018v2p5noDAVSmu_2);
   fChain->SetBranchAddress("rawDeepTau2018v2p5noDAVSe_2", &rawDeepTau2018v2p5noDAVSe_2, &b_rawDeepTau2018v2p5noDAVSe_2);
   fChain->SetBranchAddress("rawPNetVSjet_2", &rawPNetVSjet_2, &b_rawPNetVSjet_2);
   fChain->SetBranchAddress("rawPNetVSmu_2", &rawPNetVSmu_2, &b_rawPNetVSmu_2);
   fChain->SetBranchAddress("rawPNetVSe_2", &rawPNetVSe_2, &b_rawPNetVSe_2);
   fChain->SetBranchAddress("decayMode_2", &decayMode_2, &b_decayMode_2);
   fChain->SetBranchAddress("decayModePNet_2", &decayModePNet_2, &b_decayModePNet_2);
   fChain->SetBranchAddress("probDM0PNet_2", &probDM0PNet_2, &b_probDM0PNet_2);
   fChain->SetBranchAddress("probDM1PNet_2", &probDM1PNet_2, &b_probDM1PNet_2);
   fChain->SetBranchAddress("probDM2PNet_2", &probDM2PNet_2, &b_probDM2PNet_2);
   fChain->SetBranchAddress("probDM10PNet_2", &probDM10PNet_2, &b_probDM10PNet_2);
   fChain->SetBranchAddress("probDM11PNet_2", &probDM11PNet_2, &b_probDM11PNet_2);
   fChain->SetBranchAddress("n_jets", &n_jets, &b_n_jets);
   fChain->SetBranchAddress("n_prebjets", &n_prebjets, &b_n_prebjets);
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
   fChain->SetBranchAddress("seeding_n_jets", &seeding_n_jets, &b_seeding_n_jets);
   fChain->SetBranchAddress("seeding_mjj", &seeding_mjj, &b_seeding_mjj);
   fChain->SetBranchAddress("seeding_jdeta", &seeding_jdeta, &b_seeding_jdeta);
   fChain->SetBranchAddress("seeding_sjdphi", &seeding_sjdphi, &b_seeding_sjdphi);
   fChain->SetBranchAddress("seeding_dijetpt", &seeding_dijetpt, &b_seeding_dijetpt);
   fChain->SetBranchAddress("seeding_jpt_1", &seeding_jpt_1, &b_seeding_jpt_1);
   fChain->SetBranchAddress("seeding_jeta_1", &seeding_jeta_1, &b_seeding_jeta_1);
   fChain->SetBranchAddress("seeding_jphi_1", &seeding_jphi_1, &b_seeding_jphi_1);
   fChain->SetBranchAddress("seeding_jpt_2", &seeding_jpt_2, &b_seeding_jpt_2);
   fChain->SetBranchAddress("seeding_jeta_2", &seeding_jeta_2, &b_seeding_jeta_2);
   fChain->SetBranchAddress("seeding_jphi_2", &seeding_jphi_2, &b_seeding_jphi_2);
   fChain->SetBranchAddress("aco_mu_pi", &aco_mu_pi, &b_aco_mu_pi);
   fChain->SetBranchAddress("aco_mu_rho", &aco_mu_rho, &b_aco_mu_rho);
   fChain->SetBranchAddress("aco_mu_a1", &aco_mu_a1, &b_aco_mu_a1);
   fChain->SetBranchAddress("alphaAngle_mu_pi_1", &alphaAngle_mu_pi_1, &b_alphaAngle_mu_pi_1);
   fChain->SetBranchAddress("alphaAngle_mu_pi_2", &alphaAngle_mu_pi_2, &b_alphaAngle_mu_pi_2);
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
   fChain->SetBranchAddress("PVBS_cov00", &PVBS_cov00, &b_PVBS_cov00);
   fChain->SetBranchAddress("PVBS_cov10", &PVBS_cov10, &b_PVBS_cov10);
   fChain->SetBranchAddress("PVBS_cov11", &PVBS_cov11, &b_PVBS_cov11);
   fChain->SetBranchAddress("PVBS_cov20", &PVBS_cov20, &b_PVBS_cov20);
   fChain->SetBranchAddress("PVBS_cov21", &PVBS_cov21, &b_PVBS_cov21);
   fChain->SetBranchAddress("PVBS_cov22", &PVBS_cov22, &b_PVBS_cov22);
   fChain->SetBranchAddress("sv_cov00_1", &sv_cov00_1, &b_sv_cov00_1);
   fChain->SetBranchAddress("sv_cov10_1", &sv_cov10_1, &b_sv_cov10_1);
   fChain->SetBranchAddress("sv_cov11_1", &sv_cov11_1, &b_sv_cov11_1);
   fChain->SetBranchAddress("sv_cov20_1", &sv_cov20_1, &b_sv_cov20_1);
   fChain->SetBranchAddress("sv_cov21_1", &sv_cov21_1, &b_sv_cov21_1);
   fChain->SetBranchAddress("sv_cov22_1", &sv_cov22_1, &b_sv_cov22_1);
   fChain->SetBranchAddress("sv_cov00_2", &sv_cov00_2, &b_sv_cov00_2);
   fChain->SetBranchAddress("sv_cov10_2", &sv_cov10_2, &b_sv_cov10_2);
   fChain->SetBranchAddress("sv_cov11_2", &sv_cov11_2, &b_sv_cov11_2);
   fChain->SetBranchAddress("sv_cov20_2", &sv_cov20_2, &b_sv_cov20_2);
   fChain->SetBranchAddress("sv_cov21_2", &sv_cov21_2, &b_sv_cov21_2);
   fChain->SetBranchAddress("sv_cov22_2", &sv_cov22_2, &b_sv_cov22_2);
   fChain->SetBranchAddress("ip_cov00_1", &ip_cov00_1, &b_ip_cov00_1);
   fChain->SetBranchAddress("ip_cov10_1", &ip_cov10_1, &b_ip_cov10_1);
   fChain->SetBranchAddress("ip_cov11_1", &ip_cov11_1, &b_ip_cov11_1);
   fChain->SetBranchAddress("ip_cov20_1", &ip_cov20_1, &b_ip_cov20_1);
   fChain->SetBranchAddress("ip_cov21_1", &ip_cov21_1, &b_ip_cov21_1);
   fChain->SetBranchAddress("ip_cov22_1", &ip_cov22_1, &b_ip_cov22_1);
   fChain->SetBranchAddress("ip_cov00_2", &ip_cov00_2, &b_ip_cov00_2);
   fChain->SetBranchAddress("ip_cov10_2", &ip_cov10_2, &b_ip_cov10_2);
   fChain->SetBranchAddress("ip_cov11_2", &ip_cov11_2, &b_ip_cov11_2);
   fChain->SetBranchAddress("ip_cov20_2", &ip_cov20_2, &b_ip_cov20_2);
   fChain->SetBranchAddress("ip_cov21_2", &ip_cov21_2, &b_ip_cov21_2);
   fChain->SetBranchAddress("ip_cov22_2", &ip_cov22_2, &b_ip_cov22_2);
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
   fChain->SetBranchAddress("leadTkPtOverTauPt_2", &leadTkPtOverTauPt_2, &b_leadTkPtOverTauPt_2);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("LHEReweightingWeight_SM", &LHEReweightingWeight_SM, &b_LHEReweightingWeight_SM);
   fChain->SetBranchAddress("LHEReweightingWeight_PS", &LHEReweightingWeight_PS, &b_LHEReweightingWeight_PS);
   fChain->SetBranchAddress("LHEReweightingWeight_MM", &LHEReweightingWeight_MM, &b_LHEReweightingWeight_MM);
   fChain->SetBranchAddress("nLHEjets", &nLHEjets, &b_nLHEjets);
   fChain->SetBranchAddress("npNLOjets", &npNLOjets, &b_npNLOjets);
   fChain->SetBranchAddress("LHE_Vpt", &LHE_Vpt, &b_LHE_Vpt);
   fChain->SetBranchAddress("wt_cp_sm", &wt_cp_sm, &b_wt_cp_sm);
   fChain->SetBranchAddress("wt_cp_ps", &wt_cp_ps, &b_wt_cp_ps);
   fChain->SetBranchAddress("wt_cp_mm", &wt_cp_mm, &b_wt_cp_mm);
   fChain->SetBranchAddress("pion_E_split_1", &pion_E_split_1, &b_pion_E_split_1);
   fChain->SetBranchAddress("pion_E_split_2", &pion_E_split_2, &b_pion_E_split_2);
   fChain->SetBranchAddress("gen_boson_pT", &gen_boson_pT, &b_gen_boson_pT);
   fChain->SetBranchAddress("gen_boson_mass", &gen_boson_mass, &b_gen_boson_mass);
   fChain->SetBranchAddress("gen_boson_eta", &gen_boson_eta, &b_gen_boson_eta);
   fChain->SetBranchAddress("gen_boson_phi", &gen_boson_phi, &b_gen_boson_phi);
   fChain->SetBranchAddress("gen_taunus_pT", &gen_taunus_pT, &b_gen_taunus_pT);
   fChain->SetBranchAddress("gen_taunus_phi", &gen_taunus_phi, &b_gen_taunus_phi);
   fChain->SetBranchAddress("genPartFlav_1", &genPartFlav_1, &b_genPartFlav_1);
   fChain->SetBranchAddress("genPartFlav_2", &genPartFlav_2, &b_genPartFlav_2);
   fChain->SetBranchAddress("genPart_pt_1", &genPart_pt_1, &b_genPart_pt_1);
   fChain->SetBranchAddress("genPart_eta_1", &genPart_eta_1, &b_genPart_eta_1);
   fChain->SetBranchAddress("genPart_phi_1", &genPart_phi_1, &b_genPart_phi_1);
   fChain->SetBranchAddress("genPart_pdgId_1", &genPart_pdgId_1, &b_genPart_pdgId_1);
   fChain->SetBranchAddress("genPart_pt_2", &genPart_pt_2, &b_genPart_pt_2);
   fChain->SetBranchAddress("genPart_eta_2", &genPart_eta_2, &b_genPart_eta_2);
   fChain->SetBranchAddress("genPart_phi_2", &genPart_phi_2, &b_genPart_phi_2);
   fChain->SetBranchAddress("genPart_pdgId_2", &genPart_pdgId_2, &b_genPart_pdgId_2);
   fChain->SetBranchAddress("genVisTau_pt_1", &genVisTau_pt_1, &b_genVisTau_pt_1);
   fChain->SetBranchAddress("genVisTau_eta_1", &genVisTau_eta_1, &b_genVisTau_eta_1);
   fChain->SetBranchAddress("genVisTau_phi_1", &genVisTau_phi_1, &b_genVisTau_phi_1);
   fChain->SetBranchAddress("genVisTau_mass_1", &genVisTau_mass_1, &b_genVisTau_mass_1);
   fChain->SetBranchAddress("genVisTau_pt_2", &genVisTau_pt_2, &b_genVisTau_pt_2);
   fChain->SetBranchAddress("genVisTau_eta_2", &genVisTau_eta_2, &b_genVisTau_eta_2);
   fChain->SetBranchAddress("genVisTau_phi_2", &genVisTau_phi_2, &b_genVisTau_phi_2);
   fChain->SetBranchAddress("genVisTau_mass_2", &genVisTau_mass_2, &b_genVisTau_mass_2);
   fChain->SetBranchAddress("gen_decayMode_1", &gen_decayMode_1, &b_gen_decayMode_1);
   fChain->SetBranchAddress("gen_decayMode_2", &gen_decayMode_2, &b_gen_decayMode_2);
   fChain->SetBranchAddress("genIP_1_x", &genIP_1_x, &b_genIP_1_x);
   fChain->SetBranchAddress("genIP_1_y", &genIP_1_y, &b_genIP_1_y);
   fChain->SetBranchAddress("genIP_1_z", &genIP_1_z, &b_genIP_1_z);
   fChain->SetBranchAddress("genIP_2_x", &genIP_2_x, &b_genIP_2_x);
   fChain->SetBranchAddress("genIP_2_y", &genIP_2_y, &b_genIP_2_y);
   fChain->SetBranchAddress("genIP_2_z", &genIP_2_z, &b_genIP_2_z);
   fChain->SetBranchAddress("GenVsReco_PVBS_dxy", &GenVsReco_PVBS_dxy, &b_GenVsReco_PVBS_dxy);
   fChain->SetBranchAddress("GenVsReco_PVBS_dz", &GenVsReco_PVBS_dz, &b_GenVsReco_PVBS_dz);
   fChain->SetBranchAddress("GenVsReco_PV_dxy", &GenVsReco_PV_dxy, &b_GenVsReco_PV_dxy);
   fChain->SetBranchAddress("GenVsReco_PV_dz", &GenVsReco_PV_dz, &b_GenVsReco_PV_dz);
   fChain->SetBranchAddress("w_DY_soup", &w_DY_soup, &b_w_DY_soup);
   fChain->SetBranchAddress("w_WJ_soup", &w_WJ_soup, &b_w_WJ_soup);
   fChain->SetBranchAddress("w_DY_NLO_soup", &w_DY_NLO_soup, &b_w_DY_NLO_soup);
   fChain->SetBranchAddress("w_Pileup", &w_Pileup, &b_w_Pileup);
   fChain->SetBranchAddress("w_Top_pt_Reweighting", &w_Top_pt_Reweighting, &b_w_Top_pt_Reweighting);
   fChain->SetBranchAddress("w_Zpt_Reweighting_Imperial", &w_Zpt_Reweighting_Imperial, &b_w_Zpt_Reweighting_Imperial);
   fChain->SetBranchAddress("w_ggH_QuarkMass_Effects", &w_ggH_QuarkMass_Effects, &b_w_ggH_QuarkMass_Effects);
   fChain->SetBranchAddress("w_Electron_ID", &w_Electron_ID, &b_w_Electron_ID);
   fChain->SetBranchAddress("w_Electron_Reco", &w_Electron_Reco, &b_w_Electron_Reco);
   fChain->SetBranchAddress("w_Muon_ID", &w_Muon_ID, &b_w_Muon_ID);
   fChain->SetBranchAddress("w_Muon_Isolation", &w_Muon_Isolation, &b_w_Muon_Isolation);
   fChain->SetBranchAddress("w_Tau_ID_PNet", &w_Tau_ID_PNet, &b_w_Tau_ID_PNet);
   fChain->SetBranchAddress("w_Trigger", &w_Trigger, &b_w_Trigger);
   fChain->SetBranchAddress("w_SingleMuon_leg1", &w_SingleMuon_leg1, &b_w_SingleMuon_leg1);
   fChain->SetBranchAddress("w_CrossTrigger_Tau_leg2", &w_CrossTrigger_Tau_leg2, &b_w_CrossTrigger_Tau_leg2);
   fChain->SetBranchAddress("w_CrossTrigger_Muon_leg1", &w_CrossTrigger_Muon_leg1, &b_w_CrossTrigger_Muon_leg1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
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
