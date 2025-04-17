
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <iostream>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "/phenix/mpcex/lajoie/ISU/new_unfolding/RooUnfold/src/RooUnfoldResponse.h"
#include "TF1.h"
#include <TGraphAsymmErrors.h>
#include "TSystem.h"
#include "TVector3.h"

#define INCLUDE_CONSTITUENTS 
#define CONSTITUENT_MATCH_DR 0.025

#include <vector>
#include <algorithm>
#include <string.h>

//fastjet tools
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>

// LHAPDF
#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF; 

// ----------------------------------------------------------------------------
// EEC calculations
// ----------------------------------------------------------------------------
/*  Make sure the EEC path is pointed to wherever
 *  the repo is located. We'll need <utility> for
 *  std::pair.
 */

#include <utility>
#include "./include/PHEnergyCorrelator.h"

// ----------------------------------------------------------------------------

TRandom *myRand;

static const int maxRecoJets = 100;
static const int kMaxMatchedJets = 100;
static const int kMaxTruthJets = 200;

float vertex;
float centrality; 
int ip12_clock_cross; 

Int_t nRecoJets;
Int_t nMatchedJets; 
Int_t nTruthJets;

bool arm0FiredTrigger; 
bool arm1FiredTrigger; 

float r_pT[maxRecoJets];
float r_ml_pT[maxRecoJets];
float r_eta[maxRecoJets];
float r_phi[maxRecoJets];
float r_cf[maxRecoJets];
int   r_nc[maxRecoJets];
float r_zg[maxRecoJets]; 
float r_oang[maxRecoJets]; 

float tm_pT[kMaxMatchedJets];
float tm_zg[kMaxMatchedJets];
int   m_t_idx[kMaxMatchedJets]; 
int   m_r_idx[kMaxMatchedJets]; 
float dR[kMaxMatchedJets]; 

float t_pT[kMaxTruthJets];
float t_eta[kMaxTruthJets];
float t_phi[kMaxTruthJets];
float t_zg[kMaxTruthJets];
float t_cf[kMaxTruthJets];
float evt_Qsqr; 
float evt_x1;
float evt_x2;
float evt_s; 
float evt_t;
float evt_u; 
int process_id; 
int parton_id1; 
int parton_id2; 
int parton_id3; 
int parton_id4; 

float parton_id3_px; 
float parton_id3_py; 
float parton_id3_pz; 

float parton_id4_px; 
float parton_id4_py; 
float parton_id4_pz; 

float r_oang_noUE; 

int saved_embed = -1; 
TString saved_Suffix = "";

#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class std::vector<std::vector<float> >+;
#pragma link C++ class std::vector<std::vector<int> >+;
#endif

#ifdef INCLUDE_CONSTITUENTS

std::vector< std::vector<float> > *re_cs_charge; 
std::vector< std::vector<float> > *re_cs_z;
std::vector< std::vector<float> > *re_cs_jT; 
std::vector< std::vector<float> > *re_cs_dR; 
std::vector< std::vector<int> > *re_cs_UE; 
std::vector< std::vector<float> > *re_cs_eta; 
std::vector< std::vector<float> > *re_cs_phi; 

std::vector< std::vector<float> > *tr_cs_charge; 
std::vector< std::vector<float> > *tr_cs_z;
std::vector< std::vector<float> > *tr_cs_jT; 
std::vector< std::vector<float> > *tr_cs_dR; 
std::vector< std::vector<float> > *tr_cs_eta; 
std::vector< std::vector<float> > *tr_cs_phi; 

#endif

#define MAX_BINS 21
#define MAX_CKIN_RUNS 20

#include "bins.h"
#include "SimStats.h"

// DSSV14
extern "C" void dssvgupdate_(double *X, double *Q2, double *DUV, double *DDV, double *DUBAR, double *DDBAR, double *DSTR, double *DGLU);
extern "C" void dssvini_(); 

// MSTW08
//extern "C" double getonepdfnlo_(int *ih, double *x, double *q, int *f); 
//extern "C" void initalphas_(int *IORD, double *FR2, double *MUR, double *ASMUR, double *MC, double *MB, double *MT);
//extern "C" double alphas_(double *MUR); 

int gembed = 0; 

int GetBin(float pT, int NPTBINS, double *PTBINS){

  if( (pT <  PTBINS[0]) || (pT >=  PTBINS[NPTBINS]) ){
    //cout << "out of range pT =  " << pT << " " << PTBINS[0] << " " << PTBINS[NPTBINS] << endl;
    return -999;
  }

  int retVal = -1;
  for(int i = 0; i < NPTBINS; i++){
    if(pT >= PTBINS[i] && pT < PTBINS[i+1]){ 
      return i;
    }
  }//i

}//GetBin

void LinkMatchedVars(TTree *tReco, TTree *tMatch, TTree *tTrue){

  tReco->SetBranchAddress("nRecoJets", &nRecoJets);
  tMatch->SetBranchAddress("nMatchedJets", &nMatchedJets);
  tMatch->SetBranchAddress("vertex", &vertex);
  tMatch->SetBranchAddress("centrality", &centrality);
  tMatch->SetBranchAddress("ip12_clock_cross", &ip12_clock_cross);
  tMatch->SetBranchAddress("arm0FiredTrigger", &arm0FiredTrigger);
  tMatch->SetBranchAddress("arm1FiredTrigger", &arm1FiredTrigger);

  tMatch->SetBranchAddress("t_pT", tm_pT);
  tMatch->SetBranchAddress("t_z_g", tm_zg);
  tMatch->SetBranchAddress("t_idx", m_t_idx);
  tMatch->SetBranchAddress("r_idx", m_r_idx);
  tMatch->SetBranchAddress("dR", dR);

  tReco->SetBranchAddress("pT", r_pT);
  tReco->SetBranchAddress("ml_pT", r_ml_pT);
  tReco->SetBranchAddress("eta", r_eta);
  tReco->SetBranchAddress("phi", r_phi);
  tReco->SetBranchAddress("cf", r_cf);
  tReco->SetBranchAddress("nc", r_nc);
  tReco->SetBranchAddress("z_g", r_zg);
  tReco->SetBranchAddress("oang", r_oang);

  tTrue->SetBranchAddress("nTruthJets", &nTruthJets);
  tTrue->SetBranchAddress("pT", t_pT);
  tTrue->SetBranchAddress("eta", t_eta);
  tTrue->SetBranchAddress("phi", t_phi);
  tTrue->SetBranchAddress("evt_Qsqr", &evt_Qsqr);
  tTrue->SetBranchAddress("evt_x1", &evt_x1);
  tTrue->SetBranchAddress("evt_x2", &evt_x2);
  tTrue->SetBranchAddress("evt_s", &evt_s);
  tTrue->SetBranchAddress("evt_t", &evt_t);
  tTrue->SetBranchAddress("evt_u", &evt_u);
  tTrue->SetBranchAddress("process_id", &process_id);
  tTrue->SetBranchAddress("parton_id1", &parton_id1);
  tTrue->SetBranchAddress("parton_id2", &parton_id2);
  tTrue->SetBranchAddress("parton_id3", &parton_id3);
  tTrue->SetBranchAddress("parton_id4", &parton_id4);  
  tTrue->SetBranchAddress("parton_id3_px", &parton_id3_px);
  tTrue->SetBranchAddress("parton_id3_py", &parton_id3_py);
  tTrue->SetBranchAddress("parton_id3_pz", &parton_id3_pz);
  tTrue->SetBranchAddress("parton_id4_px", &parton_id4_px);
  tTrue->SetBranchAddress("parton_id4_py", &parton_id4_py);
  tTrue->SetBranchAddress("parton_id4_pz", &parton_id4_pz);
  tTrue->SetBranchAddress("z_g", t_zg);
  tTrue->SetBranchAddress("cf", t_cf);

#ifdef INCLUDE_CONSTITUENTS

  tReco->SetBranchAddress("cs_charge",&re_cs_charge); 
  tReco->SetBranchAddress("cs_z",&re_cs_z); 
  tReco->SetBranchAddress("cs_jT",&re_cs_jT); 
  tReco->SetBranchAddress("cs_dR",&re_cs_dR); 
  tReco->SetBranchAddress("cs_UE",&re_cs_UE); 
  tReco->SetBranchAddress("cs_eta",&re_cs_eta); 
  tReco->SetBranchAddress("cs_phi",&re_cs_phi); 

  tTrue->SetBranchAddress("cs_charge",&tr_cs_charge); 
  tTrue->SetBranchAddress("cs_z",&tr_cs_z); 
  tTrue->SetBranchAddress("cs_jT",&tr_cs_jT); 
  tTrue->SetBranchAddress("cs_dR",&tr_cs_dR); 
  tTrue->SetBranchAddress("cs_eta",&tr_cs_eta); 
  tTrue->SetBranchAddress("cs_phi",&tr_cs_phi); 

#endif

  return;

}//LinkMatchedVars

double JetCharge( int jetIdx ){

  // Calculate jet charge - arXiv:1209.2421v2

  const double kappa = 0.5; 
  double jc = 0; 

  for(unsigned int i=0; i<re_cs_z->at(jetIdx).size(); i++){
    if(re_cs_charge->at(jetIdx).at(i)==0.0) continue;   

    double ph = sqrt(pow(re_cs_z->at(jetIdx).at(i)*r_pT[jetIdx],2) + 
		     pow(re_cs_jT->at(jetIdx).at(i),2) ); 
    double htheta = 2.0*atan(exp(-re_cs_eta->at(jetIdx).at(i)));
    double cpt = ph*sin(htheta); 

    jc += (re_cs_charge->at(jetIdx).at(i))*pow(cpt, kappa); 

  }

  jc /= pow(r_pT[jetIdx],kappa); 

  return jc; 

}

double JetChargeTrue( int jetIdx ){

  // Calculate jet charge - arXiv:1209.2421v2

  const double kappa = 0.5; 
  double jc = 0; 

  for(unsigned int i=0; i<tr_cs_z->at(jetIdx).size(); i++){
    if(tr_cs_charge->at(jetIdx).at(i)==0.0) continue;   
    
    // 90% track eficiency
    if(myRand->Uniform()>0.9) continue;  

    double ph = sqrt(pow(tr_cs_z->at(jetIdx).at(i)*t_pT[jetIdx],2) + 
		     pow(tr_cs_jT->at(jetIdx).at(i),2) ); 
    
    // 300MeV low-momentum cutoff
    if(ph<0.300) continue; 

    double htheta = 2.0*atan(exp(-tr_cs_eta->at(jetIdx).at(i)));
    double cpt = ph*sin(htheta); 

    jc += (tr_cs_charge->at(jetIdx).at(i))*pow(cpt, kappa); 

  }

  jc /= pow(t_pT[jetIdx],kappa); 

  return jc; 

}

int getArm(float phi) {

  if (phi >= -TMath::PiOver2() && phi < TMath::PiOver2()) return 0; 
  if (phi >=  TMath::PiOver2() && phi < 3.0*TMath::PiOver2()) return 1; 
  return -1;

}

float dphiReduce(float dphi){
  if (dphi < -TMath::Pi()) 
    dphi += TMath::TwoPi();
  if (dphi >=  TMath::Pi()) 
    dphi -= TMath::TwoPi();
  return dphi; 
}//end phireduce()

int GetMaxPtIndex(float *pT, float *cf, float *phi, int nRecoJets, int runno, int NPTBINS, double *PTBINS, int AcceptFlag){

  float maxPt = 0.0;
  int indexMax = -1;
  for(int i = 0; i < nRecoJets; i++){

    if((runno==8) && (getArm(phi[i])==0)) continue; // East arm only
    //if((runno==8) && (cf[i]>=0.75)) continue;       // PPG184 cut
    if((runno==8) && ((cf[i]<0.3)||(cf[i]>=0.7))) continue;       // std cut

    //03/29/2022 - Run15 CfNC variation run with nC>=4, cf<0.65
    //restore the stricter cut
    if((runno==15) && (saved_embed==0) && (saved_Suffix=="_cfNC")) {
      if((r_cf[i]<0.3)||(r_cf[i]>=0.6)) continue;       // std variation cut
      if(r_nc[i]<5) continue;                           // std variation cut
    }

    // If AcceptFlag >=0 then use it to select the arm
    // That will be used for jet finding. 

    if(AcceptFlag>=0){
      if(getArm(phi[i])!=AcceptFlag) continue; 
    }

    if(pT[i] > maxPt){
      maxPt = pT[i];
      indexMax = i;
    }

  }

  return indexMax;

}//GetMaxPt


// Copy of soft drop algorithm from MakeJets
// Used here to calculate z_g without UE particles

#define ZG_CUT 0.1

Double_t SubstructureNoUE(int idx_max, float R){
   
  bool splitFound = false; 
  int sSize = 0;
  float pT_1 = -999.0; 
  float pT_2 = -999.0; 
  fastjet::PseudoJet BigSub;
  fastjet::PseudoJet SmallSub;

  // Cambridge-Aachen for jet substructure
  fastjet::JetDefinition *C_A = new fastjet::JetDefinition(fastjet::cambridge_algorithm, R, fastjet::E_scheme, fastjet::Best);
  
  // Copy constituents into FastJet structures 

  std::vector<fastjet::PseudoJet> SubCons; 

  for(unsigned int j=0; j<re_cs_z->at(idx_max).size(); j++){

    if(re_cs_UE->at(idx_max).at(j)!=0) continue; 

    float px =  re_cs_jT->at(idx_max).at(j)*cos(re_cs_phi->at(idx_max).at(j));  
    float py =  re_cs_jT->at(idx_max).at(j)*sin(re_cs_phi->at(idx_max).at(j));     
    float pz =  re_cs_z->at(idx_max).at(j)*r_pT[idx_max];

    float ptot = sqrt(px*px + py*py + pz*pz); 

    fastjet::PseudoJet pseudoFromConstit(px,py,pz,ptot); 
    SubCons.push_back(pseudoFromConstit); 
      
  }

  // z_g values:
  // -1 : uninitialized
  // -2 : singleton (one jet)
  // -3 : zero jets

  float zg_found = -1.0; 
  r_oang_noUE = -1.0; 

  do{
    fastjet::ClusterSequence cs(SubCons, *C_A);
    std::vector<fastjet::PseudoJet> subJets = cs.exclusive_jets_up_to(2);
    sSize = subJets.size();

    if(sSize==2){

      BigSub          = subJets[0];
      SmallSub        = subJets[1];
      pT_1            = BigSub.perp();
      pT_2            = SmallSub.perp();
      if(pT_1 < pT_2){
	fastjet::PseudoJet tTempJet = SmallSub;
	SmallSub      = BigSub;
	BigSub        = tTempJet;
	float tTempPt  = pT_2;
	pT_2          = pT_1;
	pT_1          = tTempPt;
      }

      float z_g             = pT_2 /(pT_1 + pT_2);
      if(z_g > ZG_CUT){
	zg_found = z_g;
	splitFound = true; 

	TVector3 bigVect(BigSub.px(),BigSub.py(),BigSub.pz()); 
	TVector3 smallVect(SmallSub.px(),SmallSub.py(),SmallSub.pz()); 

	r_oang_noUE = bigVect.Angle(smallVect);

      }
      else{
	SubCons = BigSub.constituents(); 
      }
	      
    }
    else{
      if(sSize==1)
        zg_found = -2;
      else
	zg_found = -3; 
    }

  }while((sSize==2)&&!splitFound); 

  delete C_A; 

  return zg_found;

}

double GetPolPartonicCrossSec(int procid, int id1, int id2, int id3, int id4, double alpha_s, double s, double t, double u)
{

  // partonic elments taken from
  // arXiv:hep-ph/9604220v1
  // R. Gastmanns and T. T. Wu, The Ubiquitous Photon: The Helicity Method for QED
  // and QCD (Clarnedon Press, Oxford, 1990).

  double parton_sigma = 0.0; 

  switch(procid)
  {
    
    // qq->qq
    case 11:
      if(id1==id2)
        parton_sigma = ((s*s - u*u)/(t*t)) + ((s*s - t*t)/(u*u)) - (2.0/3.0)*((s*s)/(t*u)); // OK
      else
        parton_sigma = ((s*s - u*u)/(t*t)); // OK
      break; 

    // qqbar->qqbar
    case 12:
      if((id1==-id2)&&(id3=-id4)&&((id1==id3)||(id1==id4)))
        parton_sigma = ((s*s - u*u)/(t*t)) - ((u*u + t*t)/(s*s)) + (2.0/3.0)*((u*u)/(s*t)); // OK
      else if((id1==-id2)&&(id3=-id4)&&((id1!=id3)&&(id1!=id4)))
        parton_sigma = - ((u*u + t*t)/(s*s)); // OK
      else if((id1!=-id2)&&(id3!=-id4)&&((id1==id3)||(id1==id4)))
	parton_sigma = ((s*s - u*u)/(t*t)); // OK
      break; 

    // qqbar->gg
    case 13:
      parton_sigma = (-8.0/3.0)*((t*t + u*u)/(t*u)) + 6.0*((t*t + u*u)/(s*s)); // OK
      break; 

    // qg->qg
    case 28:
      parton_sigma = (9.0/4.0)*(s*s - u*u)*( (1.0/(t*t)) - (4.0/(9.0*u*s)) ); // OK
      break; 

    // gg->qqbar
    case 53:
      parton_sigma = (9.0/8.0)*(t*t + u*u)*( (3.0/(4.0*s*s)) - (1.0/(3.0*u*t)) ); // OK 
      break; 

    // gg->gg
    case 68:
      parton_sigma = (9.0*9.0/16.0)*(-6.0 + (4.0*s*s/(u*t)) + ((u*s)/(t*t)) + ((2.0*u*t)/(s*s)) ); // OK 
      break; 

    default:

      break; 

  }

  // prefactor
  parton_sigma *= 4.0*TMath::Pi()*alpha_s*alpha_s/(9.0*s*s); 

  //cout << "GetPolPartonicCrossSec: parton_sigma = " << parton_sigma << endl; 

  return parton_sigma; 

}

double GetPartonicCrossSec(int procid, int id1, int id2, int id3, int id4, double alpha_s, double s, double t, double u)
{

  double parton_sigma = 0.0; 

  // R. Gastmanns and T. T. Wu, The Ubiquitous Photon: The Helicity Method for QED
  // and QCD (Clarnedon Press, Oxford, 1990).

  //cout << "process_id = " << procid << endl; 
  //cout << " partion id's : " << parton_id1 << " + " << parton_id2 << " -> " << parton_id3 << " + " << parton_id4 << endl;  

  switch(procid)
  {
    
    // qq->qq
    case 11:
      if(id1==id2)
        parton_sigma = ((s*s*s*s + t*t*t*t + u*u*u*u)/(t*t*u*u)) - ((8.0*s*s)/(3.0*t*u)); // OK (10.9)
      else
        parton_sigma = ((s*s + u*u)/(t*t)); // OK (10.7)
      break; 

    // qqbar->qqbar
    case 12:
      // q qbar -> q qbar
      if((id1==-id2)&&(id3=-id4)&&((id1==id3)||(id1==id4))) // OK (10.10)
        parton_sigma = (1.0/(t*t*s*s))*( s*s*s*s + t*t*t*t + u*u*u*u - ((8.0*s*t*u*u)/3.0) ); 
      // q qbar -> q' qbar'
      else if((id1==-id2)&&(id3=-id4)&&((id1!=id3)&&(id1!=id4))) // OK (10.8)
        parton_sigma = ((u*u + t*t)/(s*s));
      // q qbar' -> q qbar'
      else if((id1!=-id2)&&(id3!=-id4)&&((id1==id3)||(id1==id4))) // OK (10.7?)
	parton_sigma = ((s*s + u*u)/(t*t));
      break; 

    // qqbar->gg - OK (10.21)
    case 13:
      parton_sigma = (8.0/3.0)*(9.0/4.0)*(t*t + u*u)*(((4.0)/(9.0*t*u))-(1.0/(s*s))); 
      break; 

    // qg->qg - OK (10.14)
    case 28:
      parton_sigma = (9.0/4.0)*(s*s + u*u)*( (1.0/(t*t)) - (4.0/(9.0*s*u)) ); 
      break; 

    // gg->qqbar - OK (10.24)
    case 53:
      parton_sigma = (9.0/4.0)*(3.0/8.0)*(t*t + u*u)*( (4.0/(9.0*t*u)) - (1.0/(s*s)) ); 
      break; 

    // gg->gg - OK (10.25)
    case 68:
      parton_sigma = (9.0*9.0/(4.0*8.0))*(s*s*s*s + t*t*t*t + u*u*u*u)*(s*s + t*t + u*u)/(s*s*t*t*u*u); 
      break; 

    default:

      break; 

  }

  // prefactor
  parton_sigma *= 4.0*TMath::Pi()*alpha_s*alpha_s/(9.0*s*s); 

  //cout << "GetPartonicCrossSec: parton_sigma = " << parton_sigma << endl; 

  return parton_sigma; 

}

double GetPolEventWeight(int r_spinPat)
{

  double devt_x1 = evt_x1; 
  double devt_x2 = evt_x2; 
  double devt_Qsqr = evt_Qsqr; 
  double devt_Q = sqrt(evt_Qsqr); 

  static bool first = true; 
  static const LHAPDF::PDF* pdf = NULL; 

  // Initialize DSSV grids, MSTW2008 or LHAPDF
  if(first) {

    dssvini_();

    pdf = LHAPDF::mkPDF("NNPDF30_nlo_as_0118", 0);
    //pdf = LHAPDF::mkPDF("CT14MC1nlo", 0);
    //pdf = LHAPDF::mkPDF("MSTW2008nlo68cl", 0);

    // init MSTW08 to evaluate alpha_s
    // int IORD = 1; // NLO
    // double FR2 = 1.0; 
    // double MUR = 1.0; 
    // double ASMUR = 0.5;
    // double MC = 1.4; 
    // double MB = 4.75; 
    // double MT = 1.0e10; 
    // initalphas_(&IORD,&FR2,&MUR,&ASMUR,&MC,&MB,&MT);

    first = false; 
  }

  // DSSV14 PDF's

  double DUV1, DDV1, DUBAR1, DDBAR1, DSTR1, DGLU1; 
  dssvgupdate_(&devt_x1, &devt_Qsqr, &DUV1, &DDV1, &DUBAR1, &DDBAR1, &DSTR1, &DGLU1);

  double DUV2, DDV2, DUBAR2, DDBAR2, DSTR2, DGLU2; 
  dssvgupdate_(&devt_x2, &devt_Qsqr, &DUV2, &DDV2, &DUBAR2, &DDBAR2, &DSTR2, &DGLU2);

  // MSTW2008
  //double alpha_s = alphas_(&devt_Q); 
  
  // LHAPDF
  double alpha_s = pdf->alphasQ2(devt_Qsqr);  

  // need to look up partonic matrix elements
  double partonic = GetPartonicCrossSec(process_id,parton_id1,parton_id2,parton_id3,parton_id4,alpha_s,evt_s,evt_t,evt_u);
  double pol_partonic = GetPolPartonicCrossSec(process_id,parton_id1,parton_id2,parton_id3,parton_id4,alpha_s,evt_s,evt_t,evt_u);
	    
  double Weight = 0.0;
  double pdfW1 = 0.0; 
  double pdfW2 = 0.0; 
  int ih = 0; 

  double DV1 = 0.0;
  double UV1 = 0.0;
  double DBAR1 = 0.0;
  double UBAR1 = 0.0;
  double STR1 = 0.0;
  double GLU1 = 0.0;

  if(parton_id1==1){
    //int pid = 7; 
    //DV1 = getonepdfnlo_(&ih,&devt_x1,&devt_Q,&pid);
    // Valence is quarks minus antiquarks
    DV1 =  pdf->xfxQ2(parton_id1, devt_x1, devt_Qsqr) - 
           pdf->xfxQ2(-parton_id1, devt_x1, devt_Qsqr); 
    pdfW1 = (DDV1/DV1);
  }
  else if(parton_id1==2){
    //int pid = 8; 
    //UV1 = getonepdfnlo_(&ih,&devt_x1,&devt_Q,&pid);
    // Valence is quarks minus antiquarks
    UV1 =  pdf->xfxQ2(parton_id1, devt_x1, devt_Qsqr) - 
           pdf->xfxQ2(-parton_id1, devt_x1, devt_Qsqr); 
    pdfW1 = (DUV1/UV1);
  }
  else if(parton_id1==-1){
    //DBAR1 = getonepdfnlo_(&ih,&devt_x1,&devt_Q,&parton_id1);
    DBAR1 =  pdf->xfxQ2(parton_id1, devt_x1, devt_Qsqr); 
    pdfW1 = (DDBAR1/DBAR1);
  }
  else if(parton_id1==-2){
    //UBAR1 = getonepdfnlo_(&ih,&devt_x1,&devt_Q,&parton_id1);
    UBAR1 =  pdf->xfxQ2(parton_id1, devt_x1, devt_Qsqr); 
    pdfW1 = (DUBAR1/UBAR1);
  }
  else if(parton_id1==3){
    //STR1 = getonepdfnlo_(&ih,&devt_x1,&devt_Q,&parton_id1);
    STR1 =  pdf->xfxQ2(parton_id1, devt_x1, devt_Qsqr); 
    pdfW1 = (DSTR1/STR1);
  }
  else if(parton_id1==21){
    //int pid = 0; 
    //GLU1 = getonepdfnlo_(&ih,&devt_x1,&devt_Q,&pid);
    GLU1 =  pdf->xfxQ2(parton_id1, devt_x1, devt_Qsqr); 
    pdfW1 = (DGLU1/GLU1);
  }

  double DV2 = 0.0;
  double UV2 = 0.0;
  double DBAR2 = 0.0;
  double UBAR2 = 0.0;
  double STR2 = 0.0;
  double GLU2 = 0.0;

  if(parton_id2==1){
    //int pid = 7; 
    //DV2 = getonepdfnlo_(&ih,&devt_x2,&devt_Q,&pid);
    // Valence is quarks minus antiquarks
    DV2 =  pdf->xfxQ2(parton_id2, devt_x2, devt_Qsqr) - 
           pdf->xfxQ2(-parton_id2, devt_x2, devt_Qsqr); 
    pdfW2 = (DDV2/DV2);
  }
  else if(parton_id2==2){
    //int pid = 8; 
    //UV2 = getonepdfnlo_(&ih,&devt_x2,&devt_Q,&pid); 
    // Valence is quarks minus antiquarks
    UV2 =  pdf->xfxQ2(parton_id2, devt_x2, devt_Qsqr) - 
           pdf->xfxQ2(-parton_id2, devt_x2, devt_Qsqr);
    pdfW2 = (DUV2/UV2);
  }
  else if(parton_id2==-1){
    //DBAR2 = getonepdfnlo_(&ih,&devt_x2,&devt_Q,&parton_id2);
    DBAR2 =  pdf->xfxQ2(parton_id2, devt_x2, devt_Qsqr); 
    pdfW2 = (DDBAR2/DBAR2);
  }
  else if(parton_id2==-2){
    //UBAR2 = getonepdfnlo_(&ih,&devt_x2,&devt_Q,&parton_id2);
    UBAR2 =  pdf->xfxQ2(parton_id2, devt_x2, devt_Qsqr); 
    pdfW2 = (DUBAR2/UBAR2); 
  }
  else if(parton_id2==3){
    //STR2 = getonepdfnlo_(&ih,&devt_x2,&devt_Q,&parton_id2);
    STR2 =  pdf->xfxQ2(parton_id2, devt_x2, devt_Qsqr); 
    pdfW2 = (DSTR2/STR2);
  }
  else if(parton_id2==21){
    //int pid = 0; 
    //GLU2 = getonepdfnlo_(&ih,&devt_x2,&devt_Q,&pid);
    GLU2 =  pdf->xfxQ2(parton_id2, devt_x2, devt_Qsqr); 
    pdfW2 = (DGLU2/GLU2);
  }

  // Positivity checks

  if(fabs(pdfW1)>1.0){
    cout << " WARNING: |pdfW1|>1, pdfW1 = " << pdfW1 << " process = " << process_id 
    	 << " id1 = " << parton_id1 << " id2 = " << parton_id2 
    	 << " id3 = " << parton_id3 << " id4 = " << parton_id4 
    	 << " x1 = " << devt_x1 << " Q^2 = " << devt_Qsqr << endl; 
    // if(parton_id1==-1)
    //   cout << " DDBAR1 = " << DDBAR1 << " DBAR1 = " << DBAR1 << endl; 
    // if(parton_id1==1)
    //   cout << " DDV1 = " << DDV1 << " DV1 = " << DV1 << endl; 
    cout << endl; 
    pdfW1 = copysign(1.0,pdfW1);  
  }

  if(fabs(pdfW2)>1.0){
    cout << " WARNING: |pdfW2|>1, pdfW2 = " << pdfW2 << " process = " << process_id 
    	 << " id1 = " << parton_id1 << " id2 = " << parton_id2 
    	 << " id3 = " << parton_id3 << " id4 = " << parton_id4 
    	 << " x2 = " << devt_x2 << " Q^2 = " << devt_Qsqr << endl; 
    // if(parton_id2==-1)
    //   cout << " DDBAR2 = " << DDBAR2 << " DBAR2 = " << DBAR2 << endl; 
    // if(parton_id2==1)
    //   cout << " DDV2 = " << DDV2 << " DV2 = " << DV2 << endl; 
    cout << endl; 
    pdfW2 = copysign(1.0,pdfW2);  
  }

  if(partonic>0.0){

    // Watch for processes that have a partonic asymmetry of -1 always. 
    if((fabs(pol_partonic/partonic)>1.0) && (process_id!=53) && (process_id!=13) &&
       !((process_id==12) && (parton_id1==-parton_id2)&&(parton_id3=-parton_id4)&&((parton_id1!=parton_id3)&&(parton_id1!=parton_id4)) ) ){
      cout << " WARNING: |pol_partonic/partonic|>1, = " << (pol_partonic/partonic) << " process = " << process_id 
	   << " id1 = " << parton_id1 << " id2 = " << parton_id2 
	   << " id3 = " << parton_id3 << " id4 = " << parton_id4 
	   << " x2 = " << devt_x2 << " Q^2 = " << devt_Qsqr << endl; 
      Weight = pdfW1*pdfW2; 
    }
    else
      Weight = pdfW1*pdfW2*(pol_partonic/partonic); 

  }
  else{
    cout << " WARNING: partonic = " << partonic << " process = " << process_id 
	 << " id1 = " << parton_id1 << " id2 = " << parton_id2 
	 << " id3 = " << parton_id3 << " id4 = " << parton_id4 << endl << endl; 
    Weight = 0.0; 
  }

  return Weight; 

}

int ConnectJetToParton(int iJet){

  TVector3 p3(parton_id3_px,parton_id3_py,parton_id3_pz); 
  TVector3 p4(parton_id4_px,parton_id4_py,parton_id4_pz); 

  double dist3 = sqrt(pow(t_eta[iJet] - p3.Eta(),2) + pow(dphiReduce(p3.Phi() - t_phi[iJet]),2)); 
  double dist4 = sqrt(pow(t_eta[iJet] - p4.Eta(),2) + pow(dphiReduce(p4.Phi() - t_phi[iJet]),2)); 

  if((dist3<dist4)&&(dist3<0.3)) 
    return parton_id3; 
  else if((dist4<dist3)&&(dist4<0.3)) 
    return parton_id4; 
  else
    return -9999; 
}

void jetAnaSim(int runno=12, float R = 0.3, int embed = 0, float centLow = 0.0, float centHigh = 20.0, 
	       double *R_even_out = NULL, double *R_odd_out = NULL, double *R_out = NULL, const int weightALL = 0, const int weightNLO = 1,
	       const int useML = 0, int AcceptFlag = -1, int halfstats = 0, int p_or_h_flag = 1, std::string inSuffix=""){

  // --------------------------------------------------------------------------
  // EEC calculations
  // --------------------------------------------------------------------------

// define flags to turn off certain calculations
#define doTrueEEC 1
#define doTrueEECChargedOnly 0
#define doRecoEEC 0
#define doRecoEECChargedOnly 0

// define flags to turn on/off certain binnings
#define doJetCFBins 0
#define doJetChargeBins 0

  // pt jet bins
  std::vector< std::pair<float, float> > ptJetBins;
  ptJetBins.push_back( std::make_pair(5., 10.) );
  ptJetBins.push_back( std::make_pair(10., 15.) );
  ptJetBins.push_back( std::make_pair(15., 20.) );

  // cf jet bins
  std::vector< std::pair<float, float> > cfJetBins;
  cfJetBins.push_back( std::make_pair(0., 1.0) );

  // jet charge bins
  std::vector< std::pair<float, float> > chJetBins;
  chJetBins.push_back( std::make_pair(-100., 0.0) );
  chJetBins.push_back( std::make_pair(0.0, 100.) );

  // now declare calculators
  //   - 1st argument: quantity used for weights (pt, et, or e)
  //   - 2nd argument: power to raise weights to
  PHEC::Calculator trueEEC( PHEC::Type::Pt, 1.0 );
  PHEC::Calculator recoEEC( PHEC::Type::Pt, 1.0 );

  // set histogram tags
  trueEEC.SetHistTag( "TrueJet" );
  recoEEC.SetHistTag( "RecoJet" );

  // set pt, cf, charge jet bins
  trueEEC.SetPtJetBins( ptJetBins );
  recoEEC.SetPtJetBins( ptJetBins );
  if (doJetCFBins) {
    trueEEC.SetCFJetBins( cfJetBins );
    recoEEC.SetCFJetBins( cfJetBins );
  }
  if (doJetChargeBins) {
    trueEEC.SetChargeBins( chJetBins );
    recoEEC.SetChargeBins( chJetBins );
  }

  // turn on spin sorting
  trueEEC.SetDoSpinBins( true );
  recoEEC.SetDoSpinBins( true );

  // run initialization routine to generate 
  // desired histograms
  //   - 1st argument: turn on/off 2-point histograms 
  //   - 2nd argument: turn on/off 3-point histograms (TODO)
  //   - 3rd argument: turn on/off lambda EEC histograms (TODO)
  trueEEC.Init(true, false, false);
  recoEEC.Init(true, false, false);

  // --------------------------------------------------------------------------

  //TString PPAdd = "_Pythia8"; 
  TString PPAdd = ""; 
  TString HIAdd = "EmbedPythia"; 

  saved_embed = embed; 
  saved_Suffix = inSuffix; 

  bool Run12EmbedJewel = false; 
  if((runno==12)&&(embed==2)) Run12EmbedJewel = true; // 1 = Pythia, 2 = Jewel

  if(Run12EmbedJewel){
    cout << "Run12EmbedJewel = TRUE" << endl; 
  }

  if( (runno!=12) && (runno!=13) && (runno!=8) && (runno!=15) ){
    cout << "ERROR: BAD runno = " << runno << endl;
    return; 
  }
  
  gembed = embed; 

  int cmsE = 510; 
  if(runno==13){
    cmsE = 510; 
  }
  else if((runno==12)||(runno==8)||(runno==15)){
    cmsE = 200; 
  }
 
  if(embed && useML){
    cout << "ERROR: Machine Learning not currently supported for embedded" << endl;  
    return; 
  }

  if(useML && !((runno==12)||(runno==13))){
    cout << "ERROR: Machine Learning not currently supported for run = " << runno << endl;  
    return; 
  }

  TString HIString = ""; 
  if(embed){
    if(runno==8) HIString = "dAu"; 
    if(runno==12) HIString = "CuAu"; 
    if(runno==15) HIString = "pAu"; 
  }

  TString Suffix = inSuffix; 

  TString RString = "";
  if(fabs(R-0.35)<0.01)
    RString = "R035"; 
  else if(fabs(R-0.3)<0.01)
    RString = "R03"; 
  else if (fabs(R-0.2)<0.01)
    RString = "R02";
  else{
    cout << "ERROR: BAD R = " << R << endl;
    return; 
  }

  unsigned int NPTBINS = 0;
  double PTBINS[MAX_BINS];
  unsigned int NPTBINS_TRUE = 0;
  double PTBINS_TRUE[MAX_BINS];

  unsigned int nCkinRuns = 0; 
  float vtxCut = 0.0;  
  int ckins[MAX_CKIN_RUNS] = {0.0}; 
  float crossSec[MAX_CKIN_RUNS] = {0.0}; 
  int nFiles[MAX_CKIN_RUNS] = {0.0};             
  unsigned int evPerFile[MAX_CKIN_RUNS] = {0.0}; 
  float tpT_cutoff[MAX_CKIN_RUNS] = {0.0};   
  float rpT_cutoff[MAX_CKIN_RUNS] = {0.0}; 

  if(runno==13){

    nCkinRuns = nRuns_r13; 
    for(unsigned int i=0; i<nCkinRuns; i++){
      ckins[i] = ckins_r13[i]; 
      crossSec[i] = crossSec_r13[i]; 
      nFiles[i] = nFiles_r13[i]; 
      evPerFile[i] = evPerFile_r13[i];

      if(p_or_h_flag==1)
        tpT_cutoff[i] = tpT_cutoff_r13[i];
      else
        tpT_cutoff[i] = tpT_p_cutoff_r13[i];

      rpT_cutoff[i] = rpT_cutoff_r13[i]; 
    }
    
    NPTBINS = NPTBINS_RECO_13; 
    for(unsigned int i = 0; i<NPTBINS_RECO_13+1; i++){
      PTBINS[i] = PTBINS_RECO_13[i];  
    }
    NPTBINS_TRUE = NPTBINS_TRUE_13; 
    for(unsigned int i = 0; i<NPTBINS_TRUE_13+1; i++){
      PTBINS_TRUE[i] = PTBINS_TRUE_13[i];  
    }
 
    vtxCut = vtxCut_r13; 

  }
  else if(runno==12){
     
    if(embed==0){
      nCkinRuns = nRuns_r12; 
      for(unsigned int i=0; i<nCkinRuns; i++){
	ckins[i] = ckins_r12[i]; 
	crossSec[i] = crossSec_r12[i]; 
	nFiles[i] = nFiles_r12[i]; 
	evPerFile[i] = evPerFile_r12[i]; 

	if(p_or_h_flag==1)
	  tpT_cutoff[i] = tpT_cutoff_r12[i];
	else
	  tpT_cutoff[i] = tpT_p_cutoff_r12[i];

        rpT_cutoff[i] = rpT_cutoff_r12[i]; 
      }
    
      NPTBINS = NPTBINS_RECO_12; 
      for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	PTBINS[i] = PTBINS_RECO_12[i];  
      }
      NPTBINS_TRUE = NPTBINS_TRUE_12; 
      for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
      }

      vtxCut = vtxCut_r12;
    }
    else{

      if(Run12EmbedJewel){

	nCkinRuns = nRuns_r12ej; 
	for(unsigned int i=0; i<nCkinRuns; i++){
	  ckins[i] = ckins_r12ej[i]; 
	  crossSec[i] = crossSec_r12ej[i]; 
	  nFiles[i] = nFiles_r12ej[i]; 
	  evPerFile[i] = evPerFile_r12ej[i]; 
          tpT_cutoff[i] = tpT_cutoff_r12ej[i]; 
          rpT_cutoff[i] = rpT_cutoff_r12ej[i]; 
	}
    
	NPTBINS = NPTBINS_RECO_12; 
	for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	  PTBINS[i] = PTBINS_RECO_12[i];  
	}
	NPTBINS_TRUE = NPTBINS_TRUE_12; 
	for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	  PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
	}

	vtxCut = vtxCut_r12ej;

      }
      else{

	nCkinRuns = nRuns_r12e; 
	for(unsigned int i=0; i<nCkinRuns; i++){
	  ckins[i] = ckins_r12e[i]; 
	  crossSec[i] = crossSec_r12e[i]; 
	  nFiles[i] = nFiles_r12e[i]; 
	  evPerFile[i] = evPerFile_r12e[i];
          tpT_cutoff[i] = tpT_cutoff_r12e[i]; 
          rpT_cutoff[i] = rpT_cutoff_r12e[i]; 
	}
    
	NPTBINS = NPTBINS_RECO_12; 
	for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	  PTBINS[i] = PTBINS_RECO_12[i];  
	}
	NPTBINS_TRUE = NPTBINS_TRUE_12; 
	for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	  PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
	}

	vtxCut = vtxCut_r12e;

      }


    }

  }
  else if(runno==8){

    if(embed==0){
      nCkinRuns = nRuns_r8; 
      for(unsigned int i=0; i<nCkinRuns; i++){
	ckins[i] = ckins_r8[i]; 
	crossSec[i] = crossSec_r8[i]; 
	nFiles[i] = nFiles_r8[i]; 
	evPerFile[i] = evPerFile_r8[i]; 

	if(p_or_h_flag==1)
          tpT_cutoff[i] = tpT_cutoff_r8[i];
	else
          tpT_cutoff[i] = tpT_p_cutoff_r8[i];

        rpT_cutoff[i] = rpT_cutoff_r8[i]; 
      }

      NPTBINS = NPTBINS_RECO_12; 
      for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	PTBINS[i] = PTBINS_RECO_12[i];  
      }
      NPTBINS_TRUE = NPTBINS_TRUE_12; 
      for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
      }

      vtxCut = vtxCut_r8; 
    }
    else{
      nCkinRuns = nRuns_r8e; 
      for(unsigned int i=0; i<nCkinRuns; i++){
	ckins[i] = ckins_r8e[i]; 
	crossSec[i] = crossSec_r8e[i]; 
	nFiles[i] = nFiles_r8e[i]; 
	evPerFile[i] = evPerFile_r8e[i];
        tpT_cutoff[i] = tpT_cutoff_r8e[i];  
        rpT_cutoff[i] = rpT_cutoff_r8e[i];  
      }

      NPTBINS = NPTBINS_RECO_12; 
      for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	PTBINS[i] = PTBINS_RECO_12[i];  
      }
      NPTBINS_TRUE = NPTBINS_TRUE_12; 
      for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
      }

      vtxCut = vtxCut_r8e; 
    }
 

  }
  else if(runno==15){

    if(embed==0){
      nCkinRuns = nRuns_r15; 
      for(unsigned int i=0; i<nCkinRuns; i++){
	ckins[i] = ckins_r15[i]; 

	if(PPAdd=="Pythia8")
	  crossSec[i] = crossSecP8_r15[i]; 
	else
	  crossSec[i] = crossSec_r15[i]; 

	nFiles[i] = nFiles_r15[i]; 
	evPerFile[i] = evPerFile_r15[i]; 

	if(p_or_h_flag==1)
	  tpT_cutoff[i] = tpT_cutoff_r15[i];
	else
	  tpT_cutoff[i] = tpT_p_cutoff_r15[i];

        rpT_cutoff[i] = rpT_cutoff_r15[i]; 
      }

      NPTBINS = NPTBINS_RECO_12; 
      for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	PTBINS[i] = PTBINS_RECO_12[i];  
      }
      NPTBINS_TRUE = NPTBINS_TRUE_12; 
      for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
      }

      vtxCut = vtxCut_r15; 
    }
    else{
      nCkinRuns = nRuns_r15e; 
      for(unsigned int i=0; i<nCkinRuns; i++){
	ckins[i] = ckins_r15e[i]; 

	if(HIAdd=="EmbedHerwig"){	  
	  crossSec[i] = crossSecH_r15e[i];
	  nFiles[i] = nFiles_r15e[i]; 
	}
	else if(HIAdd=="EmbedPythia8"){	  
	  crossSec[i] = crossSecP8_r15e[i];
	  nFiles[i] = nFiles_r15e[i]; 
	}
	else{
	  crossSec[i] = crossSec_r15e[i]; 
	  nFiles[i] = nFiles_r15e[i]; 
	}

	evPerFile[i] = evPerFile_r15e[i];
        tpT_cutoff[i] = tpT_cutoff_r15e[i];  
        rpT_cutoff[i] = rpT_cutoff_r15e[i];  

      }

      NPTBINS = NPTBINS_RECO_12; 
      for(unsigned int i = 0; i<NPTBINS_RECO_12+1; i++){
	PTBINS[i] = PTBINS_RECO_12[i];  
      }
      NPTBINS_TRUE = NPTBINS_TRUE_12; 
      for(unsigned int i = 0; i<NPTBINS_TRUE_12+1; i++){
	PTBINS_TRUE[i] = PTBINS_TRUE_12[i];  
      }

      vtxCut = vtxCut_r15e; 
    }


  }

  const int sNPTBINS = NPTBINS;
  const unsigned int snCkinRuns = nCkinRuns;  

  // Vertex weighting file
  TFile *fVtx = new TFile(Form("vertexRatio_Run%i.root",runno), "read");
  TH1F *hVtxRatio = (TH1F*)fVtx->Get("hVertexRatio");
  if(!hVtxRatio){
    cout << " ERROR - missing vertex weights!" << endl; 
    return;
  }
  // no weighting for embedded events (use data vertex)
  if(embed>0){
    hVtxRatio->Reset(); 
    for(unsigned int i=0; i<(unsigned)hVtxRatio->GetNbinsX(); i++){
      hVtxRatio->SetBinContent(i+1,1.0); 
    }
  }

  // Output file

  TString MLString = ""; 
  if(useML) MLString = "_ML"; 

  // Update the output name based on arm selection
  TString oSuffix = Suffix; 
  if(AcceptFlag==0) oSuffix = "_Arm0" + Suffix; 
  if(AcceptFlag==1) oSuffix = "_Arm1" + Suffix; 
  
  if(doRecoEEC) oSuffix += "_RecoEEC"; 
  if(doTrueEEC) oSuffix += "_TrueEEC"; 

  TString HalfString = ""; 
  if(halfstats==1) HalfString = "_half"; 

  TFile *fOut = NULL; 
  if(embed==0){
    if(p_or_h_flag==1)
      fOut = new TFile(Form("jetAnaSimRun%i%s_pp_R%4.2f%s%s.root", runno, HalfString.Data(), R, MLString.Data(),oSuffix.Data()), "Recreate");
    else
      fOut = new TFile(Form("jetAnaSimRun%i%s_pp_R%4.2f%s%s_pjets.root", runno, HalfString.Data(), R, MLString.Data(),oSuffix.Data()), "Recreate");
  }
  else{
    if(p_or_h_flag==1)
      fOut = new TFile(Form("jetAnaSimRun%i%s_%s_R%4.2f_%i_%i%s%s.root", runno, HalfString.Data(), HIString.Data(), R, 
			  (int)centLow, (int)centHigh, MLString.Data(),oSuffix.Data()), "Recreate");
    else
      fOut = new TFile(Form("jetAnaSimRun%i%s_%s_R%4.2f_%i_%i%s%s_pjets.root", runno, HalfString.Data(), HIString.Data(), R, 
			  (int)centLow, (int)centHigh, MLString.Data(),oSuffix.Data()), "Recreate");
  }

  // Histograms

  TH1D *hCSWeight = new TH1D("hCSWeight", "cross section weight",nCkinRuns,0.0,(float)nCkinRuns);
  TH1D *hSimLumi = new TH1D("hSimLumi", "effective simulated Lumi (mb-1)",nCkinRuns,0.0,(float)nCkinRuns);

  TH1F *hReco = new TH1F("hReco", "Highest pT Reco Jets", NPTBINS, PTBINS);
  hReco->Sumw2();
   
  TH1F *hReco_ERTfired = new TH1F("hReco_ERTfired", "Highest pT Reco Jets", NPTBINS, PTBINS);
  hReco_ERTfired->Sumw2();

  TH2D *hRecoJetCharge_ERTfired = new TH2D("hRecoJetCharge_ERTfired", "Jet Charge vs. Reco Pt", NPTBINS, PTBINS, NJCBINS, JCBINS );
  hRecoJetCharge_ERTfired->Sumw2();

  TH1F *hTrue_Jets = new TH1F("hTrue_Jets", "Truth Jets", NPTBINS, PTBINS);
  hTrue_Jets->Sumw2();

  TH1F *hTrue_ERTfired = new TH1F("hTrue_ERTfired", "Truth Jets", NPTBINS, PTBINS);
  hTrue_ERTfired->Sumw2();

  TH1F *hTrueSamePat = new TH1F("hTrueSamePat", "Same spin: ++/--", NPTBINS, PTBINS);
  hTrueSamePat->Sumw2();
  TH1F *hTrueOppPat = new TH1F("hTrueOppPat", "Opp spin: +-/-+", NPTBINS, PTBINS);
  hTrueOppPat->Sumw2();

  TH1F *hReco_bin[sNPTBINS];
  TH1F *hReco_bin_ERTfired[sNPTBINS];
  for(unsigned int i = 0; i < (unsigned int) sNPTBINS; i++){
    hReco_bin[i] = new TH1F(Form("hRecoBin_%i", i), Form("Bin %i", i+1), 100, PTBINS[i], PTBINS[i+1]);
    hReco_bin_ERTfired[i] = new TH1F(Form("hRecoBinERTfired_%i", i), Form("Bin %i", i+1), 100, PTBINS[i], PTBINS[i+1]);
  }//i

  TH1D *hVertex = new TH1D("hVertex", "Vertex Distribution; Vertex (cm)", 30, -30.0, 30.0);
  TH1D *hVertexCut = new TH1D("hVertexCut", "Vertex Distribution; Vertex (cm)", 30, -30.0, 30.0);
  TH1D *hVertexJets = new TH1D("hVertexJets", "Vertex Distribution; Vertex (cm)", 30, -30.0, 30.0);

  TH2D *hZg = new TH2D("hZg", "z_{g} vs. reco p_{T}", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZg->Sumw2(); 
  TH2D *hOang = new TH2D("hOang", "oang vs. reco p_{T}", NPTBINS, PTBINS, NOANGBINS, OANGBINS);
  hOang->Sumw2(); 

  TH2D *hZgNoUE = new TH2D("hZgNoUE", "z_{g} vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZgNoUE->Sumw2(); 
  TH2D *hZgNoUETrig = new TH2D("hZgNoUETrig", "z_{g} vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NZGBINS, ZGBINS);
  hZgNoUETrig->Sumw2(); 

  TH2D *hOangNoUE = new TH2D("hOangNoUE", "oang vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NOANGBINS, OANGBINS);
  hOangNoUE->Sumw2(); 
  TH2D *hOangNoUETrig = new TH2D("hOangNoUETrig", "oang vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NOANGBINS, OANGBINS);
  hOangNoUETrig->Sumw2(); 

  TH2D *hFFZ = new TH2D("hFFZ", "z vs. reco p_{T}", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZ->Sumw2(); 

  TH2D *hFFZNoUE = new TH2D("hFFZNoUE", "z vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZNoUE->Sumw2(); 
  TH2D *hFFZNoUETrig = new TH2D("hFFZNoUETrig", "z vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NFFZBINS, FFZBINS);
  hFFZNoUETrig->Sumw2(); 

  TH2D *hFFXI = new TH2D("hFFXI", "z vs. reco p_{T}", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXI->Sumw2(); 

  TH2D *hFFXINoUE = new TH2D("hFFXINoUE", "z vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXINoUE->Sumw2(); 
  TH2D *hFFXINoUETrig = new TH2D("hFFXINoUETrig", "z vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NFFXIBINS, FFXIBINS);
  hFFXINoUETrig->Sumw2(); 

  TH2D *hFFJT = new TH2D("hFFJT", "jT/pT vs. reco p_{T}", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJT->Sumw2(); 

  TH2D *hFFJTNoUE = new TH2D("hFFJTNoUE", "jT/pT vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJTNoUE->Sumw2(); 
  TH2D *hFFJTNoUETrig = new TH2D("hFFJTNoUETrig", "jT/pT vs. reco p_{T} (no UE)", NPTBINS, PTBINS, NFFJTBINS, FFJTBINS);
  hFFJTNoUETrig->Sumw2(); 

  TH2D *hFFdR = new TH2D("hFFdR", "dR vs. reco p_{T}", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdR->Sumw2(); 

  TH2D *hFFdRNoUE = new TH2D("hFFdRNoUE", "dR vs. reco p_{T} (noUE)", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdRNoUE->Sumw2(); 
  TH2D *hFFdRNoUETrig = new TH2D("hFFdRNoUETrig", "dR vs. reco p_{T} (noUE)", NPTBINS, PTBINS, NFFDRBINS, FFDRBINS);
  hFFdRNoUETrig->Sumw2();   

  // jet yield in run number vs. IP12 crossing

  TH2D *JetYield = new TH2D("JetYield","yield in runnumber vs. ip12 crossing",121,-0.5,120.5,2208,387538.5,396367.5); 
  TH2D *JetYieldNar = new TH2D("JetYieldNar","yield in runnumber vs. ip12 crossing",121,-0.5,120.5,20,391859.5,391879.5); 
  
  // charged fraction, etc.

  int ncfbins = 100; 
  double cfbins[101]; 
  for(unsigned int i=0; i<101; i++) cfbins[i] = (1.0/100.0)*i;

  TH2D *hCF = new TH2D("hCF", "Charged Fraction",NPTBINS, PTBINS, ncfbins, cfbins);
  hCF->Sumw2();

  int ncbins = 21; 
  double cbins[22]; 
  for(unsigned int i=0; i<22; i++) cbins[i] = -0.5 + 1.0*i;

  TH2D *hNC = new TH2D("hNC", "Jet NC",NPTBINS,PTBINS,ncbins,cbins);
  hNC->Sumw2();

  TH1F *hPhi = new TH1F("hPhi", "Jet Phi",100,-TMath::PiOver2(),3.0*TMath::PiOver2());
  hPhi->Sumw2();

  // Test histogram

  TH1D *hEmbed = new TH1D("hEmbed","",2,-0.5,1.5); 

  // Random number generator (for spin pattern)
  myRand = new TRandom3(); 
  myRand->SetSeed(0); // new random requence

  // check for consistency/relative luminosity
  double nEvenOdd[2] = {0.0,0.0}; 
  double nSpinPat[2][7] = {{0.0,0.0,0.0,0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0,0.0,0.0,0.0}}; 

  // Weight function to reweight PYTHIA to partonic NLO

  double par_510[10]  = {1.02908,-0.0227364,0.00284752,-7.91447e-5,1.00291e-6,-6.02391e-09,9.36645e-12,7.75787e-14,-4.11498e-16,6.03721e-19}; 
  double par_200[10]  = {-0.694587,0.188099,-0.00729453,0.000175712,-2.3791e-06,1.31339e-08,0.0,0.0,0.0,0.0}; 
  double noWeight[10] = {1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
  TF1 *NLOWeight = new TF1("NLOWeight","pol9",5.0,255.0); 
  for(unsigned int i=0; i<10; i++){
    if((runno==13) && (weightNLO==1) && (embed==0) && (p_or_h_flag==0) ) 
      NLOWeight->SetParameter(i,par_510[i]);
    else if((runno==12) && (weightNLO==1) && (embed==0) && (p_or_h_flag==0))
      NLOWeight->SetParameter(i,par_200[i]);
    else
      NLOWeight->SetParameter(i,noWeight[i]);
  }

  // Add in the ability to combine centrality ranges to create a minbias unfolding
  // for Run-8/12

  int nMBCombine = 1; 
  float centlow8[4] = {0.0,20.0,40.0,60.0}; 
  float centhigh8[4] = {20.0,40.0,60.0,88.0}; 
  float centlow12[4] = {0.0,20.0,40.0,60.0}; 
  float centhigh12[4] = {20.0,40.0,60.0,93.0}; 
  float centlow15[4] = {0.0,20.0,40.0,60.0}; 
  float centhigh15[4] = {20.0,40.0,60.0,84.0}; 
  float centlow[4] = {0.0,0.0,0.0,0.0}; 
  float centhigh[4] = {0.0,0.0,0.0,0.0}; 
  float centWeightFactor8[4] = {9.66417e+07, 9.66911e+07, 9.61213e+07, 9.66911e+07}; 
  float centWeightFactor12[4] = {9.66417e+07, 9.66911e+07, 9.61213e+07, 9.66911e+07}; // TODO: update for Run-12
  float centWeightFactor15[4] = {2.45504e+08, 2.45494e+08,  2.45504e+08, 2.45504e+08}; 
  float centWeightFactor[4] = {0.0,0.0,0.0,0.0}; 

  // Default configuration - just use what you are given 

  centlow[0] = centLow; 
  centhigh[0] = centHigh; 
  centWeightFactor[0] = 1.0; 
  
  // Test Run-8/12 for MB combine values

  if((runno==8) && (centLow==0.0) && (centHigh==88.0) && (embed>0) ){
    nMBCombine = 4; 
    for(int i=0; i<nMBCombine; i++){
      centlow[i] = centlow8[i]; 
      centhigh[i] = centhigh8[i]; 
      centWeightFactor[i] = centWeightFactor8[i]; 
    }
  }
  else if( (runno==12) && (centLow==0.0) && (centHigh==93.0) && (embed>0) ){
    nMBCombine = 4; 
    for(int i=0; i<nMBCombine; i++){
      centlow[i] = centlow12[i]; 
      centhigh[i] = centhigh12[i]; 
      centWeightFactor[i] = centWeightFactor12[i]; 
    }
  }
  // else if( (runno==15) && (centLow==0.0) && (centHigh==84.0) ){
  //   nMBCombine = 4; 
  //   for(int i=0; i<nMBCombine; i++){
  //     centlow[i] = centlow15[i]; 
  //     centhigh[i] = centhigh15[i]; 
  //     centWeightFactor[i] = centWeightFactor15[i]; 
  //   }
  // }

  // Loop over each CKIN, create a TChain and process events

  double totalCs = 0.0; 
  double csWeightedEvents = 0.0; 
  double fullyWeightedEvents = 0.0; 
  double fullyWeightedCutEvents = 0.0; 

  for(int ic = 0 ; ic<nMBCombine; ic++){

    // Everything is weighted relative to one CKIN range
  
    TString p_or_h = "_hadrons"; 
    if(p_or_h_flag==0) p_or_h = "_partons";

    TString fileLoc =  Form("/phenix/crs/phnxreco/lajoie/SimOut/Run%i_PP%s/%i", runno, PPAdd.Data(), cmsE);
    //if(runno==15) fileLoc =  Form("/phenix/crs/phnxreco/jrunchey/ISUSim/FullSim/Run%i_PP/%i", runno, cmsE);

    if(embed) { 
      if(runno==8) fileLoc =  Form("/phenix/crs/phnxreco/lajoie/SimOut/Run%i_DAU%s_%i_%i/%i", runno, HIAdd.Data(), (int)centlow[ic], (int)centhigh[ic], cmsE);
      if(runno==15) fileLoc =  Form("/phenix/crs/phnxreco/lajoie/SimOut/Run%i_PAU%s_%i_%i/%i", runno, HIAdd.Data(), (int)centlow[ic], (int)centhigh[ic], cmsE);
      if((runno==12) && Run12EmbedJewel) 
	fileLoc =  Form("/phenix/crs/phnxreco/lajoie/SimOut/Run%i_CUAU%s_%i_%i/%i", runno, HIAdd.Data(), (int)centlow[ic], (int)centhigh[ic], cmsE);
      if((runno==12) && !Run12EmbedJewel) 
	fileLoc =  Form("/phenix/crs/phnxreco/lajoie/SimOut/Run%i_CUAU%s_%i_%i/%i", runno, HIAdd.Data(), (int)centlow[ic], (int)centhigh[ic], cmsE);
    }

    TString fNormName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_%s%s.root", fileLoc.Data(), runno, ckins[normIdx], p_or_h.Data(),RString.Data(),Suffix.Data()); 

    // special name for Jon's files
    //if(runno==15) fNormName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_%s%s_NoPbGl.root", 
    //			     fileLoc.Data(), runno, ckins[normIdx], p_or_h.Data(),RString.Data(),kSuffix.Data()); 

    if(embed) {
      fNormName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_embed_%s_%i_%i%s.root", fileLoc.Data(), runno, ckins[normIdx], p_or_h.Data(),RString.Data(),
			       (int)centlow[ic], (int)centhigh[ic],Suffix.Data()); 

      //if(runno==15) fNormName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_embed_%s_%i_%i%s_NoPbGl.root", 
      //			       fileLoc.Data(), runno, ckins[normIdx], p_or_h.Data(),RString.Data(),(int)centlow[ic],(int)centhigh[ic],kSuffix.Data());
    }

    TFile *fNorm = new TFile(fNormName);
    TH1I *hNormCkin = (TH1I *)fNorm->Get("hEvents"); 
    double nEventsNorm = 0.0; 
    if(hNormCkin){
      nEventsNorm = (double) hNormCkin->GetBinContent(1); 
    }
    else{
      cout << " ERROR - null normalization histogram!" << endl; 
      return;
    }
    fNorm->Close(); 
    delete fNorm; 
  
    for(unsigned int index = 0; index<nCkinRuns; index++){
 
      TString fInName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_%s%s.root", fileLoc.Data(), runno, ckins[index], p_or_h.Data(),RString.Data(),Suffix.Data()); 

      // special name for Jon's files
      //if(runno==15) fInName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_%s%s_NoPbGl.root", 
      //			     fileLoc.Data(), runno, ckins[index], p_or_h.Data(),RString.Data(),kSuffix.Data()); 

      if(embed) {
	fInName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_embed_%s_%i_%i%s.root", fileLoc.Data(), runno, ckins[index], p_or_h.Data(),RString.Data(),
			       (int)centlow[ic], (int)centhigh[ic],Suffix.Data()); 

	//if(runno==15) fInName = Form("%s/CHAINS/simChain_Run%i_CKIN%i%s_embed_%s_%i_%i%s_NoPbGl.root", 
	//		       fileLoc.Data(), runno, ckins[index], p_or_h.Data(),RString.Data(),(int)centlow[ic],(int)centhigh[ic],kSuffix.Data()); 
      }

      TFile *fIn = new TFile(fInName);

      double nEventsThisCkin = 0.0; 
      TH1I *hNormThisCkin = (TH1I *)fIn->Get("hEvents"); 
      if(hNormThisCkin){
	nEventsThisCkin = (double) hNormThisCkin->GetBinContent(1); 
      }
      else{
	cout << " ERROR - null events histogram!" << endl; 
	return;
      }

      // everything scaled relative to one CKIN range
      // equal number of events in each centrality range
      float csWeight = (nEventsNorm*crossSec[index])/(nEventsThisCkin*crossSec[normIdx])
	*(centWeightFactor[0]/centWeightFactor[ic]);

      // account for the wider centrality range in the peripheral bins 
      if( (runno==8) && (nMBCombine==4) && (ic==3) ) 
	csWeight *= (88.0-60.0)/20.0; 

      if( (runno==12) && (nMBCombine==4) && (ic==3) ) 
	csWeight *= (93.0-60.0)/20.0; 

      if( (runno==15) && (nMBCombine==4) && (ic==3) ) 
	csWeight *= (84.0-60.0)/20.0; 
 
      float simLumi = nEventsThisCkin/crossSec[index]; 
      hCSWeight->SetBinContent(index+1,csWeight); 
      hSimLumi->SetBinContent(index+1,simLumi); 

      totalCs += crossSec[index]; 
      csWeightedEvents += nEventsThisCkin*csWeight; 
    
      // picoJet TTRees:
      TTree *tMatch = (TTree *)fIn->Get("pico_mJets"); 
      TTree *tTrue = (TTree *)fIn->Get("pico_tJets"); 
      TTree *tReco = (TTree *)fIn->Get("pico_rJets"); 
   
      if(!tMatch || !tTrue || !tReco){
	cout << " Error reading in tMatch, tTrue or TReco trees...." << endl; 
	return; 
      }

      cout << " Processing CKIN = " << ckins[index] << endl;

      LinkMatchedVars(tReco, tMatch, tTrue);

      // MATCHED AND TRUTH JETS

      unsigned long long nEntries = tMatch->GetEntries();

      cout << "Examining " << nEntries << " events (jet events) ...." << endl; 

      for(unsigned long long iEnt = 0; iEnt < nEntries; iEnt++){

	if( (halfstats==1) && (iEnt%2 == 0)){
	  // skip half the statistics (even numbers)
	  continue; 
	}

	tReco->GetEntry(iEnt);
	tMatch->GetEntry(iEnt);
	tTrue->GetEntry(iEnt);

	hEmbed->Fill((double)embed); 
      
	// Set the event weight - 
	// If NLO weighting is used then we need to decide how to 
	// do that. We use the highest energy truth jet or Q^2

	float maxTruthPT = 0.0; 
	int max_truth_idx = -1; 
	for(int iJet = 0; iJet < nTruthJets; iJet++){
	
	  if(t_pT[iJet] <= 0) continue;
	  //if(t_nc[iJet] < 3) continue; // constituents cut
	  if(fabs(t_eta[iJet]) > 0.15) continue;

	  // do not consider jets with unphysical pT
	  if(t_pT[iJet]>=PTBINS_TRUE[NPTBINS_TRUE]) continue; 
	  if(t_pT[iJet]<PTBINS_TRUE[0]) continue; 

	  int tjbin = GetBin(t_pT[iJet], NPTBINS_TRUE, PTBINS_TRUE);
	  if(tjbin<0) continue; 
	
	  if(t_pT[iJet]>maxTruthPT){
	    maxTruthPT = t_pT[iJet]; 
	    max_truth_idx = iJet; 
	  }

	}

	// Implement tpT cutoffs to avoid overweighted statisticsl fluctuations
	if(max_truth_idx>=0){
	  if(t_pT[max_truth_idx]>tpT_cutoff[index]) continue; 
	}

	// Factor in the NLO weight based on the true jet pT
	// If no qualifying jet found, use the ckin of the event 
	double NLOWeightFactor = 1.0; 
	double NLO_pt = sqrt(evt_Qsqr); 
	if(max_truth_idx>=0){ 
	  NLO_pt = t_pT[max_truth_idx]; 
	}

	if((runno==12) && (weightNLO==1) && (embed==0)){
	  // limit range for 200 GeV p+p NLO
	  if(NLO_pt<5.0) 
	    NLOWeightFactor = NLOWeight->Eval(5.0); 
	  else if (NLO_pt>60.0)
	    NLOWeightFactor = NLOWeight->Eval(60.0); 
	  else
	    NLOWeightFactor = NLOWeight->Eval(NLO_pt); 
	}
	else{
	  NLOWeightFactor = NLOWeight->Eval(NLO_pt); 
	}
 
	float evWeight = hVtxRatio->Interpolate(vertex)*csWeight*NLOWeightFactor; 
	fullyWeightedEvents += evWeight; 

	// Simulated random values:
	int r_spinPat = 0; 
	if((runno==13)||((runno==15)&&(embed==0))) r_spinPat = myRand->Integer(4); // pp 0-3, long or transverse
	if(((runno==15)&&(embed==1)))r_spinPat = myRand->Integer(2) + 4; // pAu 4-5

	// must use value from file to get trig eff. correct
	int even_odd =  ip12_clock_cross%2; 
	if((even_odd<0)||(even_odd>1)) even_odd = myRand->Integer(2); 
 
	hVertex->Fill(vertex,evWeight);
        hTrue_Jets->Fill(t_pT[max_truth_idx],evWeight); 

	// --

	if(fabs(vertex)>vtxCut) continue;

	fullyWeightedCutEvents += evWeight; 
	nEvenOdd[even_odd]+=evWeight; 
	nSpinPat[even_odd][r_spinPat]+=evWeight; 
            
	hVertexCut->Fill(vertex,evWeight);

	if(nRecoJets > 0){

	  hVertexJets->Fill(vertex,evWeight);

	  int indexMax = -1;
	  if(useML)
	    indexMax = GetMaxPtIndex(r_ml_pT, r_cf, r_phi, nRecoJets, runno, NPTBINS, PTBINS, AcceptFlag);
	  else
	    indexMax = GetMaxPtIndex(r_pT, r_cf, r_phi, nRecoJets, runno, NPTBINS, PTBINS, AcceptFlag);

	  // do not consider jets with unphysical reco pT
	  if(indexMax>=0){
	    if(useML){
	      if(r_ml_pT[indexMax]>=PTBINS[NPTBINS]) continue; 
	      if(r_ml_pT[indexMax]<PTBINS[0]) continue; 
	    }
	    else{
	      if(r_pT[indexMax]>rpT_cutoff[index]) continue;
	      if(r_pT[indexMax]>=PTBINS[NPTBINS]) continue; 
	      if(r_pT[indexMax]<PTBINS[0]) continue; 
	    }
	  }
	  else{
	    continue; 
	  }

	  float r_maxPt = -9999.0; 
	  if(useML)
	    r_maxPt = r_ml_pT[indexMax];
	  else
	    r_maxPt = r_pT[indexMax];

	  if(fabs(r_eta[indexMax])>0.15) continue; 

	  int binReco = GetBin(r_maxPt, NPTBINS, PTBINS);
	  if(binReco<0) continue; 

	  hReco->Fill(r_maxPt,evWeight);
	  hCF->Fill(r_maxPt, r_cf[indexMax],evWeight); 
	  hPhi->Fill(r_phi[indexMax],evWeight); 
	  hNC->Fill(r_maxPt, r_nc[indexMax],evWeight); 

	  // NoUE histograms are filled without trigger requirement
	  // This increases statistics for the NoUE reference
	  // Unfold using NoTriggerEff 
	  float r_z_g_noUE = 0.0; // Used for trig also
	  if(embed){
	    r_z_g_noUE = SubstructureNoUE(indexMax,R);
	    if(r_z_g_noUE>0.0) {
	      hZgNoUE->Fill(r_maxPt,r_z_g_noUE,evWeight); 
	      hOangNoUE->Fill(r_maxPt,r_oang_noUE,evWeight); 
	    }
	    
	    for(unsigned int i=0; i<re_cs_z->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue; 
	      if(re_cs_UE->at(indexMax).at(i)==0) {
		hFFZNoUE->Fill(r_maxPt, re_cs_z->at(indexMax).at(i),evWeight);
		hFFXINoUE->Fill(r_maxPt, -log(re_cs_z->at(indexMax).at(i)),evWeight);
	      }
	    }

	    for(unsigned int i=0; i<re_cs_jT->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
	      //if(re_cs_UE->at(indexMax).at(i)==0) hFFJTNoUE->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i),evWeight);
	      if(re_cs_UE->at(indexMax).at(i)==0) hFFJTNoUE->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i)/r_maxPt,evWeight);
	    }

	    for(unsigned int i=0; i<re_cs_dR->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
	      if(re_cs_UE->at(indexMax).at(i)==0) hFFdRNoUE->Fill(r_maxPt, re_cs_dR->at(indexMax).at(i),evWeight);
	    }
	  }

	  hReco_bin[binReco]->Fill(r_maxPt,evWeight);

	  bool trigFired = false;
	  if((arm0FiredTrigger && (getArm(r_phi[indexMax])==0)) || (arm1FiredTrigger && (getArm(r_phi[indexMax])==1))) trigFired = true;  
	  if((runno==12) && (embed>0)) trigFired = true; 

	  if(trigFired){

	    hReco_ERTfired->Fill(r_maxPt,evWeight);

	    double jetCharge = JetCharge(indexMax); 
	    hRecoJetCharge_ERTfired->Fill(r_maxPt,jetCharge);

	    hReco_bin_ERTfired[binReco]->Fill(r_maxPt,evWeight);
 	
	    // get the index of the matched truth jet
	  
	    int matched_truth_idx = -1; 
	    for(int ij2 = 0; ij2<nMatchedJets; ij2++){
	      if(m_r_idx[ij2]==indexMax){
		matched_truth_idx = m_t_idx[ij2];
		hTrue_ERTfired->Fill(t_pT[m_t_idx[ij2]],evWeight); 
		break; 
	      }
	    }

	    // Triggered NoUE histograms 
	    // Unfold using standard matrix 
	    if(embed){
	      if(r_z_g_noUE>0.0) {
		hZgNoUETrig->Fill(r_maxPt,r_z_g_noUE,evWeight);
		hOangNoUETrig->Fill(r_maxPt,r_oang_noUE,evWeight);
	      }
	    
	      for(unsigned int i=0; i<re_cs_z->at(indexMax).size(); i++){
		if(re_cs_charge->at(indexMax).at(i)==0.0) continue; 
		if(re_cs_UE->at(indexMax).at(i)==0) {
		  hFFZNoUETrig->Fill(r_maxPt, re_cs_z->at(indexMax).at(i),evWeight);
		  hFFXINoUETrig->Fill(r_maxPt, -log(re_cs_z->at(indexMax).at(i)),evWeight);
		}
	      }

	      for(unsigned int i=0; i<re_cs_jT->at(indexMax).size(); i++){
		if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
		if(re_cs_UE->at(indexMax).at(i)==0) hFFJTNoUETrig->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i)/r_maxPt,evWeight);
	      }

	      for(unsigned int i=0; i<re_cs_dR->at(indexMax).size(); i++){
		if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
		if(re_cs_UE->at(indexMax).at(i)==0) hFFdRNoUETrig->Fill(r_maxPt, re_cs_dR->at(indexMax).at(i),evWeight);
	      }
	    }

	    if(r_zg[indexMax]>=0.1 && r_zg[indexMax]<=0.5) {
	      hZg->Fill(r_maxPt,r_zg[indexMax],evWeight); 
	      hOang->Fill(r_maxPt,r_oang[indexMax],evWeight); 
	    }

	    for(unsigned int i=0; i<re_cs_z->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue; 
	      hFFZ->Fill(r_maxPt, re_cs_z->at(indexMax).at(i),evWeight); 
	      hFFXI->Fill(r_maxPt, -log(re_cs_z->at(indexMax).at(i)),evWeight); 
	    }

	    for(unsigned int i=0; i<re_cs_jT->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
	      hFFJT->Fill(r_maxPt, re_cs_jT->at(indexMax).at(i)/r_maxPt,evWeight); 
	    }

	    for(unsigned int i=0; i<re_cs_dR->at(indexMax).size(); i++){
	      if(re_cs_charge->at(indexMax).at(i)==0.0) continue;
	      hFFdR->Fill(r_maxPt, re_cs_dR->at(indexMax).at(i),evWeight); 
	    }

	    if(r_spinPat < 6){

	      if(((runno==12)||(runno==15))){

                // ------------------------------------------------------------------
                // EEC calculation over max pt reco jets
                // ------------------------------------------------------------------
                /* N.B. cuts on jet pt, eta, and CNF are baked into the requirement
                 *   that indexMax >= 0. A max reco jet has to satisfy these cuts.
                 */
                if (doRecoEEC) {

                  // collect jet and spin information into a handy struct
                  PHEC::Type::Jet jet_reco(
                    r_cf[indexMax],
                    r_pT[indexMax],
                    r_eta[indexMax],
                    r_phi[indexMax],
                    jetCharge,
                    r_spinPat
                  );

                  // loop through pairs of constituents
                  for (
                    std::size_t iRecoCstA = 0;
                    iRecoCstA < re_cs_z->at(indexMax).size();
                    ++iRecoCstA
                  ) {

		    // keep only charged cst.s
                    if (doRecoEECChargedOnly && (re_cs_charge->at(indexMax).at(iRecoCstA) == 0.0)) {
                      continue;
                    }

                    for (
                      std::size_t iRecoCstB = 0;
                      iRecoCstB <= iRecoCstA;
                      ++iRecoCstB
                    ) {

		      // keep only charged cst.s
                      if (doRecoEECChargedOnly && (re_cs_charge->at(indexMax).at(iRecoCstB) == 0.0)) {
                         continue;
                      }

                      // collect cst information into a handy struct
                      PHEC::Type::Cst cstA_reco(
                        re_cs_z->at(indexMax).at(iRecoCstA),
                        re_cs_jT->at(indexMax).at(iRecoCstA),
                        re_cs_eta->at(indexMax).at(iRecoCstA),
                        re_cs_phi->at(indexMax).at(iRecoCstA),
                        re_cs_charge->at(indexMax).at(iRecoCstA)
                      );
                      PHEC::Type::Cst cstB_reco(
                        re_cs_z->at(indexMax).at(iRecoCstB),
                        re_cs_jT->at(indexMax).at(iRecoCstB),
                        re_cs_eta->at(indexMax).at(iRecoCstB),
                        re_cs_phi->at(indexMax).at(iRecoCstB),
                        re_cs_charge->at(indexMax).at(iRecoCstB)
                      );

                      // run 2-point calculation for pair
                      recoEEC.CalcEEC( jet_reco, std::make_pair(cstA_reco, cstB_reco), evWeight );

                    }  // end 2nd cst loop
                  }  // end 1st cst loop
                }  // end max reco jet eec calculation

                // ------------------------------------------------------------------


                // ----------------------------------------------------------------------
                // EEC calculation over max pt truth jets
                // ----------------------------------------------------------------------
                /* N.B. cuts on jet pt and eta are baked into the requirement
                 *   that matched_truth_idx >= 0. A max truth jet has to satisfy
                 *   these cuts.
                 */
                if (doTrueEEC && (matched_truth_idx>=0)) {

                  // calculate truth jet charge
                  const double jetChargeTrue = JetChargeTrue(matched_truth_idx);

                  // collect jet and spin information into a handy struct
                  PHEC::Type::Jet jet_true(
		    t_cf[matched_truth_idx],
                    t_pT[matched_truth_idx],
                    t_eta[matched_truth_idx],
                    t_phi[matched_truth_idx],
                    jetChargeTrue,
                    r_spinPat
                  );

                  // loop through pairs of constituents
                  for (
                    std::size_t iTruthCstA = 0;
                    iTruthCstA < tr_cs_z->at(matched_truth_idx).size();
                    ++iTruthCstA
                  ) {

		    // keep only charged cst.s
		    if (doTrueEECChargedOnly && (tr_cs_charge->at(matched_truth_idx).at(iTruthCstA) == 0.0)) {
		      continue;
		    }

                    for (
                      std::size_t iTruthCstB = 0;
                      iTruthCstB <= iTruthCstA;
                      ++iTruthCstB
                    ) {

		      // keep only charged cst.s
		      if (doTrueEECChargedOnly && (tr_cs_charge->at(matched_truth_idx).at(iTruthCstB) == 0.0)) {
			continue;
		      }

                      // collect cst information into a handy struct
                      PHEC::Type::Cst cstA_true(
                        tr_cs_z->at(matched_truth_idx).at(iTruthCstA),
                        tr_cs_jT->at(matched_truth_idx).at(iTruthCstA),
                        tr_cs_eta->at(matched_truth_idx).at(iTruthCstA),
                        tr_cs_phi->at(matched_truth_idx).at(iTruthCstA),
                        tr_cs_charge->at(matched_truth_idx).at(iTruthCstA)
                      );
                      PHEC::Type::Cst cstB_true(
                        tr_cs_z->at(matched_truth_idx).at(iTruthCstB),
                        tr_cs_jT->at(matched_truth_idx).at(iTruthCstB),
                        tr_cs_eta->at(matched_truth_idx).at(iTruthCstB),
                        tr_cs_phi->at(matched_truth_idx).at(iTruthCstB),
                        tr_cs_charge->at(matched_truth_idx).at(iTruthCstB)
                      );

                      // run 2-point calculation for pair
                      trueEEC.CalcEEC( jet_true, std::make_pair(cstA_true, cstB_true), evWeight );

                    }  // end 2nd cst loop
                  }  // end 1st cst loop
                }  // end max truth jet eec calculation

                // ------------------------------------------------------------------

	      }

	    }
 
	  }

	}

      }

      fIn->Close();
      delete fIn; 

    }

  }
  

  // --------------------------------------------------------------------------
  // EEC calculations
  // --------------------------------------------------------------------------

  if (doTrueEEC) trueEEC.End( fOut );
  if (doRecoEEC) recoEEC.End( fOut );

  // --------------------------------------------------------------------------

  fOut->Write(); 
  fOut->Close(); 

  delete myRand; 

  cout << endl; 
  cout << "Total cross section = " << totalCs << " mb" << endl ; 
  cout << "Total cross section weighted events = " << csWeightedEvents << endl; 
  cout << "Total fully weighted events = " << fullyWeightedEvents << endl; 
  cout << "Total fully weighted events after vertex cut = " << fullyWeightedCutEvents << endl; 
  cout << "Effective integrated luminosity = " << fullyWeightedCutEvents/totalCs << " mb^-1" << endl; 
  cout << endl; 

  return;

}

void spinMC(int nIter = 100, int weightALL = 0)
{

  const int nBins = NPTBINS_RECO_13; 

  TFile fOut("spinMC.root","recreate"); 

  double EODiff[nBins]; 
  double EODiffS[nBins]; 
  double ALLVal[nBins];
  double ALLValS[nBins];

  double tEODiff[nBins]; 
  double tEODiffS[nBins]; 
  double tALLVal[nBins];
  double tALLValS[nBins];

  double TruthALL[nBins]; 
  double TruthALLS[nBins]; 
  
  double RawALL[nBins];
  double RawALLS[nBins];

  double R_even = 1.0; 
  double R_odd = 1.0; 
  double R = 1.0; 
  
  TTree *tSpinMC  = new TTree("tSpinMC", "A_LL Spin MC");
  tSpinMC -> Branch("weighted",&weightALL,"weighted/I");
  tSpinMC -> Branch("R_even",&R_even,"R_even/D");
  tSpinMC -> Branch("R_odd",&R_odd,"R_odd/D");
  tSpinMC -> Branch("R_comb",&R,"R/D");

  for(int i=0; i<nBins; i++){
    tSpinMC -> Branch(Form("EODiff%i",i),&EODiff[i],Form("EODiff%i/D",i));
    tSpinMC -> Branch(Form("EODiffS%i",i),&EODiffS[i],Form("EODiffS%i/D",i));
    tSpinMC -> Branch(Form("ALLVal%i",i),&ALLVal[i],Form("ALLVal%i/D",i));
    tSpinMC -> Branch(Form("ALLValS%i",i),&ALLValS[i],Form("ALLValS%i/D",i));
    tSpinMC -> Branch(Form("tEODiff%i",i),&tEODiff[i],Form("tEODiff%i/D",i));
    tSpinMC -> Branch(Form("tEODiffS%i",i),&tEODiffS[i],Form("tEODiffS%i/D",i));
    tSpinMC -> Branch(Form("tALLVal%i",i),&tALLVal[i],Form("tALLVal%i/D",i));
    tSpinMC -> Branch(Form("tALLValS%i",i),&tALLValS[i],Form("tALLValS%i/D",i));
    tSpinMC -> Branch(Form("TruthALL%i",i),&TruthALL[i],Form("TruthALL%i/D",i));
    tSpinMC -> Branch(Form("TruthALLS%i",i),&TruthALLS[i],Form("TruthALLS%i/D",i));
    tSpinMC -> Branch(Form("RawALL%i",i),&RawALL[i],Form("RawALL%i/D",i));
    tSpinMC -> Branch(Form("RawALLS%i",i),&RawALLS[i],Form("RawALLS%i/D",i));
  }

  // simulate, unfold, calculate new A_LL

  for(int i=0; i<nIter; i++){

    R_even = 1.0; 
    R_odd = 1.0; 
    R = 1.0; 

    // run jetAnaSim
    jetAnaSim(13, 0.3, 0, 0.0, 20.0, &R_even, &R_odd, &R, weightALL, 1); 
    
    // unfold w/o trigger efficiency correction
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,0,0,-1,1,2,2\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,0,1,-1,1,2,2\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,1,0,-1,1,2,2\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,1,1,-1,1,2,2\\)"); 

    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,0,-1,-1,1,2,2\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,1,-1,-1,1,2,2\\)"); 

    // calculate the A_LL's
    gSystem->Exec(Form("root -l -b -q calcALL.C\\(-1,1,0,%f,%f\\)",R_even,R_odd)); 
    gSystem->Exec(Form("root -l -b -q calcALLCombMatrix.C\\(-1,1,0,%f\\)",R)); 

    //Open the output TGraphs, 
    //histogram the even-odd difference and 
    //the combined A_LL. 

    TFile finEvenOdd("ALL_Unfold_closure.root","read");
    TFile finComb("ALL_Unfold_Comb_closure.root","read");

    TGraphAsymmErrors *ALL_even = (TGraphAsymmErrors *) finEvenOdd.Get("gALL_even"); 
    TGraphAsymmErrors *ALL_odd  = (TGraphAsymmErrors *) finEvenOdd.Get("gALL_odd"); 
    TGraphAsymmErrors *ALL = (TGraphAsymmErrors *) finComb.Get("ALL"); 
    
    for(int j=0; j<nBins; j++) {

      double all_even, all_odd, all_comb; 
      double xT; 
      
      double all_even_high_err; 
      double all_odd_high_err; 
      double all_comb_high_err; 

      ALL_even->GetPoint(j,xT,all_even); 
      all_even_high_err = ALL_even->GetErrorYhigh(j); 
      ALL_odd->GetPoint(j,xT,all_odd); 
      all_odd_high_err = ALL_odd->GetErrorYhigh(j); 
      ALL->GetPoint(j,xT,all_comb); 
      all_comb_high_err = ALL->GetErrorYhigh(j); 

      EODiff[j] = (all_even-all_odd); 
      EODiffS[j] = (all_even-all_odd)/sqrt(pow(all_even_high_err,2)+pow(all_odd_high_err,2)); 
      ALLVal[j] = all_comb; 
      ALLValS[j] = (all_comb/all_comb_high_err); 

    }

    finEvenOdd.Close(); 
    finComb.Close();

    // unfold with trigger efficiency correction
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,0,0,-1,1,2,1\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,0,1,-1,1,2,1\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,1,0,-1,1,2,1\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,1,1,-1,1,2,1\\)"); 

    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,0,-1,-1,1,2,1\\)"); 
    gSystem->Exec("root -l -b -q doUnfolding.C\\(13,\\'p\\',0.3,0,0.0,20.0,1,-1,-1,1,2,1\\)"); 

    // calculate the A_LL's
    gSystem->Exec(Form("root -l -b -q calcALL.C\\(-1,1,0,%f,%f\\)",R_even,R_odd)); 
    gSystem->Exec(Form("root -l -b -q calcALLCombMatrix.C\\(-1,1,0,%f\\)",R)); 

    //Open the output TGraphs, 
    //histogram the even-odd difference and 
    //the combined A_LL. 

    TFile tfinEvenOdd("ALL_Unfold_closure.root","read");
    TFile tfinComb("ALL_Unfold_Comb_closure.root","read");

    TGraphAsymmErrors *tALL_even = (TGraphAsymmErrors *) tfinEvenOdd.Get("gALL_even"); 
    TGraphAsymmErrors *tALL_odd  = (TGraphAsymmErrors *) tfinEvenOdd.Get("gALL_odd"); 
    TGraphAsymmErrors *tALL = (TGraphAsymmErrors *) tfinComb.Get("ALL"); 
    
    for(int j=0; j<nBins; j++) {

      double all_even, all_odd, all_comb; 
      double xT; 
      
      double all_even_high_err; 
      double all_odd_high_err; 
      double all_comb_high_err; 

      tALL_even->GetPoint(j,xT,all_even); 
      all_even_high_err = tALL_even->GetErrorYhigh(j); 
      tALL_odd->GetPoint(j,xT,all_odd); 
      all_odd_high_err = tALL_odd->GetErrorYhigh(j); 
      tALL->GetPoint(j,xT,all_comb); 
      all_comb_high_err = tALL->GetErrorYhigh(j); 

      tEODiff[j] = (all_even-all_odd); 
      tEODiffS[j] = (all_even-all_odd)/sqrt(pow(all_even_high_err,2)+pow(all_odd_high_err,2)); 
      tALLVal[j] = all_comb; 
      tALLValS[j] = (all_comb/all_comb_high_err); 

    }
   
    tfinEvenOdd.Close(); 
    tfinComb.Close();
    
    // Extract A_LL from weighted truth spectra

    TFile tfinAna("jetAnaSimRun13_pp_R0.30.root","read");
    TH1F *hTrueSamePat = (TH1F *)tfinAna.Get("hTrueSamePat"); 
    TH1F *hTrueOppPat = (TH1F *)tfinAna.Get("hTrueOppPat"); 

    TH1F *numerator = (TH1F *) hTrueSamePat->Clone("numerator");
    numerator->Add(hTrueOppPat,-1.0*R); 
    TH1F *denominator = (TH1F *) hTrueSamePat->Clone("denominator"); 
    denominator->Add(hTrueOppPat,1.0*R); 
    numerator->Divide(denominator); 

    for(int j=1; j<nBins+1; j++) {
      TruthALL[j-1] = numerator->GetBinContent(j); 
      if(numerator->GetBinError(j)>0.0)
        TruthALLS[j-1] = TruthALL[j-1]/numerator->GetBinError(j); 
      else
        TruthALLS[j-1] = -9999.0; 
    }

    // Extract the raw A_LL

    TH1F *hRecoSamePat = (TH1F *)tfinAna.Get("hRecoSamePat_ERTfired"); 
    TH1F *hRecoOppPat = (TH1F *)tfinAna.Get("hRecoOppPat_ERTfired"); 

    TH1F *reco_numerator = (TH1F *) hRecoSamePat->Clone("reco_numerator");
    reco_numerator->Add(hRecoOppPat,-1.0*R); 
    TH1F *reco_denominator = (TH1F *) hRecoSamePat->Clone("reco_denominator"); 
    reco_denominator->Add(hRecoOppPat,1.0*R); 
    reco_numerator->Divide(reco_denominator); 

    for(int j=1; j<nBins+1; j++) {
      RawALL[j-1] = reco_numerator->GetBinContent(j);
      if(reco_numerator->GetBinError(j)>0.0)
        RawALLS[j-1] = TruthALL[j-1]/reco_numerator->GetBinError(j); 
      else
        RawALLS[j-1] = -9999.0; 
    }

    tfinAna.Close(); 

    // Fill the tree

    tSpinMC->Fill(); 

    // wait to avoid TRandom->SetSeed(0) getting an identical value
    gSystem->Exec("sleep 1");
 
  }

  fOut.Write(); 
  fOut.Close(); 

}
