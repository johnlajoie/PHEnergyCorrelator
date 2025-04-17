/// ============================================================================
/*! \file    PHCorrelatorCalculator.h
 *  \authors Derek Anderson, Alex Clarke
 *  \date    09.21.2024
 *
 *  Driver class to run n-point energy-energy correlator
 *  calculations on inputs.
 */
/// ============================================================================

#ifndef PHCORRELATORCALCULATOR_H
#define PHCORRELATORCALCULATOR_H

// c++ utilities
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>
// root libraries
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector3.h>
// analysis componenets
#include "PHCorrelatorAnaTools.h"
#include "PHCorrelatorAnaTypes.h"
#include "PHCorrelatorHistManager.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! ENC Calculator
  // ==========================================================================
  class Calculator {

    private:

      // data members (calc options)
      double       m_weight_power;
      Type::Weight m_weight_type;

      // data members (bins)
      std::vector< std::pair<float, float> > m_ptjet_bins;
      std::vector< std::pair<float, float> > m_cfjet_bins;
      std::vector< std::pair<float, float> > m_chrg_bins;

      // data member (hist manager)
      HistManager m_manager;

      // ---------------------------------------------------------------------=
      //! Get weight of a constituent
      // ----------------------------------------------------------------------
      double GetCstWeight(TLorentzVector& cst, TLorentzVector& jet) {

        // grab relevant cst & jet values
        double numer = 1.0;
        double denom = 1.0;
        switch (m_weight_type) {

          case Type::E:
            numer = cst.E();
            denom = jet.E();
            break;

          case Type::Et:
            numer = cst.Et();
            denom = jet.Et();
            break;

          case Type::Pt:
            numer = cst.Pt();
            denom = jet.Pt();
            break;

          default:
            numer = cst.Pt();
            denom = jet.Pt();
            break;

        }  // end switch

        // raise cst, jet values to specified value (defualt is 1.0)
        numer = pow(numer, m_weight_power);
        denom = pow(denom, m_weight_power);

        // calculate weight and exit
        const double weight = numer / denom;
        return weight;

      }  // end 'GetCstWeight(TLorentzVector&, TLorentzVector&)'

      // ----------------------------------------------------------------------
      //! Get hist index/indices
      // ----------------------------------------------------------------------
      /*! Returns a vector of indices of histograms to be filled. If not doing
       *  spin sorting, vector will only have FOUR entries, which are
       *    [0] = integrated pt, charge (+cf bin)
       *    [1] = pt bin, integrated charge (+cf bin)
       *    [2] = charge bin, integrated pt (+cf bin)
       *    [3] = pt bin, charge bin (+cf bin)
       *
       *  If doing spin sorting, the vector will have either 4, 8, or 16
       *  entries.
       *    - `size() == 16`: pp case; entries correspond to spin-integrated,
       *       blue-only, yellow-only, and blue-and-yellow indices.
       *    - `size() == 8`: pAu case; entries correspond to spin-integrated
       *       and blue-only indices.
       *    - `size() == 4`: an unexpected spin pattern was provided, only
       *       spin-integrated case returned.
       *
       *  The order of the vector of indices will always be
       *    [0 - 3]   = spin integrated
       *    [4 - 7]   = blue beam index
       *    [8 - 11]  = yellow beam index
       *    [12 - 15] = blue and yellow index
       */
      std::vector<Type::HistIndex> GetHistIndices(const Type::Jet& jet) {

        // for pt and cf, index will correspond to what bin
        // the jet falls in PLUS integrated pt/charge bin
        Type::HistIndex base_index(0, 0, 0, 0);

        // determine pt bin
        if (m_manager.GetDoPtJetBins()) {
          for (std::size_t ipt = 0; ipt < m_ptjet_bins.size(); ++ipt) {
            if ((jet.pt >= m_ptjet_bins[ipt].first) && (jet.pt < m_ptjet_bins[ipt].second)) {
              base_index.pt = ipt;
            }
          }  // end pt bin loop
        }

        // determine cf bin
        //   - n.b. there is NO integrated bin for cf
        if (m_manager.GetDoCFJetBins()) {
          for (std::size_t icf = 0; icf < m_cfjet_bins.size(); ++icf) {
            if ((jet.cf >= m_cfjet_bins[icf].first) && (jet.cf < m_cfjet_bins[icf].second)) {
              base_index.cf = icf;
            }
          }  // end charge bin loop
        }

        // determine charge bin
        if (m_manager.GetDoChargeBins()) {
          for (std::size_t ich = 0; ich < m_chrg_bins.size(); ++ich) {
            if ((jet.charge >= m_chrg_bins[ich].first) && (jet.charge < m_chrg_bins[ich].second)) {
              base_index.chrg = ich;
            }
          }  // end charge bin loop
        }

        // By default, only add spin-integrated bin
        std::vector<std::size_t> spin_indices;
        spin_indices.push_back( HistManager::Int );

        // if needed, determine spin bins
        if (m_manager.GetDoSpinBins()) {
            switch (jet.pattern) {

              // blue up, yellow up (pp)
              case Type::PPBUYU:
                spin_indices.push_back( HistManager::BU );
                spin_indices.push_back( HistManager::YU );
                spin_indices.push_back( HistManager::BUYU );
                break;

              // blue down, yellow up (pp)
              case Type::PPBDYU:
                spin_indices.push_back( HistManager::BD );
                spin_indices.push_back( HistManager::YU );
                spin_indices.push_back( HistManager::BDYU );
                break;

              // blue up, yellow down (pp)
              case Type::PPBUYD:
                spin_indices.push_back( HistManager::BU );
                spin_indices.push_back( HistManager::YD );
                spin_indices.push_back( HistManager::BUYD );
                break;

              // blue down, yellow down (pp)
              case Type::PPBDYD:
                spin_indices.push_back( HistManager::BD );
                spin_indices.push_back( HistManager::YD );
                spin_indices.push_back( HistManager::BDYD );
                break;

              // blue up (pAu)
              case Type::PABU:
                spin_indices.push_back( HistManager::BU );
                break;

              // blue down (pAu)
              case Type::PABD:
                spin_indices.push_back( HistManager::BD );
                break;

              // by default, only add integrated
              default:
                break;

          }
        }

        // now assemble list of indices to fill
        std::vector<Type::HistIndex> indices;
        for (std::size_t isp = 0; isp < spin_indices.size(); ++isp) {

          // integrated everything (except cf)
          indices.push_back(
            Type::HistIndex(
              m_ptjet_bins.size(),
              base_index.cf,
              m_chrg_bins.size(),
              spin_indices[isp]
            )
          );

          // binned pt, integrated charge
          indices.push_back(
            Type::HistIndex(
              base_index.pt,
              base_index.cf,
              m_chrg_bins.size(),
              spin_indices[isp]
            )
          );

          // binned charge, integrated pt
          indices.push_back(
            Type::HistIndex(
              m_ptjet_bins.size(),
              base_index.cf,
              base_index.chrg,
              spin_indices[isp]
            )
          );

          // binned everything
          indices.push_back(
            Type::HistIndex(
              base_index.pt,
              base_index.cf,
              base_index.chrg,
              spin_indices[isp]
            )
          );
        }
        return indices;

      }  // end 'GetHistIndices(Type::Jet&)'

    public:

      /* TODO
       *   - DoE3CCalc(jet, {cst, cst, cst})
       *   - DoLECCalc(jet, {lambda, cst})
       *   - FindLambda({cst...})
       */

      // ----------------------------------------------------------------------
      //! Getters
      // ----------------------------------------------------------------------
      HistManager& GetManager() {return m_manager;}

      // ----------------------------------------------------------------------
      //! Setters
      // ----------------------------------------------------------------------
      void SetWeightPower(const double power)       {m_weight_power = power;}
      void SetWeightType(const Type::Weight weight) {m_weight_type  = weight;}
      void SetHistTag(const std::string& tag)       {m_manager.SetHistTag(tag);}

      // ----------------------------------------------------------------------
      //! Set jet pt bins
      // ----------------------------------------------------------------------
      void SetPtJetBins(const std::vector< std::pair<float, float> >& bins) {

        // copy bins to member
        m_ptjet_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_ptjet_bins.begin());

        // update hist manager and exit
        m_manager.DoPtJetBins( m_ptjet_bins.size() );
        return;

      }  // end 'SetPtJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Set jet CF bins
      // ----------------------------------------------------------------------
      void SetCFJetBins(const std::vector< std::pair<float, float> >& bins) {

        // copy bins to member
        m_cfjet_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_cfjet_bins.begin());

        // update hist manager and exit
        m_manager.DoCFJetBins( m_cfjet_bins.size() );
        return;

      }  // end 'SetCFJetBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Set jet charge bins
      // ----------------------------------------------------------------------
      void SetChargeBins(const std::vector< std::pair<float, float> >& bins) {

        // copy bins to member
        m_chrg_bins.resize( bins.size() );
        std::copy(bins.begin(), bins.end(), m_chrg_bins.begin());

        // update hist manager and exit
        m_manager.DoChargeBins( m_chrg_bins.size() );
        return;

      }  // end 'SetChargeBins(std::vector<std::pair<float, float>>&)'

      // ----------------------------------------------------------------------
      //! Turn on/off spin binning
      // ----------------------------------------------------------------------
      void SetDoSpinBins(const bool spin) {

        // and set corresponding flag in manager
        m_manager.DoSpinBins(spin);
        return;

      }  // end 'SetDoSpinBins(bool)'

      // ----------------------------------------------------------------------
      //! Initialize calculator
      // ----------------------------------------------------------------------
      void Init(const bool do_eec, const bool do_e3c = false, const bool do_lec = false) {

        // turn on/off relevant histograms
        m_manager.SetDoEECHists(do_eec);
        m_manager.SetDoE3CHists(do_e3c);
        m_manager.SetDoLECHists(do_lec);

        // then generate necessary histograms
        m_manager.GenerateHists();
        return;

      } // end 'Init(bool, bool, bool)'

      // ----------------------------------------------------------------------
      //! Do EEC calculation
      // ----------------------------------------------------------------------
      /*! Note that the optional third argument, `evtweight`, is there
       *  to allow for weighting by ckin, spin, etc. By default, it's
       *  set to 1.
       */ 
      void CalcEEC(
        const Type::Jet& jet,
        const std::pair<Type::Cst, Type::Cst>& csts,
        const double evt_weight = 1.0
      ) {

        // calculate jet, cst kinematics --------------------------------------

        // get jet 4-momenta
        TLorentzVector vecJet4  = Tools::GetJetLorentz(jet, false);
        TLorentzVector unitJet4 = Tools::GetJetLorentz(jet, true);

        // get cst 4-momenta
        std::pair<TLorentzVector, TLorentzVector> vecCst4 = std::make_pair(
          Tools::GetCstLorentz(csts.first, jet.pt, false),
          Tools::GetCstLorentz(csts.second, jet.pt, false)
        );

       
	/*

	// Collins/Boer-Mulders Analysis
	  
        // get average of cst 3-vectors
        TVector3 vecAvgCst3 = Tools::GetWeightedAvgVector(
          vecCst4.first.Vect(),
          vecCst4.second.Vect(),
          false
        );
        TVector3 unitAvgCst3 = Tools::GetWeightedAvgVector(
          vecCst4.first.Vect(),
          vecCst4.second.Vect(),
          true
        );
	*/

	// Dihadron FF Analysis

        // (0) get beam and spin directions
        //   first  = blue beam/spin
        //   second = yellow beam/spin
        std::pair<TVector3, TVector3> vecBeam3 = Tools::GetBeams();
        std::pair<TVector3, TVector3> vecSpin3 = Tools::GetSpins( jet.pattern );

	// Define the vectors for the angle calculations

	TVector3 PC = vecCst4.first.Vect() + vecCst4.second.Vect();
	TVector3 PC_unit = PC.Unit(); 
	TVector3 RC = 0.5*(vecCst4.first.Vect() - vecCst4.second.Vect());
	
	// blue beam is PB, yellow is PA

	TVector3 PB = vecBeam3.first;
	TVector3 PB_unit = PB.Unit();
	TVector3 SB = vecSpin3.first; 

	TVector3 PA = vecBeam3.second;
	TVector3 PA_unit = PA.Unit(); 
	TVector3 SA = vecSpin3.second; 
	
	// Blue Polarized

	double cThetaSB = (PB_unit.Cross(PC)*(1.0/(PB_unit.Cross(PC).Mag()))).Dot(PB_unit.Cross(SB)*(1.0/PB_unit.Cross(SB).Mag())); 
	double sThetaSB = (PC.Cross(SB)).Dot(PB_unit)*(1.0/((PB_unit.Cross(PC).Mag())*(PB_unit.Cross(SB).Mag()))); 

	// Yellow Polarized

	double cThetaSA = (PA_unit.Cross(PC)*(1.0/(PA_unit.Cross(PC).Mag()))).Dot(PA_unit.Cross(SA)*(1.0/PA_unit.Cross(SA).Mag())); 
	double sThetaSA = (PC.Cross(SA)).Dot(PA_unit)*(1.0/((PA_unit.Cross(PC).Mag())*(PA_unit.Cross(SA).Mag()))); 

	// Dihadron

	double cThetaRC = (PC_unit.Cross(PA)*(1.0/(PC_unit.Cross(PA).Mag()))).Dot(PC_unit.Cross(RC)*(1.0/PC_unit.Cross(RC).Mag())); 
	double sThetaRC = (PA.Cross(RC)).Dot(PC_unit)*(1.0/((PC_unit.Cross(PA).Mag())*(PC_unit.Cross(RC).Mag()))); 

	// Convert to angles in the full range [0,2pi]

	double ThetaSB = (sThetaSB>0.0) ? acos(cThetaSB) : -acos(cThetaSB); 
	if (ThetaSB < 0)               ThetaSB += TMath::TwoPi();
	if (ThetaSB >= TMath::TwoPi()) ThetaSB -= TMath::TwoPi();

	double ThetaSA = (sThetaSA>0.0) ? acos(cThetaSA) : -acos(cThetaSA); 
	if (ThetaSA < 0)               ThetaSA += TMath::TwoPi();
	if (ThetaSA >= TMath::TwoPi()) ThetaSA -= TMath::TwoPi();

	double ThetaRC = (sThetaRC>0.0) ? acos(cThetaRC) : -acos(cThetaRC);
	if (ThetaRC < 0)               ThetaRC += TMath::TwoPi();
	if (ThetaRC >= TMath::TwoPi()) ThetaRC -= TMath::TwoPi();

	// The angle differences in the full range [0,2pi]

	double ThetaSB_RC = ThetaSB - ThetaRC; 
	if (ThetaSB_RC < 0)               ThetaSB_RC += TMath::TwoPi();
	if (ThetaSB_RC >= TMath::TwoPi()) ThetaSB_RC -= TMath::TwoPi();

	double ThetaSA_RC = ThetaSA - ThetaRC; 
	if (ThetaSA_RC < 0)               ThetaSA_RC += TMath::TwoPi();
	if (ThetaSA_RC >= TMath::TwoPi()) ThetaSA_RC -= TMath::TwoPi();
	
	/*

        // collins & boer-mulders angle calculations --------------------------

        // (0) get beam and spin directions
        //   first  = blue beam/spin
        //   second = yellow beam/spin
        std::pair<TVector3, TVector3> vecBeam3 = Tools::GetBeams();
        std::pair<TVector3, TVector3> vecSpin3 = Tools::GetSpins( jet.pattern );

        // (1) get vectors normal to the jet-beam plane
        std::pair<TVector3, TVector3> normJetBeam3 = std::make_pair(
          ( vecBeam3.first.Cross(unitJet4.Vect()) ).Unit(),
          ( vecBeam3.second.Cross(unitJet4.Vect()) ).Unit()
        );

        // (2) get phiSpin: angles between the jet-beam plane and spin
        //   - n.b. for spin pattern >= 4, the yellow spin is randomized
	//   - angle between jet plane and spin 

        // get vectors normal to the spin-beam plane
        std::pair<TVector3, TVector3> normJetSpin = std::make_pair(
          ( vecBeam3.first.Cross(vecSpin3.first) ).Unit(),
          ( vecBeam3.second.Cross(vecSpin3.second) ).Unit()
        );

        double phiSpinBlue = atan2( normJetBeam3.first.Cross(normJetSpin.first).Mag(), normJetBeam3.first.Dot(normJetSpin.first) );
        double phiSpinYell = atan2( normJetBeam3.second.Cross(normJetSpin.second).Mag(), normJetBeam3.second.Dot(normJetSpin.second) );

	// define the zero of phiSpin as the jet-beam plane and the sense of rotation 
	// the result will be in [0,2pi]

	if(normJetBeam3.first.Dot(vecSpin3.first)<0.0) phiSpinBlue = TMath::TwoPi() - phiSpinBlue; 
	if(normJetBeam3.second.Dot(vecSpin3.second)<0.0) phiSpinYell = TMath::TwoPi() - phiSpinYell; 

        // (3) get vector normal to hadron average-jet plane
        TVector3 normHadJet3 = ( unitJet4.Vect().Cross(unitAvgCst3) ).Unit();

        // (4) get phiHadron: angle between the jet-beam plane and the
        //   - angle between jet-hadron plane
        //   - constrain to range [0,2pi)
        double phiHadBlue = atan2( normJetBeam3.first.Cross(normHadJet3).Mag(), normJetBeam3.first.Dot(normHadJet3) );
        double phiHadYell = atan2( normJetBeam3.second.Cross(normHadJet3).Mag(), normJetBeam3.second.Dot(normHadJet3) );
	
	// define the zero of phiHad as the jet-beam plane and the sense of rotation
	// the result will be in [0,2pi]

	if(normJetBeam3.first.Dot(unitAvgCst3)<0.0) phiHadBlue = TMath::TwoPi() - phiHadBlue; 
	if(normJetBeam3.second.Dot(unitAvgCst3)<0.0) phiHadYell = TMath::TwoPi() - phiHadYell; 

        // (5) double phiHadron for boer-mulders,
        //   - constrain to [0, 2pi)
        double phiHadBlue2 = 2.0 * phiHadBlue;
        double phiHadYell2 = 2.0 * phiHadYell;
	if (phiHadBlue2 < 0)               phiHadBlue2 += TMath::TwoPi();
	if (phiHadBlue2 >= TMath::TwoPi()) phiHadBlue2 -= TMath::TwoPi();
        if (phiHadYell2 < 0)               phiHadYell2 += TMath::TwoPi();
	if (phiHadYell2 >= TMath::TwoPi()) phiHadYell2 -= TMath::TwoPi();

        // (6) now calculate phiColl: phiSpin - phiHadron,
        //   - constrain to [0, 2pi)
        double phiCollBlue = phiSpinBlue - phiHadBlue;
        double phiCollYell = phiSpinYell - phiHadYell;
	if (phiCollBlue < 0)               phiCollBlue += TMath::TwoPi();
	if (phiCollBlue >= TMath::TwoPi()) phiCollBlue -= TMath::TwoPi();
	if (phiCollYell < 0)               phiCollYell += TMath::TwoPi();
	if (phiCollYell >= TMath::TwoPi()) phiCollYell -= TMath::TwoPi();

        // (7) now calculate phiBoer: phiSpin - (2 * phiHadron),
        double phiBoerBlue = phiSpinBlue - phiHadBlue2;
        double phiBoerYell = phiSpinYell - phiHadYell2;

        // (8) constrain phiBoerBlue to [0, 2pi)
	if (phiBoerBlue < 0)               phiBoerBlue += TMath::TwoPi();
	if (phiBoerBlue >= TMath::TwoPi()) phiBoerBlue -= TMath::TwoPi();

        // (9) constrain to phiBoerYell to [0, 2pi)
	if (phiBoerYell < 0)               phiBoerYell += TMath::TwoPi();
	if (phiBoerYell >= TMath::TwoPi()) phiBoerYell -= TMath::TwoPi();
 
	*/

        // calculate eec quantities -------------------------------------------

        // get EEC weights
        std::pair<double, double> cst_weights = std::make_pair(
          GetCstWeight(vecCst4.first, vecJet4),
          GetCstWeight(vecCst4.second, vecJet4)
        );

        // and then calculate RL (dist b/n cst.s for EEC) and overall EEC weight
        const double dist    = Tools::GetCstDist(csts);
        const double weight  = cst_weights.first * cst_weights.second * evt_weight;

        // fill histograms ---------------------------------------------------=

        // fill histograms if needed
        if (m_manager.GetDoEECHists()) {

          // grab hist indices
          std::vector<Type::HistIndex> indices = GetHistIndices(jet);

	  // Collins/Boer-Mulders

	  /*
          // collect quantities to be histogrammed
          Type::HistContent content(weight, dist);
          if (m_manager.GetDoSpinBins()) {
            content.phiCollB = phiCollBlue;
            content.phiCollY = phiCollYell;
            content.phiBoerB = phiBoerBlue;
            content.phiBoerY = phiBoerYell;
            content.spinB    = vecSpin3.first.Y();
            content.spinY    = vecSpin3.second.Y();
            content.pattern  = jet.pattern;
          }
	  */

	  // Dihadron FF

          // collect quantities to be histogrammed
          Type::HistContent content(weight, dist);
          if (m_manager.GetDoSpinBins()) {
            content.phiCollB = ThetaSB_RC;
            content.phiCollY = ThetaSA_RC;
            content.phiBoerB = 0.0;
            content.phiBoerY = 0.0;
            content.spinB    = vecSpin3.first.Y();
            content.spinY    = vecSpin3.second.Y();
            content.pattern  = jet.pattern;
          }

          // fill spin-integrated histograms
          for (std::size_t idx = 0; idx < Const::NBinsPerSpin(); ++idx) {
            m_manager.FillEECHists(indices[idx], content);
          }

          // if needed, fill spin sorted histograms
          if (m_manager.GetDoSpinBins() && (indices.size() > Const::BlueSpinStart())) {

            // fill blue spins
            for (
              std::size_t idx = Const::BlueSpinStart();
              idx < Const::YellSpinStart();
              ++idx
            ) {
              m_manager.FillEECHists(indices[idx], content);
            }

            // fill yellow and both spins
            if (indices.size() > Const::YellSpinStart()) {
              for (
                std::size_t idx = Const::YellSpinStart();
                idx < indices.size();
                ++idx
              ) {
                m_manager.FillEECHists(indices[idx], content);
              }
            }
          }  // end spin hist filling
        }  // end hist filing
        return;

      }  // end 'CalcEEC(Type::Jet&, std::pair<Type::Cst, Type::Cst>&, double)'

      // ----------------------------------------------------------------------
      //! End calculations
      // ----------------------------------------------------------------------
      void End(TFile* file) {

        // save histograms to file
        m_manager.SaveHists(file);
        return;

      }  // end 'End(TFile*)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Calculator()  {

        m_weight_power = 1.0;
        m_weight_type  = Type::Pt;

      }  // end default ctor

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Calculator() {};

      // ----------------------------------------------------------------------
      //! ctor accepting arguments
      // ----------------------------------------------------------------------
      Calculator(const Type::Weight weight, const double power = 1.0) {

        m_weight_power = power;
        m_weight_type  = weight;

      }  // end ctor(Type::Weight, double)

  };  // end PHEnergyCorrelator::Calculator

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
