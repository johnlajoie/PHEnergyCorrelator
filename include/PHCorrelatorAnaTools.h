/// ============================================================================
/*! \file    PHCorrelatorAnaTools.h
 *  \authors Derek Anderson
 *  \date    09.23.2024
 *
 *  Useful tools for analysis/histogramming.
 */
/// ============================================================================

#ifndef PHCORRELATORANATOOLS_H
#define PHCORRELATORANATOOLS_H

// c++ utilities
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
// root libraries
#include <TH1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
// analysis components
#include "PHCorrelatorAnaTypes.h"
#include "PHCorrelatorConstants.h"



namespace PHEnergyCorrelator {
  namespace Tools {

    // ------------------------------------------------------------------------
    //! Wrapper function for std::pow
    // ------------------------------------------------------------------------
    double Exponentiate(const double arg) {

      return std::pow(Const::Base(), arg);

    }  // end 'Exponentiate(double&)'



    // ------------------------------------------------------------------------
    //! Wrapper function for std::log
    // ------------------------------------------------------------------------
    double Log(const double arg) {

      return std::log10(arg) / std::log10(Const::Base());

    }  // end 'Log(double&)'



    // ------------------------------------------------------------------------
    //! Get distance between 2 constituents
    // ------------------------------------------------------------------------
    double GetCstDist(const std::pair<Type::Cst, Type::Cst>& csts) {

      const double dist = hypot(
        csts.first.eta - csts.second.eta,
        remainder(csts.first.phi - csts.second.phi, TMath::TwoPi())
      );
      return dist;

    }  // end 'GetCstDist(std::pair<Type::Cst, Type::Cst>&)'



    // ------------------------------------------------------------------------
    //! Get variance from a standard error + counts
    // ------------------------------------------------------------------------
    double GetVariance(const double err, const double counts) {

      const double sqvar = err * sqrt(counts);
      const double var   = sqvar * sqvar;
      return var;

    }  // end 'GetVariance(double, double)'



    // ------------------------------------------------------------------------
    //! Get hadron angle controlled for wraparounds
    // ------------------------------------------------------------------------
    /*! Takes an angle in (0, pi) and constrains it to (-pi/2, pi/2).
     */
    double GetWrappedHadronAngle(double angle) {

      if (angle > TMath::PiOver2()) {
        angle -= TMath::Pi();
      }
      return angle;

    }  // end 'GetWrappedHadronAngle(double)'



    // ------------------------------------------------------------------------
    //! Get doubled hadron angle controlled for wraparounds
    // ------------------------------------------------------------------------
    /*! Takes an angle and constrains it to (-pi/2, pi/2).
     */
    double GetWrappedDoubledHadronAngle(double angle) {

      // angles for convenience
      const double piDiv2     = TMath::PiOver2();
      const double pi3Div2    = 3.0 * TMath::PiOver2();
      const double minPiDiv2  = -1.0 * TMath::PiOver2();
      const double minPi3Div2 = -3.0 * TMath::PiOver2();

      // check what adjustment is needed
      const bool doPiSub    = ( (angle > piDiv2) && (angle <= pi3Div2) );
      const bool doTwoPiSub = ( angle > pi3Div2 );
      const bool doPiAdd    = ( (angle < minPiDiv2) && (angle >= minPi3Div2) );
      const bool doTwoPiAdd = ( angle < minPi3Div2 );

      // apply adjustment
      if (doPiSub) {
        angle -= TMath::Pi();
      } else if (doTwoPiSub) {
        angle -= TMath::TwoPi();
      } else if (doPiAdd) {
        angle += TMath::Pi();
      } else if (doTwoPiAdd) {
        angle += TMath::TwoPi();
      }
      return angle;

    }  // end 'GetWrappedDoubledHadronAngle(double)'



    // ------------------------------------------------------------------------
    //! Get spin - hadron angle controlled for wraparounds
    // ------------------------------------------------------------------------
    /*! Takes an angle and constrains it to (0, pi).
     */
    double GetWrappedSpinHadronAngle(double angle) {

      // check what adjustment is needed
      const bool doPiSub = ( angle > TMath::Pi() );
      const bool doPiAdd = ( angle < 0.0 );

      // apply adjustment
      if (doPiSub) {
        angle -= TMath::Pi();
      } else if (doPiAdd) {
        angle += TMath::Pi();
      }
      return angle;

    }  // end 'GetWrappedSpinHadronAngle(double)'



    // ------------------------------------------------------------------------
    //! Divide a range into a certain number of bins
    // ------------------------------------------------------------------------
    std::vector<double> GetBinEdges(
      const std::size_t num,
      const double start,
      const double stop,
      const Type::Axis axis = Type::Log
    ) {

      // throw error if start/stop are out of order
      // or if num is zero
      if (num <= 0)     assert(num > 0);
      if (start > stop) assert(start <= stop);

      // set start/stop, calculate bin steps
      const double start_use = (axis == Type::Log) ? Log(start) : start;
      const double stop_use  = (axis == Type::Log) ? Log(stop)  : stop;
      const double step      = (stop_use - start_use) / num;
 
      // instantiate vector to hold bins
      std::vector<double> bins;

      // and fill vector
      double edge = start_use;
      for (std::size_t inum = 0; inum < num; ++inum) {
        bins.push_back( edge );
        edge += step;
      }
      bins.push_back( edge );

      // if need be, transform back from log
      if (axis == Type::Log) {
        std::transform(
          bins.begin(),
          bins.end(),
          bins.begin(),
          Exponentiate
        ); 
      }
      return bins;

    }  // end 'GetBinEdges(uint32_t, double, double)'



    // ------------------------------------------------------------------------
    //! Get jet 4-vector in cartesian coordinates from jet info
    // ------------------------------------------------------------------------
    /*! Returns jet lorentz vector based on a Type::Jet. Note that the
     *  energy of the vector is set to be the magnitude of the jet 3-
     *  vector.
     *
     *  If `norm` is set to true, vector will be normalized by
     *  its magnitude
     */
    TLorentzVector GetJetLorentz(const Type::Jet& jet, const bool norm = false) {

      // calculate momentum components
      const double th = 2.0 * atan( exp(-1.0 * jet.eta) );
      const double pz = jet.pt / tan(th);
      const double px = (pz / cos(th)) * cos(jet.phi);
      const double py = (pz / cos(th)) * sin(jet.phi);

      // normalize if need be
      TVector3 vector(px, py, pz);
      if (norm) {
        vector *= (1.0 / vector.Mag());
      }

      // return 4-vector
      return TLorentzVector(vector.Px(), vector.Py(), vector.Pz(), vector.Mag());

    }  // end 'GetJetLorentz(Type::Jet&, bool)'



    // ------------------------------------------------------------------------
    //! Get constituent 4-vector in cartesian coordinates from cst info
    // ------------------------------------------------------------------------
    /*! Returns cst lorentz vector based on a Type::Cst and the pt of a jet.
     *  Note that the energy of the vector is set to be the magnitude of the
     *  jet 3-vector.
     *
     *  If `norm` is set to true, vector will be normalized by
     *  its magnitude
     */
    TLorentzVector GetCstLorentz(
      const Type::Cst& cst,
      const double ptJet,
      const bool norm = false
    ) {

      // get total momentum
      const double ptCst = cst.z * ptJet;
      const double pCst  = sqrt( (ptCst * ptCst) + (cst.jt * cst.jt) );

      // calculate momentum components
      const double th = 2.0 * atan( exp(-1.0 * cst.eta) );
      const double px = pCst * sin(th) * cos(cst.phi);
      const double py = pCst * sin(th) * sin(cst.phi);
      const double pz = pCst * cos(th);

      // normalize if need be
      TVector3 vector(px, py, pz);
      if (norm) {
        vector *= (1.0 / vector.Mag());
      }

      // return 3-vector
      return TLorentzVector(vector.Px(), vector.Py(), vector.Pz(), vector.Mag());

    }  // end 'GetCstLorentz(Types::Cst&, double, bool)'



    // ------------------------------------------------------------------------
    //! Get magnitude-weighted average of two 3-vectors
    // ------------------------------------------------------------------------
    TVector3 GetWeightedAvgVector(
      const TVector3& va,
      const TVector3& vb,
      const bool norm = false
    ) {

      // calculate weights
      const double wa = va.Mag() / (va.Mag() + vb.Mag());
      const double wb = vb.Mag() / (va.Mag() + vb.Mag());

      // scale vectors
      const TVector3 sva = va * wa;
      const TVector3 svb = vb * wb;

      // sum and norm (if needed)
      TVector3 sum = sva + svb;
      if (norm) {
        sum *= (1.0 / sum.Mag());
      }
      return sum;

    }  // end 'GetWeightedAvgVector(TVector3& x 2, bool)'



    // ------------------------------------------------------------------------
    //! Get beam directions
    // ------------------------------------------------------------------------
    std::pair<TVector3, TVector3> GetBeams() {

      return std::make_pair(Const::BlueBeam(), Const::YellowBeam());

    }  // end 'GetBeams()'



    // ------------------------------------------------------------------------
    //! Get spins based on a provided spin pattern
    // ------------------------------------------------------------------------
    /*! Returns a pair of spin vectors based on a provided spin pattern.
     *  The 1st element will always be the blue spin, and the 2nd the
     *  yellow.
     */
    std::pair<TVector3, TVector3> GetSpins(const int pattern) { 

      TVector3 blue(0.0, 0.0, 0.0);
      TVector3 yellow(0.0, 0.0, 0.0);
      switch (pattern) {

        // blue up, yellow up (pp)
        case Type::PPBUYU:
          blue   = Const::SpinUp();
          yellow = Const::SpinUp();
          break;

        // blue down, yellow up (pp)
        case Type::PPBDYU:
          blue   = Const::SpinDown();
          yellow = Const::SpinUp();
          break;

        // blue up, yellow down (pp)
        case Type::PPBUYD:
          blue   = Const::SpinUp();
          yellow = Const::SpinDown();
          break;

        // blue down, yellow down (pp)
        case Type::PPBDYD:
          blue   = Const::SpinDown();
          yellow = Const::SpinDown();
          break;

        // blue up (pAu)
        case Type::PABU:
          blue   = Const::SpinUp();
          yellow = Const::SpinNull();
          break;

        // blue down (pAu)
        case Type::PABD:
          blue   = Const::SpinDown();
          yellow = Const::SpinNull();
          break;

        // by default, return both as null vectors
        default:
          blue   = Const::SpinNull();
          yellow = Const::SpinNull();
          break;

      }
      return std::make_pair(blue, yellow);

    }  // end 'GetSpins(int)'

  }  // end Tools namespace
}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
