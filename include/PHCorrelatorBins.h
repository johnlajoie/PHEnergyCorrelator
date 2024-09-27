/// ============================================================================
/*! \file    PHCorrelatorBins.h
 *  \authors Derek Anderson
 *  \date    09.24.2024
 *
 *  Classes to handle binning of histograms filed
 *  during ENC calculations.
 */
/// ============================================================================

#ifndef PHCORRELATORBINS_H
#define PHCORRELATORBINS_H

// c++ utilities
#include <map>
#include <string>
#include <vector>
// root libraries
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
// analysis components
#include "PHCorrelatorTools.h"
#include "PHCorrelatorTypes.h"



namespace PHEnergyCorrelator {

  // ==========================================================================
  //! Binning definition
  // ==========================================================================
  /*! A small class to consolidate data
   *  for defining histogram bins.
   */ 
  class Binning {

    private:

      // data members
      double              m_start;
      double              m_stop;
      std::size_t         m_num;
      std::vector<double> m_bins;

    public:

      // ----------------------------------------------------------------------
      //! Uniform bin getters
      //-----------------------------------------------------------------------
      double      GetStart() const {return m_start;}
      double      GetStop()  const {return m_stop;}
      std::size_t GetNum()   const {return m_num;}

      // ----------------------------------------------------------------------
      //! Variable bin getter
      // ----------------------------------------------------------------------
      std::vector<double> GetBins() const {return m_bins;}

      // ----------------------------------------------------------------------
      //! default ctor/dtor
      // ----------------------------------------------------------------------
      Binning()  {};
      ~Binning() {};

      // ----------------------------------------------------------------------
      //! ctor accepting uniform parameters
      // ----------------------------------------------------------------------
      Binning(
        const std::size_t num,
        const double start,
        const double stop,
        const Type::Axis axis = Type::Norm
      ) {

        m_num   = num;
        m_start = start;
        m_stop  = stop;
        m_bins  = Tools::GetBinEdges(m_num, m_start, m_stop, axis);

      }  // end ctor(uint32_t, double, double, Type::Axis)

      // ----------------------------------------------------------------------
      //! ctor accepting non-uniform parameters
      // ----------------------------------------------------------------------
      Binning(const std::vector<double> edges) {

        m_bins  = edges;
        m_num   = edges.size() - 1;
        m_start = edges.front();
        m_stop  = edges.back();

      }  // end ctor(std::vector<double>)

  };  // end Binning



  // ==========================================================================
  //! Bin Database
  // ==========================================================================
  /*! A class to centralize binning for various quantities like R_{L}, etc.
   *  Methods are provided to update exisiting/add new bin definitions
   *  on the fly.
   */
  class Bins {

    private:

      // data members
      std::map<std::string, Binning> m_bins;

    public:

      // ----------------------------------------------------------------------
      //! Add a binning
      // ----------------------------------------------------------------------
      void Add(const std::string& name, const Binning& binning) {

        // throw error if binning already exists
        if (m_bins.count(name) >= 1) {
          assert(m_bins.count(name) == 0);
        }

        // otherwise insert new binning
        m_bins[name] = binning;
        return;

      }  // end 'Add(std::string&, Binning&)'

      // ----------------------------------------------------------------------
      //! Change a binning
      // ----------------------------------------------------------------------
      void Set(const std::string& variable, const Binning& binning) {

        // throw error if binning doesn't exist
        if (m_bins.count(variable) == 0) {
          assert(m_bins.count(variable) >= 1);
        }

        // otherwise update binning
        m_bins[variable] = binning;
        return;

      }  // end 'Set(std::string&, Binning&)'

      // ----------------------------------------------------------------------
      //! Get a binning
      // ----------------------------------------------------------------------
      Binning Get(const std::string& variable) {

        // throw error if binning doesn't exist
        if (m_bins.count(variable) == 0) {
          assert(m_bins.count(variable) == 0);
        }

        // otherwise return binning
        return m_bins[variable];

      }  // end 'Get(std::string&)'

      // ----------------------------------------------------------------------
      //! default ctor
      // ----------------------------------------------------------------------
      Bins() {

        /* TODO
         *  - Add cos(angle)
         *  - Add angle
         *  - Add xi
         */
        m_bins["energy"]  = Binning(202, -1., 100.);
        m_bins["side"]    = Binning(75, 1e-5, 1., Type::Log);
        m_bins["logside"] = Binning(75, -5., 0.);

      }  // end 'ctor()'

      // ----------------------------------------------------------------------
      //! default dtor
      // ----------------------------------------------------------------------
      ~Bins() {};

  };  // end Bins

}  // end PHEnergyCorrelator namespace

#endif

// end ========================================================================
