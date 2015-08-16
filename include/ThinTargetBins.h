#ifndef THINTARGETBINS_H
#define THINTARGETBINS_H

#include <utility>
#include <map>
#include <vector>

namespace NeutrinoFluxReweight{
  
  /*! \class ThinTargetBins
   *  \brief A class to manage the bin definitions for MIPP Numi Yields
   */
  class ThinTargetBins{
  public:
    
    ThinTargetBins();
    static ThinTargetBins* getInstance();
    
    //! Read a NA49 data pip xml file name to parse the bins
    void pC_pi_from_xml(const char* filename);

    //!Barton:
    void barton_pC_pi_from_xml(const char* filename);

    //! Read a NA49 data K xml file name to parse the bins
    void pC_k_from_xml(const char* filename);
    
    //!MIPP k/pi:
    void mipp_pC_k_pi_from_xml(const char* filename);
    
    //! Read a pion incident 
    void meson_incident_from_xml(const char* filename);
    
    //! Return the Bin ID for this data
    int BinID_pC_pi(double xf, double pt,int pdgcode);
    
    //! Return the Bin ID for this data
    int barton_BinID_pC_pi(double xf, double pt,int pdgcode);
    
    //! Return the Bin ID for this data
    int BinID_pC_k(double xf, double pt,int pdgcode);
    
    //! Return the MIPP Thin Target Bin ID for this data
    int mipp_BinID_pC_k(double pz, double pt,int pdgcode);
    
    //! Return Pion incident bin
    int meson_inc_BinID(double xf, double pt,int pdgcode);
    
    //Pions production:
    std::vector<double> pC_pi_xfmin, pC_pi_xfmax, pC_pi_ptmin, pC_pi_ptmax;
    std::vector<double> b_pC_pi_xfmin, b_pC_pi_xfmax, b_pC_pi_ptmin, b_pC_pi_ptmax;
    
    //Kaons production:
    std::vector<double> pC_k_xfmin, pC_k_xfmax, pC_k_ptmin, pC_k_ptmax;
    std::vector<double> mipp_pC_k_pzmin, mipp_pC_k_pzmax, mipp_pC_k_ptmin, mipp_pC_k_ptmax;
    
    //Pion incident
    std::vector<double> meson_inc_xfmin, meson_inc_xfmax, meson_inc_ptmin, meson_inc_ptmax;
    
  private:
    static ThinTargetBins* instance;
    
  };

  
  
};
#endif
