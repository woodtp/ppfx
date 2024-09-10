
#include "InteractionData.h"
#include "PDGParticleCodes.h"
#include <iostream>
#include <iomanip>

namespace NeutrinoFluxReweight{

  InteractionData::InteractionData(){

    particle = TDatabasePDG::Instance();

    gen     = 0;
    Inc_pdg = 0;
    Prod_pdg= 0;

    Inc_P     = -1000.;
    Prod_P    = -1000.;
    Inc_Mass  = -1000.;
    Prod_Mass = -1000.;

    for(int i=0; i<4; i++){
      Inc_P4[i]=0;
      Prod_P4[i]=0;
      if(i<3) Vtx[i]=0;
    }

    xF      = -1000.;
    Pz      = -1000.;
    Theta   = -1000.;
    Pt      = -1000.;

    Ecm     = -1000.;
    Betacm  = -1000.;
    Gammacm = -1000.;

    Vol = "NotDefined";
    nucleus = -1;

  }

  InteractionData::InteractionData(const int genid, const double (&incMom)[3], const int incPdg, const double (&prodMom)[3], const int prodPdg, const std::string &volname, const int nucleus_pdg, const std::string &procname, const double (&vtx)[3])
  {
    particle = TDatabasePDG::Instance();
    // Z direction along the direction of the incident particle
    // Cos between incMom and prodMom:
    // The units are in GeV

    gen      = genid;

    Inc_pdg  = incPdg;
    Prod_pdg = prodPdg;

    Inc_P  = std::sqrt(incMom[0]*incMom[0] + incMom[1]*incMom[1] + incMom[2]*incMom[2]);
    Prod_P = std::sqrt(prodMom[0]*prodMom[0] + prodMom[1]*prodMom[1] + prodMom[2]*prodMom[2]);

    const double cos_theta = (incMom[0]*prodMom[0]+incMom[1]*prodMom[1]+incMom[2]*prodMom[2])/(Inc_P*Prod_P);
    const double sin_theta = std::sqrt(1.-pow(cos_theta,2.0));

    //Theta in rads:
    Theta = std::acos(cos_theta);

    Pt = Prod_P*sin_theta;
    Pz = Prod_P*cos_theta;

    constexpr int pdg_deut = 1000010020;
    if(Inc_pdg < pdg_deut) Inc_Mass = particle->GetParticle(Inc_pdg)->Mass();
    else if(Inc_pdg == pdg_deut) { Inc_Mass = 1.875; }
    else {Inc_Mass = 2.809;}

    if(Prod_pdg < pdg_deut) Prod_Mass = particle->GetParticle(Prod_pdg)->Mass();
    else if(Prod_pdg == pdg_deut) { Prod_Mass = 1.875; }
    else { Prod_Mass = 2.809; }

    static const double NUCLEON_MASS = (particle->GetParticle(pdg::P)->Mass() + particle->GetParticle(pdg::N)->Mass())/2.;
    static const double NUCLEON_MASS2 = NUCLEON_MASS*NUCLEON_MASS;
    //Ecm, gamma:
    /**
     * Center of mass of the system of projectile
     * and nucleon (not the nucleus!)
     */

    const double inc_E_lab = std::sqrt(Inc_P*Inc_P + pow(Inc_Mass,2));
    Ecm       = std::sqrt(pow(Inc_Mass,2) + NUCLEON_MASS2 + 2.*inc_E_lab*NUCLEON_MASS);
    Betacm    = std::sqrt(pow(inc_E_lab,2)-pow(Inc_Mass,2.0))/(inc_E_lab + NUCLEON_MASS);
    Gammacm   = 1./std::sqrt(1.-pow(Betacm,2.0));

    //xF:
    const double prod_E_lab  = std::sqrt(Prod_P*Prod_P + pow(Prod_Mass,2));
    const double PL          = Gammacm*(Pz-Betacm*prod_E_lab);  // PL is measured in CM frame
    xF  = PL*2./Ecm;

    //4 momenta:
    Inc_P4[3] = inc_E_lab;
    Prod_P4[3] = prod_E_lab;
    for(int i=0; i<3; i++) {Inc_P4[i]=incMom[i]; Prod_P4[i]=prodMom[i];}


    //Volume:
    Vol = volname;

    //target nucleus z
    if(!(nucleus_pdg / 1000000000)) {
//       std::cout<<"Error! "<<nucleus_pdg<<" is not a nucleus code! Vol = "<<volname<<", Proc = "<<procname<<std::endl;
      nucleus = -1;
    }
    else {
      nucleus = (nucleus_pdg/10000)%1000;
    }

    //Process:
    Proc = procname;

    //Vertex:
    for(int i=0; i<3; i++) Vtx[i] = vtx[i];

  }

  //----------------------------------------------------------------------
  InteractionData::~InteractionData(){

  }

  std::ostream& operator<<(std::ostream& os, const InteractionData& id) {
    return id.print(os);
  }

  std::ostream& InteractionData::print(std::ostream& os) const {
    using namespace std;
    os<<"IN:"<<setw(5)<<Inc_pdg
      <<" | p3: {";
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(6)<<Inc_P4[i]<<" ";
    }
    os<<"}\nOUT:"<<setw(5)<<Prod_pdg
      <<" | p3: {"<<setiosflags(ios::fixed) << setprecision(2);
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(6)<<Prod_P4[i]<<" ";
    }
    os <<"} | v3: {";
    for(int i=0; i<3; i++) {
      os<<setiosflags(ios::fixed) << setprecision(2)<<setw(5)<<Vtx[i]<<" ";
    }
    os<<"} | xF, pT: "<<xF<<", "<<Pt;
    os<<endl;
    return os;
  }

}
