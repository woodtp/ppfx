
#include "ThinTargetpCPionReweighter.h"
#include "CentralValuesAndUncertainties.h"
#include "ThinTargetMC.h"
#include "ThinTargetBins.h"
#include "DataHistos.h"

#include <iostream>
namespace NeutrinoFluxReweight{
  
  ThinTargetpCPionReweighter::ThinTargetpCPionReweighter(int iuniv, const ParameterTable& cv_pars, const ParameterTable& univ_pars):iUniv(iuniv),cvPars(cv_pars),univPars(univ_pars){
    
    std::map<std::string, double> cv_table   = cvPars.table;
    std::map<std::string, double> univ_table = univPars.table;
    
    double univ_data_prod_xs = univ_table["prod_prtC_xsec"];
    double cv_data_prod_xs = cv_table["prod_prtC_xsec"];
    
    //the number of bins needs to be written from the xmls files 
    char namepar[100];
    for(int ii=0;ii<9801;ii++){
      
      sprintf(namepar,"ThinTarget_pC_%s_sys_%d","pip",ii);
      double data_cv  = cv_table[std::string(namepar)];
      double data_sys = univ_table[std::string(namepar)];
      sprintf(namepar,"ThinTarget_pC_%s_stats_%d","pip",ii);
      double data_sta = univ_table[std::string(namepar)];
      
      data_cv  /= cv_data_prod_xs;
      data_sys /= univ_data_prod_xs;
      data_sta /= univ_data_prod_xs;      
      vbin_data_pip.push_back(data_sta + data_sys - data_cv);
      
      sprintf(namepar,"ThinTarget_pC_%s_sys_%d","pim",ii);
      data_cv  = cv_table[std::string(namepar)];
      data_sys = univ_table[std::string(namepar)];
      sprintf(namepar,"ThinTarget_pC_%s_stats_%d","pim",ii);
      data_sta = univ_table[std::string(namepar)];
      
      data_cv  /= cv_data_prod_xs;
      data_sys /= univ_data_prod_xs;
      data_sta /= univ_data_prod_xs;      
      vbin_data_pim.push_back(data_sta + data_sys - data_cv);
      
    }    
    
  }
  
   ThinTargetpCPionReweighter::~ThinTargetpCPionReweighter(){
    
  }
  bool ThinTargetpCPionReweighter::canReweight(const InteractionData& aa){

    //checking:
    if(aa.Inc_pdg != 2212)return false;
    if(aa.Inc_P < 12.0)return false;
    if(aa.Vol != "TGT1" && aa.Vol != "BudalMonitor")return false;
    
    ThinTargetBins*  Thinbins =  ThinTargetBins::getInstance();
    int bin = Thinbins->BinID_pC_X(aa.xF,aa.Pt,aa.Prod_pdg);
    if(bin<0)return false;
    if(bin>=0)return true;
  }
  
  double ThinTargetpCPionReweighter::calculateWeight(const InteractionData& aa){
    
    double wgt = 1.0;
    
    ThinTargetBins*  Thinbins =  ThinTargetBins::getInstance();
    int bin = Thinbins->BinID_pC_X(aa.xF,aa.Pt,aa.Prod_pdg);
    
    if(bin<0)return wgt;;
      
    //Calculating the scale:
    double data_scale = calculateDataScale(aa.Inc_pdg,aa.Inc_P,aa.Prod_pdg,aa.xF,aa.Pt);
    
    if(aa.Prod_pdg == 211)vbin_data_pip[bin]  *= data_scale;
    else if(aa.Prod_pdg ==-211)vbin_data_pim[bin]  *= data_scale;
    else{
      std::cout<<"Wrong input, pdg_prod: "<< aa.Prod_pdg  <<std::endl;
      return wgt;
    }
    
    ThinTargetMC*  mc =  ThinTargetMC::getInstance();
    double mc_cv = mc->getMCval_pC_X(aa.Inc_P,aa.xF,aa.Pt,aa.Prod_pdg);
    
    mc_cv /= calculateMCProd(aa.Inc_P);
    
    if(mc_cv<1.e-12)return wgt;
    
    if(aa.Prod_pdg == 211)wgt = vbin_data_pip[bin]/mc_cv;
    if(aa.Prod_pdg == 211)wgt = vbin_data_pip[bin]/mc_cv;
    
    return wgt;
  }
  
  double ThinTargetpCPionReweighter::calculateMCProd(double inc_mom){
    double xx[13] ={12,20,31,40,50,60,70,80,90,100,110,120,158};
    double yy[13] ={153386793./197812683.,160302538./197811564.,164508480./197831250.,166391359./197784915.,
		    167860919./197822312.,168882647./197807739.,169681805./197803099.,170311264./197811098.,
		    170860912./197822002.,171309291./197834756.,171651963./197811822.,171991260./197823012.,
		  172902228./197804669.};
    
    int idx_lowp = -1;
    int idx_hip  = -1;
    for(int i=0;i<12;i++){
      if(inc_mom>=double(xx[i]) && inc_mom<double(xx[i+1])){
	idx_lowp=i;
	idx_hip =i+1;
      }
    }
    double frac_low = yy[idx_lowp];
    double frac_hi  = yy[idx_hip];
    double frac_m   =  frac_low + (inc_mom-double(xx[idx_lowp]))*(frac_hi-frac_low)/(double(xx[idx_hip])-double(xx[idx_lowp]));
    
    return frac_m*243.2435;
    
  }
  double ThinTargetpCPionReweighter::calculateDataScale(int inc_pdg, double inc_mom, int prod_pdg,double xf, double pt){
    double scaling_violation = 1.0;
    DataHistos*  dtH =  DataHistos::getInstance();
    //temporary:
    const int Nscl = 11;
    const int moms[Nscl] = {12,20,31,40,50,60,70,80,100,120,158};
    
    int idx_part = -1;
    if(prod_pdg == 211)idx_part = 0;
    if(prod_pdg ==-211)idx_part = 1;
    if(idx_part<0){
      std::cout<<"Error in the prod particle"<<std::endl;
      return 1.0;
    }
    
    int binid = dtH->hNA49Scl[idx_part][Nscl-1]->FindBin(xf,pt);
    double scl_ref158 = dtH->hNA49Scl[idx_part][Nscl-1]->GetBinContent(binid);    
    
    int idx_lowp = -1;
    int idx_hip  = -1;
    for(int i=0;i<Nscl-1;i++){
      if(inc_mom>=double(moms[i]) && inc_mom<double(moms[i+1])){
	idx_lowp=i;
	idx_hip =i+1;
      }
    }
    if(idx_lowp<0 || idx_hip<0){
      std::cout<<"Error calculating the scaling"<<std::endl;
      return 1.0;
    }
    double scl_low = dtH->hNA49Scl[idx_part][idx_lowp]->GetBinContent(binid);
    double scl_hi  = dtH->hNA49Scl[idx_part][idx_hip]->GetBinContent(binid);
    double scl_m   =  scl_low + (inc_mom-double(moms[idx_lowp]))*(scl_hi-scl_low)/(double(moms[idx_hip])-double(moms[idx_lowp]));
    if(scl_ref158<1.e-10){
      // std::cout<<"ref158 zero!!! "<<scl_ref158<<std::endl;
      return 1.0;
    }
    scaling_violation = scl_m/scl_ref158;
    return scaling_violation;
  }

}