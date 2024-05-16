#include "ReweightDriver.h"
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace NeutrinoFluxReweight
{

ReweightDriver::ReweightDriver(int iuniv, const ParameterTable &cv_pars, const ParameterTable &univ_pars, std::string fileIn) :
    iUniv(iuniv), cvPars(cv_pars), univPars(univ_pars), fileOptions(fileIn)
{
    ParseOptions();
    Configure();
}

void ReweightDriver::Configure()
{

    //Creating the vector of reweighters:

    if (doMIPPNumi)
    {
        MIPP_NUMI_PION_Universe = new MIPPNumiPionYieldsReweighter(iUniv, cvPars, univPars);
        MIPP_NUMI_KAON_Universe = new MIPPNumiKaonYieldsReweighter(iUniv, cvPars, univPars);
    }

    TARG_ATT_Universe = new TargetAttenuationReweighter(iUniv, cvPars, univPars);
    VOL_ABS_IC_Universe = new AbsorptionICReweighter(iUniv, cvPars, univPars);
    VOL_ABS_DPIP_Universe = new AbsorptionDPIPReweighter(iUniv, cvPars, univPars);
    VOL_ABS_DVOL_Universe = new AbsorptionDVOLReweighter(iUniv, cvPars, univPars);
    VOL_ABS_NUCLEON_Universe = new NucleonAbsorptionOutOfTargetReweighter(iUniv, cvPars, univPars);
    VOL_ABS_OTHER_Universe = new OtherAbsorptionOutOfTargetReweighter(iUniv, cvPars, univPars);

    THINTARGET_PC_PION_Universe = new ThinTargetpCPionReweighter(iUniv, cvPars, univPars);
    THINTARGET_PC_KAON_Universe = new ThinTargetpCKaonReweighter(iUniv, cvPars, univPars);
    THINTARGET_NC_PION_Universe = new ThinTargetnCPionReweighter(iUniv, cvPars, univPars);
    THINTARGET_PC_NUCLEON_Universe = new ThinTargetpCNucleonReweighter(iUniv, cvPars, univPars);
    THINTARGET_MESON_INCIDENT_Universe = new ThinTargetMesonIncidentReweighter(iUniv, cvPars, univPars);
    THINTARGET_NUCLEON_A_Universe = new ThinTargetnucleonAReweighter(iUniv, cvPars, univPars);
    OTHER_Universe = new OtherReweighter(iUniv, cvPars, univPars);
}

void ReweightDriver::ParseOptions()
{
    //Parsing the file input:
    using boost::property_tree::ptree;
    ptree top;
    std::string val = "";
    read_xml(fileOptions.c_str(), top, 2); // option 2 removes comment strings
    ptree &options = top.get_child("inputs.Settings");

    doMIPPNumi = "MIPPNuMIOn" == options.get<std::string>("Reweighters");
}

double ReweightDriver::calculateWeight(const InteractionChainData &icd)
{

    //Boolean flags:
    std::vector<bool> interaction_nodes(icd.interaction_chain.size(), false);

    bool has_mipp = false;
    auto reweight_MIPP = [&interaction_nodes, &icd, &has_mipp](IInteractionChainReweighting *const reweighter) -> double
    {
        interaction_nodes = reweighter->canReweight(icd);
        for (auto const &node : interaction_nodes)
        {
            if (!node)
                continue;
            has_mipp = true;
            return reweighter->calculateWeight(icd);
        }
        return 1.0;
    };

    auto const &interaction_chain = icd.interaction_chain;
    auto reweight = [&interaction_nodes, &interaction_chain](IInteractionReweighting *const reweighter) -> double
    {
        double wgt = 1.0;
        for (auto it = interaction_nodes.rbegin(); it != interaction_nodes.rend(); ++it)
        {
            if (*it)
                continue;
            auto const &aa = interaction_chain[std::distance(interaction_nodes.rbegin(), it)];
            if (!reweighter->canReweight(aa))
                continue;
            wgt *= reweighter->calculateWeight(aa);
            *it = true;
        }
        return wgt;
    };

    auto reweight_att_abs = [&icd](IInteractionChainReweighting *const reweighter) -> double
    {
        auto const &nodes = reweighter->canReweight(icd);
        //we just see for the first position (primary proton)
        if (nodes.empty() || !nodes[0])
            return 1.0;
        return reweighter->calculateWeight(icd);
    };

    /// ----- PROCESS INTERACTION NODES ----- ///

    //MIPP NuMI Pions:
    mipp_pion_wgt = doMIPPNumi ? reweight_MIPP(MIPP_NUMI_PION_Universe) : 1.0;

    //MIPP NuMI Kaons:
    mipp_kaon_wgt = doMIPPNumi && !has_mipp ? reweight_MIPP(MIPP_NUMI_KAON_Universe) : 1.0;

    //Thin Target pC->piX:
    pC_pi_wgt = reweight(THINTARGET_PC_PION_Universe);

    //Thin Target pC->KX:
    pC_k_wgt = reweight(THINTARGET_PC_KAON_Universe);

    //Thin Target nC->piX:
    nC_pi_wgt = reweight(THINTARGET_NC_PION_Universe);

    //Thin Target pC->nucleonX:
    pC_nu_wgt = reweight(THINTARGET_PC_NUCLEON_Universe);

    //Thin Target Meson Incident:
    meson_inc_wgt = reweight(THINTARGET_MESON_INCIDENT_Universe);

    //Thin Target Nucleon Incident not hanldle NA49 or Barton:
    // nuA_wgt = reweight(THINTARGET_NC_PION_Universe);
    nuA_wgt = 1.0;
    nuA_dvol_wgt = 1.0;
    nuA_othervol_wgt = 1.0;
    for (auto it = interaction_nodes.rbegin(); it != interaction_nodes.rend(); ++it)
    {
        if (*it)
            continue;
        auto const &aa = icd.interaction_chain[std::distance(interaction_nodes.rbegin(), it)];
        if (!THINTARGET_NUCLEON_A_Universe->canReweight(aa))
            continue;

        double rewval = THINTARGET_NUCLEON_A_Universe->calculateWeight(aa);

        nuA_wgt *= rewval;

        const bool is_data_vol = (aa.Vol == "TGT1") || (aa.Vol == "BudalMonitor") || (aa.Vol == "Budal_HFVS") || (aa.Vol == "Budal_VFHS");

        if (is_data_vol)
            nuA_dvol_wgt *= rewval; // Interactions occuring primarily with carbon
        else
            nuA_othervol_wgt *= rewval; // Interactions everywhere else

        *it = true;
    }

    //Any other interaction not handled yet:
    other_wgt = reweight(OTHER_Universe);

    //Target attenuation correction:
    att_wgt = reweight_att_abs(TARG_ATT_Universe);

    //Absorption correction:
    tot_abs_wgt = 1.0;

    // Correction of the pi & K absorption in volumes (Al)
    abs_ic_wgt = reweight_att_abs(VOL_ABS_IC_Universe);

    //Correction of the pi & K absorption in volumes (Fe)
    abs_dpip_wgt = reweight_att_abs(VOL_ABS_DPIP_Universe);

    //Correction of the pi & K absorption in volumes (He)
    abs_dvol_wgt = reweight_att_abs(VOL_ABS_DVOL_Universe);

    //Correction of nucleons on Al, Fe and He.
    abs_nucleon_wgt = reweight_att_abs(VOL_ABS_NUCLEON_Universe);

    //Correction of any other particle on Al, Fe and He.
    abs_other_wgt = reweight_att_abs(VOL_ABS_OTHER_Universe);

    tot_abs_wgt *= abs_ic_wgt * abs_dpip_wgt * abs_dvol_wgt * abs_nucleon_wgt * abs_other_wgt;

    const double tot_wgt = mipp_pion_wgt * mipp_kaon_wgt * pC_pi_wgt * pC_k_wgt * pC_nu_wgt * meson_inc_wgt * nuA_wgt * other_wgt * att_wgt * tot_abs_wgt;

    if (tot_wgt != tot_wgt)
    {
        std::cout << "Alert nan total wgt... check!!!" << std::endl;
        return 1.0;
    }
    return tot_wgt;
}

ReweightDriver::~ReweightDriver()
{

    if (doMIPPNumi)
    {
        delete MIPP_NUMI_PION_Universe;
        delete MIPP_NUMI_KAON_Universe;
    }
    delete TARG_ATT_Universe;
    delete VOL_ABS_IC_Universe;
    delete VOL_ABS_DPIP_Universe;
    delete VOL_ABS_DVOL_Universe;
    delete VOL_ABS_NUCLEON_Universe;
    delete VOL_ABS_OTHER_Universe;
    delete THINTARGET_PC_PION_Universe;
    delete THINTARGET_PC_KAON_Universe;
    delete THINTARGET_NC_PION_Universe;
    delete THINTARGET_PC_NUCLEON_Universe;
    delete THINTARGET_MESON_INCIDENT_Universe;
    delete THINTARGET_NUCLEON_A_Universe;
    delete OTHER_Universe;
}

} // namespace NeutrinoFluxReweight
