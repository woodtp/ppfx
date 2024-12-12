#include "ReweightDriver.h"
#include "PDGParticleCodes.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <iostream>

// #define NUA_XF_MIRRORING

namespace NeutrinoFluxReweight {

ReweightDriver::ReweightDriver(int iuniv,
                               const ParameterTable& cv_pars,
                               const ParameterTable& univ_pars,
                               std::string fileIn)
  : iUniv(iuniv)
  , cvPars(cv_pars)
  , univPars(univ_pars)
  , fileOptions(fileIn)
{
    ParseOptions();
    Configure();
}

void ReweightDriver::Configure()
{
    // Creating the vector of reweighters:

    if (doMIPPNumi) {
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
    // Parsing the file input:
    using boost::property_tree::ptree;
    ptree top;
    std::string val = "";
    read_xml(fileOptions.c_str(), top, 2); // option 2 removes comment strings
    ptree& options = top.get_child("inputs.Settings");

    doMIPPNumi = "MIPPNuMIOn" == options.get<std::string>("Reweighters");
}

double ReweightDriver::checkWeight(double wgt, const std::vector<InteractionData>& interaction_chain)
{
    auto dumpNodes = [&interaction_chain]() {
        int i = 0;
        for(const auto& aa : interaction_chain) {
            std::cout << " ======== INTERACTION NODE [" << i << "] ========\n" << aa;
            i++;
        }
        std::cout << " ======================================\n\n";
    };

    if (wgt > 20.0) {
        if (wgt > 100.0) {
            std::cout << "\n[WARNING] wgt > 100.0 (wgt = " << wgt << "). Setting weight to 0 and skipping this entry.\n";
            dumpNodes();

            mipp_pion_wgt = 0.0;
            mipp_kaon_wgt = 0.0;
            pC_pi_wgt = 0.0;
            pC_k_wgt = 0.0;
            nC_pi_wgt = 0.0;
            pC_nu_wgt = 0.0;
            meson_inc_wgt = 0.0;
            nuA_wgt = 0.0;
            other_wgt = 0.0;
            att_wgt = 0.0;
            tot_abs_wgt = 0.0;

            return 0.0;
        }
        std::cout << "\n[WARNING] wgt > 20.0 (wgt = " << wgt << ").\n";
        dumpNodes();
    }
    return wgt;
}


double ReweightDriver::calculateWeight(const InteractionChainData& icd)
{
    // Flags to keep track of nodes for which a weight has been calculated
    std::vector<bool> interaction_nodes(icd.interaction_chain.size(), false);

    /// ----- PROCESS INTERACTION NODES ----- ///

    // MIPP NuMI Pions:
    m_hasMIPP = false;
    mipp_pion_wgt = getWeightMIPP(MIPP_NUMI_PION_Universe, icd, interaction_nodes);

    // MIPP NuMI Kaons:
    mipp_kaon_wgt = getWeightMIPP(MIPP_NUMI_KAON_Universe, icd, interaction_nodes);

    // Thin Target pC->piX:
    pC_pi_wgt = getWeight(THINTARGET_PC_PION_Universe, icd.interaction_chain, interaction_nodes);

    // Thin Target pC->KX:
    pC_k_wgt = getWeight(THINTARGET_PC_KAON_Universe, icd.interaction_chain, interaction_nodes);

    // Thin Target nC->piX:
    nC_pi_wgt = getWeight(THINTARGET_NC_PION_Universe, icd.interaction_chain, interaction_nodes);

    // Thin Target pC->nucleonX:
    pC_nu_wgt = getWeight(THINTARGET_PC_NUCLEON_Universe, icd.interaction_chain, interaction_nodes);

    // Thin Target Meson Incident:
    meson_inc_wgt = getWeight(THINTARGET_MESON_INCIDENT_Universe, icd.interaction_chain, interaction_nodes);

    // Thin Target Nucleon Incident not hanldle NA49 or Barton:
    nuA_wgt = getWeight(THINTARGET_NUCLEON_A_Universe, icd.interaction_chain, interaction_nodes);

    // Any other interaction not handled yet:
    other_wgt = getWeight(OTHER_Universe, icd.interaction_chain, interaction_nodes);

    // Target attenuation correction:
    att_wgt = getWeightAttenuation(TARG_ATT_Universe, icd);

    // Absorption correction:
    // Correction of the pi & K absorption in volumes (Al)
    abs_ic_wgt = getWeightAttenuation(VOL_ABS_IC_Universe, icd);

    // Correction of the pi & K absorption in volumes (Fe)
    abs_dpip_wgt = getWeightAttenuation(VOL_ABS_DPIP_Universe, icd);

    // Correction of the pi & K absorption in volumes (He)
    abs_dvol_wgt = getWeightAttenuation(VOL_ABS_DVOL_Universe, icd);

    // Correction of nucleons on Al, Fe and He.
    abs_nucleon_wgt = getWeightAttenuation(VOL_ABS_NUCLEON_Universe, icd);

    // Correction of any other particle on Al, Fe and He.
    abs_other_wgt = getWeightAttenuation(VOL_ABS_OTHER_Universe, icd);

    tot_abs_wgt = abs_ic_wgt * abs_dpip_wgt * abs_dvol_wgt * abs_nucleon_wgt * abs_other_wgt;

    const double tot_wgt = mipp_pion_wgt * mipp_kaon_wgt * pC_pi_wgt * pC_k_wgt * pC_nu_wgt *
                           meson_inc_wgt * nuA_wgt * other_wgt * att_wgt * tot_abs_wgt;

    if (tot_wgt != tot_wgt) {
        std::cout << "Alert nan total wgt... check!!!" << std::endl;
        return 1.0;
    }
    return checkWeight(tot_wgt, icd.interaction_chain);
}

double ReweightDriver::getWeight(IInteractionReweighting* reweighter,
                                 const std::vector<InteractionData>& interaction_chain,
                                 std::vector<bool>& interaction_nodes)
{
    double wgt = 1.0;
    for (auto iNode = interaction_nodes.rbegin(); iNode != interaction_nodes.rend(); ++iNode) {
        if (*iNode) continue;
        const auto& aa = interaction_chain[interaction_nodes.rend() - 1 - iNode];
        if (!reweighter->canReweight(aa)) continue;
        wgt *= reweighter->calculateWeight(aa);
        *iNode = true;
    }
    return checkWeight(wgt, interaction_chain);
}

double ReweightDriver::getWeight(ThinTargetMesonIncidentReweighter* reweighter,
                                 const std::vector<InteractionData>& interaction_chain,
                                 std::vector<bool>& interaction_nodes)
{
    meson_inc_incoming_pip_wgt = 1.0;
    meson_inc_incoming_pim_wgt = 1.0;
    meson_inc_incoming_Kp_wgt = 1.0;
    meson_inc_incoming_Km_wgt = 1.0;
    meson_inc_incoming_K0_wgt = 1.0;
    meson_inc_outgoing_pip_wgt = 1.0;
    meson_inc_outgoing_pim_wgt = 1.0;
    meson_inc_outgoing_Kp_wgt = 1.0;
    meson_inc_outgoing_Km_wgt = 1.0;
    meson_inc_outgoing_K0_wgt = 1.0;

    double wgt = 1.0;

    for (auto iNode = interaction_nodes.rbegin(); iNode != interaction_nodes.rend(); ++iNode) {
        if (*iNode) continue;

        const auto& aa = interaction_chain[interaction_nodes.rend() - 1 - iNode];

        if (!reweighter->canReweight(aa)) continue;

        const double rewval = reweighter->calculateWeight(aa);

        switch (aa.Inc_pdg) {
            case pdg::PIP:
                meson_inc_incoming_pip_wgt *= rewval;
                break;
            case pdg::PIM:
                meson_inc_incoming_pim_wgt *= rewval;
                break;
            case pdg::KP:
                meson_inc_incoming_Kp_wgt *= rewval;
                break;
            case pdg::KM:
                meson_inc_incoming_Km_wgt *= rewval;
                break;
            case pdg::K0L:
            case pdg::K0S:
                meson_inc_incoming_K0_wgt *= rewval;
                break;
            default:
                break;
        }
        switch(aa.Prod_pdg) {
            case pdg::PIP:
                meson_inc_outgoing_pip_wgt *= rewval;
                break;
            case pdg::PIM:
                meson_inc_outgoing_pim_wgt *= rewval;
                break;
            case pdg::KP:
                meson_inc_outgoing_Kp_wgt *= rewval;
                break;
            case pdg::KM:
                meson_inc_outgoing_Km_wgt *= rewval;
                break;
            case pdg::K0L:
            case pdg::K0S:
                meson_inc_outgoing_K0_wgt *= rewval;
                break;
            default:
                break;
        }

        wgt *= rewval;
        *iNode = true;
    }

    return checkWeight(wgt, interaction_chain);
}

double ReweightDriver::getWeight(ThinTargetnucleonAReweighter* reweighter,
                                 const std::vector<InteractionData>& interaction_chain,
                                 std::vector<bool>& interaction_nodes)
{
    nuA_inC_inPS_wgt = 1.0;
    nuA_inC_OOPS_wgt = 1.0;
    nuA_outC_Ascale_wgt = 1.0;
    nuA_outC_OOPS_wgt = 1.0;
    nuA_other_wgt = 1.0;

    double wgt = 1.0;

    for (auto iNode = interaction_nodes.rbegin(); iNode != interaction_nodes.rend(); ++iNode) {
        if (*iNode) continue;

        const auto& aa = interaction_chain[interaction_nodes.rend() - 1 - iNode];

        if (!reweighter->canReweight(aa)) continue;

#ifdef UH_ICARUS
#pragma message("[INFO] UH_ICARUS is enabled. Ignoring QEL-like interactions.")
        if (aa.Inc_pdg == aa.Prod_pdg && (aa.xF > 0.95 || (aa.Pt < aa.xF - 0.5))) {
            // if QEL-like, skip
            *iNode = true;
            continue;
        }
#endif

        double rewval = reweighter->calculateWeight(aa);

        const bool is_carbon_vol = aa.Vol == "TGT1"
                                || aa.Vol == "BudalMonitor"
                                || aa.Vol == "Budal_HFVS"
                                || aa.Vol == "Budal_VFHS";

        const bool is_oops = aa.xF < 0.0
                          || aa.Pt > 2.0
                          || aa.Inc_P < 12.0;

        if (reweighter->isDataBased(aa)) {
            nuA_outC_Ascale_wgt *= rewval;
        }
        else if (is_carbon_vol) {
            if (aa.xF < 0.0 && aa.xF >= -0.25 && aa.Pt <= 2.0 && aa.Inc_P >= 12.0) {
#ifdef NUA_XF_MIRRORING
#pragma message("[WARNING] NUA_XF_MIRRORING is enabled. This was just a test and shouldn't be used.")
                // A. Wood (apwood@central.uh.edu)
                // 2024-12-12
                // This was a test to approximate what would happen if we had HP data for interactions
                // where xF < 0.0. Could remove, but leaving for now in case the UH Group wants to revisit.
                const double inc_mom[3] = { aa.Inc_P4[0], aa.Inc_P4[1], aa.Inc_P4[2] };
                const double prod_mom[3] = { aa.Prod_P4[0], aa.Prod_P4[1], aa.Prod_P4[2] };
                const double vtx_int[3] = { aa.Vtx[0], aa.Vtx[1], aa.Vtx[2] };

                InteractionData aa_tmp(
                    aa.gen,
                    inc_mom,
                    aa.Inc_pdg,
                    prod_mom,
                    aa.Prod_pdg,
                    aa.Vol,
                    aa.nucleus,
                    aa.Proc,
                    vtx_int
                );

                aa_tmp.xF = std::abs(aa_tmp.xF);

                // Re-check the modified interaction node against the previous reweighters
                // as well as the material scaling.
                if(THINTARGET_PC_PION_Universe->canReweight(aa_tmp)) {
                    rewval = THINTARGET_PC_PION_Universe->calculateWeight(aa_tmp);
                }
                else if(THINTARGET_PC_KAON_Universe->canReweight(aa_tmp)) {
                    rewval = THINTARGET_PC_KAON_Universe->calculateWeight(aa_tmp);
                }
                else if(THINTARGET_NC_PION_Universe->canReweight(aa_tmp)) {
                    rewval = THINTARGET_NC_PION_Universe->calculateWeight(aa_tmp);
                }
                else if(THINTARGET_PC_NUCLEON_Universe->canReweight(aa_tmp)) {
                    rewval = THINTARGET_PC_NUCLEON_Universe->calculateWeight(aa_tmp);
                }
                else if(THINTARGET_MESON_INCIDENT_Universe->canReweight(aa_tmp)) {
                    rewval = THINTARGET_MESON_INCIDENT_Universe->calculateWeight(aa_tmp);
                }
                else if (THINTARGET_NUCLEON_A_Universe->isDataBased(aa_tmp)) {
                    rewval = THINTARGET_NUCLEON_A_Universe->calculateWeight(aa_tmp);
                }
#endif
                nuA_inC_inPS_wgt *= rewval;
            }
            else if (is_oops) {
                nuA_inC_OOPS_wgt *= rewval;
            }
            else {
                // Contributions from exotic processes (e.g., 2212 + C -> -2112 + X)
                nuA_other_wgt *= rewval;
            }
        }
        else {
            nuA_outC_OOPS_wgt *= rewval; // Interactions everywhere else
        }

        wgt *= rewval;
        *iNode = true;
    }

    return checkWeight(wgt, interaction_chain);
}

double ReweightDriver::getWeightMIPP(IInteractionChainReweighting* reweighter,
                                     const InteractionChainData& icd,
                                     std::vector<bool>& interaction_nodes)
{
    if (!doMIPPNumi || m_hasMIPP) return 1.0;

    interaction_nodes = reweighter->canReweight(icd);
    for (const auto& node : interaction_nodes) {
        if (!node) continue;
        m_hasMIPP = true;
        return checkWeight(reweighter->calculateWeight(icd), icd.interaction_chain);
    }
    return 1.0;
}

double ReweightDriver::getWeightAttenuation(IInteractionChainReweighting* reweighter,
                                            const InteractionChainData& icd)
{
    const auto& nodes = reweighter->canReweight(icd);
    // we just see for the first position (primary proton)
    if (nodes.empty() || !nodes[0]) return 1.0;
    return checkWeight(reweighter->calculateWeight(icd), icd.interaction_chain);
}

ReweightDriver::~ReweightDriver()
{
    if (doMIPPNumi) {
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
