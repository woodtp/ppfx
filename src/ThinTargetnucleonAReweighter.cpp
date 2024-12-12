
#include "ThinTargetnucleonAReweighter.h"

#include <cstdlib>
#include <cwctype>
#include <stdexcept>

#include "MakeReweight.h"
#include "PDGParticleCodes.h"
#include "ThinTargetMC.h"

namespace NeutrinoFluxReweight {

ThinTargetnucleonAReweighter::ThinTargetnucleonAReweighter(int iuniv,
                                                           const ParameterTable& cv_pars,
                                                           const ParameterTable& univ_pars)
  : iUniv(iuniv)
  , cvPars(cv_pars)
  , univPars(univ_pars)
  , m_scalingBinID(-1)
{
    const std::string mode(getenv("MODE"));
    m_isRefOpt = "REF" == mode || "OPT" == mode;

    m_thinBins = ThinTargetBins::getInstance();

    vbin_data_pip.reserve(m_thinBins->GetNbins_material_scaling());
    vbin_data_pim.reserve(m_thinBins->GetNbins_material_scaling());
    vbin_data_kap.reserve(m_thinBins->GetNbins_material_scaling());
    vbin_data_kam.reserve(m_thinBins->GetNbins_material_scaling());

    // Currently, We are using the same number of xF ranges for nucleon inc. and meson inc.
    vbin_prt_inc_pip.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_prt_inc_pim.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_prt_inc_kap.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_prt_inc_kam.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_prt_inc_k0.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_prt_inc_p.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_prt_inc_n.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_pip.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_pim.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_kap.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_kam.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_k0.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_p.reserve(m_thinBins->GetNbins_meson_incident());
    vbin_neu_inc_n.reserve(m_thinBins->GetNbins_meson_incident());

    char namepar[100];

    data_prod_xs = univPars.getParameterValue("prod_prtC_xsec");

    // 4 particles
    const char* cinc[4] = { "pip", "pim", "kap", "kam" };
    for (int ii = 0; ii < 4; ii++) {
        for (int jj = 0; jj < m_thinBins->GetNbins_material_scaling(); jj++) {
            sprintf(namepar, "ThinTarget_material_scaling_%s_%d", cinc[ii], jj);
            double dataval = univPars.getParameterValue(std::string(namepar));
            if (ii == 0)
                vbin_data_pip.push_back(dataval);
            if (ii == 1)
                vbin_data_pim.push_back(dataval);
            if (ii == 2)
                vbin_data_kap.push_back(dataval);
            if (ii == 3)
                vbin_data_kam.push_back(dataval);
        }
    }

    // for all nucleons incident not covered by any thin target reweighters
    // 2 incident nucleons, 7 produced particles:
    const char* nuinc[2] = { "prt", "neu" };
    const char* cpro[7] = { "pip", "pim", "kap", "kam", "k0", "n", "p" };

    for (int ii = 0; ii < 2; ii++) {
        for (int jj = 0; jj < 7; jj++) {
            for (int kk = 0; kk < m_thinBins->GetNbins_meson_incident(); kk++) {
                sprintf(namepar, "ThinTarget_%s_incident_%s_%d", nuinc[ii], cpro[jj], kk);
                double dataval = univPars.getParameterValue(std::string(namepar));
                if (ii == 0 && jj == 0)
                    vbin_prt_inc_pip.push_back(dataval);
                if (ii == 0 && jj == 1)
                    vbin_prt_inc_pim.push_back(dataval);
                if (ii == 0 && jj == 2)
                    vbin_prt_inc_kap.push_back(dataval);
                if (ii == 0 && jj == 3)
                    vbin_prt_inc_kam.push_back(dataval);
                if (ii == 0 && jj == 4)
                    vbin_prt_inc_k0.push_back(dataval);
                if (ii == 0 && jj == 5)
                    vbin_prt_inc_n.push_back(dataval);
                if (ii == 0 && jj == 6)
                    vbin_prt_inc_p.push_back(dataval);
                if (ii == 1 && jj == 0)
                    vbin_neu_inc_pip.push_back(dataval);
                if (ii == 1 && jj == 1)
                    vbin_neu_inc_pim.push_back(dataval);
                if (ii == 1 && jj == 2)
                    vbin_neu_inc_kap.push_back(dataval);
                if (ii == 1 && jj == 3)
                    vbin_neu_inc_kam.push_back(dataval);
                if (ii == 1 && jj == 4)
                    vbin_neu_inc_k0.push_back(dataval);
                if (ii == 1 && jj == 5)
                    vbin_neu_inc_n.push_back(dataval);
                if (ii == 1 && jj == 6)
                    vbin_neu_inc_p.push_back(dataval);
            }
        }
    }
    // left over:
    sprintf(namepar, "ThinTarget_prtleftover_incident_%d", 0);
    bin_prtleftover_inc = univPars.getParameterValue(std::string(namepar));
    sprintf(namepar, "ThinTarget_neuleftover_incident_%d", 0);
    bin_neuleftover_inc = univPars.getParameterValue(std::string(namepar));
}

ThinTargetnucleonAReweighter::~ThinTargetnucleonAReweighter() {}

bool ThinTargetnucleonAReweighter::canReweight(const InteractionData& aa)
{
    return pdg::P == aa.Inc_pdg || pdg::N == aa.Inc_pdg;
}

bool ThinTargetnucleonAReweighter::isDataBased(const InteractionData& aa)
{
    m_scalingBinID = m_thinBins->material_scaling_BinID(aa.xF, aa.Pt, aa.Prod_pdg);

    const bool is_meson_production = pdg::PIP == aa.Prod_pdg
                                  || pdg::PIM == aa.Prod_pdg
                                  || pdg::KP  == aa.Prod_pdg
                                  || pdg::KM  == aa.Prod_pdg
                                  || pdg::K0S == aa.Prod_pdg
                                  || pdg::K0L == aa.Prod_pdg;

    const bool is_data_volume =
      m_isRefOpt ? aa.Vol == "TargetNoSplitSegment" && aa.Vol == "TargetFinHorizontal"
                 : aa.Vol == "TGT1" || aa.Vol == "BudalMonitor" || aa.Vol == "Budal_HFVS" ||
                    aa.Vol == "Budal_VFHS";

    return aa.Inc_P >= 12.0 && !is_data_volume && is_meson_production && m_scalingBinID >= 0;
}

double ThinTargetnucleonAReweighter::calculateWeight(const InteractionData& aa)
{
    if (!canReweight(aa)) {
        const std::string msg = "ThinTargetnucleonAReweighter::calculateWeight() - cannot "
                                "reweight this interaction, incident pdg = " +
                                std::to_string(aa.Inc_pdg);
        throw std::invalid_argument(msg);
    }

    double wgt = 1.0;

    const std::string tgtent = m_isRefOpt ? "TargetFinHorizontal" : "TGT1";

    const double inc_mom[3] = { aa.Inc_P4[0], aa.Inc_P4[1], aa.Inc_P4[2] };
    const double prod_mom[3] = { aa.Prod_P4[0], aa.Prod_P4[1], aa.Prod_P4[2] };
    const double vtx_int[3] = { aa.Vtx[0], aa.Vtx[1], aa.Vtx[2] };

    InteractionData aux_aa2(
      aa.gen, inc_mom, aa.Inc_pdg, prod_mom, aa.Prod_pdg, tgtent, aa.nucleus, aa.Proc, vtx_int);

    bool not_handled = false;
    if (isDataBased(aa)) {
        ThinTargetMC* mc = ThinTargetMC::getInstance();
        const double mc_prod = mc->getMCxs_pC_piK(0, aa.Inc_P);
        const double fact_gen0 = data_prod_xs / mc_prod;
        MakeReweight* makerew = MakeReweight::getInstance();

        if (pdg::P == aa.Inc_pdg) {
            if (pdg::PIP == aa.Prod_pdg || pdg::PIM == aa.Prod_pdg) {
                if (-1 == iUniv)
                    tt_pCPionRew = (makerew->cv_rw)->THINTARGET_PC_PION_Universe;
                else
                    tt_pCPionRew = (makerew->vec_rws[iUniv])->THINTARGET_PC_PION_Universe;

                if (tt_pCPionRew->canReweight(aux_aa2)) {
                    wgt = tt_pCPionRew->calculateWeight(aux_aa2);
                    if (0 == aux_aa2.gen)
                        wgt *= fact_gen0;
                }
                else
                    not_handled = true;
            }
            else if (pdg::KP == aa.Prod_pdg || pdg::KM == aa.Prod_pdg ||
                     pdg::K0L == aa.Prod_pdg || pdg::K0S == aa.Prod_pdg) {
                if (-1 == iUniv)
                    tt_pCKaonRew = (makerew->cv_rw)->THINTARGET_PC_KAON_Universe;
                else
                    tt_pCKaonRew = (makerew->vec_rws[iUniv])->THINTARGET_PC_KAON_Universe;

                if (tt_pCKaonRew->canReweight(aux_aa2)) {
                    wgt = tt_pCKaonRew->calculateWeight(aux_aa2);
                    if (0 == aux_aa2.gen)
                        wgt *= fact_gen0;
                }
                else
                    not_handled = true;
            }
            else
                not_handled = true;
        }
        else if (pdg::N == aa.Inc_pdg) {
            if (pdg::PIP == aa.Prod_pdg || pdg::PIM == aa.Prod_pdg) {
                if (-1 == iUniv)
                    tt_nCPionRew = (makerew->cv_rw)->THINTARGET_NC_PION_Universe;
                else
                    tt_nCPionRew = (makerew->vec_rws[iUniv])->THINTARGET_NC_PION_Universe;

                if (tt_nCPionRew->canReweight(aux_aa2)) {
                    wgt = tt_nCPionRew->calculateWeight(aux_aa2);
                    if (0 == aux_aa2.gen)
                        wgt *= fact_gen0;
                }
                else
                    not_handled = true;
            }
            else
                not_handled = true;
        }
        else
            not_handled = true;

        double scaling = 1.0;
        if (pdg::PIP == aa.Prod_pdg)
            scaling = vbin_data_pip[m_scalingBinID];
        else if (pdg::PIM == aa.Prod_pdg)
            scaling = vbin_data_pim[m_scalingBinID];
        else if (pdg::KP == aa.Prod_pdg)
            scaling = vbin_data_kap[m_scalingBinID];
        else if (pdg::KM == aa.Prod_pdg)
            scaling = vbin_data_kam[m_scalingBinID];
        else if (pdg::K0S == aa.Prod_pdg || pdg::K0L == aa.Prod_pdg)
            scaling = vbin_data_kap[m_scalingBinID];
        wgt *= scaling;
        if (!not_handled)
            return wgt;
    } // if(is_data_based)

    // trick... using a function for meson incident... same binning.
    int binnu = m_thinBins->meson_inc_BinID(aa.xF, aa.Pt, pdg::PIP);
    // we've modified the meson bins to include negative xF, so these are unphysical bins
    // but don't return 1. either way
    if (binnu < 0) {
        if (pdg::P == aa.Inc_pdg)
            return bin_prtleftover_inc;
        else if (pdg::N == aa.Inc_pdg)
            return bin_neuleftover_inc;
    }

    if (pdg::P == aa.Inc_pdg) {
        // add extra uncertainties for xF < 0
        // treatment here is basically 40% corr. across hadron species + 40% uncorr. across all
        // species, xF bins
#ifdef UH_ICARUS
#pragma message "[INFO] UH_ICARUS is defined. Disabling 40% correlation for negative xF."
        // A. Wood (apwood@central.uh.edu)
        // 2024-12-12
        // We found the 40% correlations to be a bit too aggressive for ICARUS.
        // See SBN DocDB 38113 for more details.
        const double negxF_corrunc =  1.0;
#else
        const double negxF_corrunc =  aa.xF < 0.0 ? bin_prtleftover_inc : 1.0;
#endif

        if (pdg::PIP == aa.Prod_pdg)
            wgt = vbin_prt_inc_pip[binnu] * negxF_corrunc;
        else if (pdg::PIM == aa.Prod_pdg)
            wgt = vbin_prt_inc_pim[binnu] * negxF_corrunc;
        else if (pdg::KP == aa.Prod_pdg)
            wgt = vbin_prt_inc_kap[binnu] * negxF_corrunc;
        else if (pdg::KM == aa.Prod_pdg)
            wgt = vbin_prt_inc_kam[binnu] * negxF_corrunc;
        else if (pdg::K0L == aa.Prod_pdg || pdg::K0S == aa.Prod_pdg)
            wgt = vbin_prt_inc_k0[binnu] * negxF_corrunc;
        else if (pdg::P == aa.Prod_pdg)
            wgt = vbin_prt_inc_p[binnu] * negxF_corrunc;
        else if (pdg::N == aa.Prod_pdg)
            wgt = vbin_prt_inc_n[binnu] * negxF_corrunc;
        else
            wgt = bin_prtleftover_inc;
    }
    else if (pdg::N == aa.Inc_pdg) {
        const double negxF_corrunc = aa.xF < 0.0 ? bin_neuleftover_inc : 1.0;

        if (pdg::PIP == aa.Prod_pdg)
            wgt = vbin_neu_inc_pip[binnu] * negxF_corrunc;
        else if (pdg::PIM == aa.Prod_pdg)
            wgt = vbin_neu_inc_pim[binnu] * negxF_corrunc;
        else if (pdg::KP == aa.Prod_pdg)
            wgt = vbin_neu_inc_kap[binnu] * negxF_corrunc;
        else if (pdg::KM == aa.Prod_pdg)
            wgt = vbin_neu_inc_kam[binnu] * negxF_corrunc;
        else if (pdg::K0L == aa.Prod_pdg || pdg::K0S == aa.Prod_pdg)
            wgt = vbin_neu_inc_k0[binnu] * negxF_corrunc;
        else if (pdg::P == aa.Prod_pdg)
            wgt = vbin_neu_inc_p[binnu] * negxF_corrunc;
        else if (pdg::N == aa.Prod_pdg)
            wgt = vbin_neu_inc_n[binnu] * negxF_corrunc;
        else
            wgt = bin_neuleftover_inc;
    }

    if (wgt < 0.)
        return 0.0001; // cap this at near-0 instead of returning 1.
    if (wgt > 10.)
        return 1.0; // ignore larger weights than this

    return wgt;
}
} // namespace NeutrinoFluxReweight
