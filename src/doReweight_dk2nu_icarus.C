#include "MakeReweight.h"
#include "NuWeight.h"
#include "PDGParticleCodes.h"

#include "TH1D.h"
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>


// PDG codes
#define PI0 111
#define PIP 211
#define PIM -211
#define K0L 130
#define K0S 310
#define KP 321
#define KM -321
#define MUP -13
#define MUM 13

// const bool USE_NUCLEUS = true;

// Some constants::

static constexpr int NbinsE = 200;
static constexpr double emin = 0.;
static constexpr double emax = 20.;
static constexpr int Nnuhel = 4;
static const char* nuhel[Nnuhel] = {"numu", "numubar", "nue", "nuebar"};
static const char* nulabel[Nnuhel] = {"#nu_{#mu}", "#bar#nu_{#mu}", "#nu_{e}", "#bar#nu_{e}"};
static constexpr int Nparent = 5;
static const char* parent_name[Nparent] = {"", "pipm", "kpm", "k0l", "mu"};
static const char* parent_label[Nparent] = {"", "#pi^{#pm}", "K^{#pm}", "K^{0}_{L}", "#mu^{#pm}"};

using namespace NeutrinoFluxReweight;

class nu_g4numi;
class Numi2Pdg;

int idx_hel(int pdgdcode);

/*!
 * Run the reweighting for a single file (inputFile) for
 * particular MIPP covariance matrix given in input.xml file.
 */
void doReweight_dk2nu(const char* inputFile, const char* outputFile, const char* optionsFile,
                      const char* cxxdet, const char* cyydet, const char* czzdet) {
    TH1::SetDefaultSumw2();

    // const char* thisDir = getenv("PPFX_DIR");
    int idet = -1;
    // bool doing_precalculated_pos = false;

    // NeutrinoFluxAuxiliar::NuWeight* nuweight;
    const bool doing_precalculated_pos =
        std::string(cyydet) == "none" && std::string(czzdet) == "none";
    if (doing_precalculated_pos) {
        idet = atoi(cxxdet);
    }
    // else{
    //   std::vector<double> vdet;
    //   double xxdet = atof(cxxdet); vdet.push_back(xxdet);
    //   double yydet = atof(cyydet); vdet.push_back(yydet);
    //   double zzdet = atof(czzdet); vdet.push_back(zzdet);
    // }

    std::cout << "Instance of MakeReweight()" << std::endl;
    MakeReweight* makerew = MakeReweight::getInstance();
    makerew->SetOptions(optionsFile);

    std::cout << "Making an output file to store histograms" << std::endl;
    TFile fOut(outputFile, "recreate");
    std::cout << "File name: " << fOut.GetName() << std::endl;

    const int Nuniverses = makerew->GetNumberOfUniversesUsed();

    std::vector<std::vector<TH1D>> hnom(Nnuhel, std::vector<TH1D>(Nparent));
    std::vector<TH1D> hcv(Nnuhel);
    std::vector<std::vector<TH1D>> hthin(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hmipp(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hatt(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hothers(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> htotal(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_pCpi(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_pCk(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nCpi(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_pCnu(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_incoming_pip(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_incoming_pim(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_incoming_Kp(Nnuhel,  std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_incoming_Km(Nnuhel,  std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_incoming_K0(Nnuhel,  std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_outgoing_pip(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_outgoing_pim(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_outgoing_Kp(Nnuhel , std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_outgoing_Km(Nnuhel , std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_mesinc_outgoing_K0(Nnuhel , std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nua(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nua_inC_inPS(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nua_inC_OOPS(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nua_outC_Ascale(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nua_outC_OOPS(Nnuhel, std::vector<TH1D>(Nuniverses));
    std::vector<std::vector<TH1D>> hthin_nua_other(Nnuhel, std::vector<TH1D>(Nuniverses));


    auto initHist = [&](const char* name, const char* title) {
        const char* xtitle = "E_{#nu} [GeV]";
        const char* ytitle = "##nu_{unoscillated} [m^{-2}]";
        return TH1D(name, Form("%s;%s;%s", title, xtitle, ytitle), NbinsE, emin, emax);
    };

    for (int ii = 0; ii < Nnuhel; ii++) {
        // auto const& nu = nulabel[ii];
        for (int kk = 0; kk < Nparent; kk++) {
            if (kk == 0)
                hnom[ii][kk] = initHist(Form("hnom_%s", nuhel[ii]),
                                        Form("Uncorrected %s flux", nulabel[ii]));
            else
                hnom[ii][kk] = initHist(Form("hnom_%s_%s", nuhel[ii], parent_name[kk]),
                                        Form("Uncorrected flux of %s from %s decays", nulabel[ii], parent_label[kk]));
        }

        hcv[ii] = initHist(Form("hcv_%s", nuhel[ii]),
                           Form("Fully PPFX-corrected %s flux (central value)", nulabel[ii]));

        for (int jj = 0; jj < Nuniverses; jj++) {
            hthin[ii][jj]        = initHist(Form("hthin_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on thin target data, univ. #%i", nulabel[ii], jj));
            hmipp[ii][jj]        = initHist(Form("hmipp_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on MIPP NuMI target data, univ. #%i", nulabel[ii], jj));
            hatt[ii][jj]         = initHist(Form("hatt_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected for attenuation in target, univ. #%i", nulabel[ii], jj));
            hothers[ii][jj]      = initHist(Form("hothers_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected for other effects, univ. #%i", nulabel[ii], jj));
            htotal[ii][jj]       = initHist(Form("htotal_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected for all effects, univ. #%i", nulabel[ii], jj));
            hthin_pCpi[ii][jj]   = initHist(Form("hthin_pCpi_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on thin target p+C#rightarrow#pi+X data, univ. #%i", nulabel[ii], jj));
            hthin_pCk[ii][jj]    = initHist(Form("hthin_pCk_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on thin target p+C#rightarrowK+X data, univ. #%i", nulabel[ii], jj));
            hthin_nCpi[ii][jj]   = initHist(Form("hthin_nCpi_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on thin target n+C#rightarrow#pi+X data, univ. #%i", nulabel[ii], jj));
            hthin_pCnu[ii][jj]   = initHist(Form("hthin_pCnu_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on thin target p+C#rightarrowN+X data, univ. #%i", nulabel[ii], jj));
            hthin_mesinc[ii][jj] = initHist(Form("hthin_mesinc_%s_%d", nuhel[ii], jj),
                                            Form("%s flux corrected based on thin target data on meson interaction, univ. #%i", nulabel[ii], jj));

            hthin_mesinc_incoming_pip[ii][jj] = initHist(Form("hthin_mesinc_incoming_pip_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on #pi^{+} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_incoming_pim[ii][jj] = initHist(Form("hthin_mesinc_incoming_pim_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on #pi^{-} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_incoming_Kp[ii][jj]  = initHist(Form("hthin_mesinc_incoming_Kp_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on K^{+} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_incoming_Km[ii][jj]  = initHist(Form("hthin_mesinc_incoming_Km_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on K^{-} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_incoming_K0[ii][jj]  = initHist(Form("hthin_mesinc_incoming_K0_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on K^{0} interaction, univ. #%i", nulabel[ii], jj));

            hthin_mesinc_outgoing_pip[ii][jj] = initHist(Form("hthin_mesinc_outgoing_pip_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on #pi^{+} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_outgoing_pim[ii][jj] = initHist(Form("hthin_mesinc_outgoing_pim_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on #pi^{-} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_outgoing_Kp[ii][jj]  = initHist(Form("hthin_mesinc_outgoing_Kp_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on K^{+} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_outgoing_Km[ii][jj]  = initHist(Form("hthin_mesinc_outgoing_Km_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on K^{-} interaction, univ. #%i", nulabel[ii], jj));
            hthin_mesinc_outgoing_K0[ii][jj]  = initHist(Form("hthin_mesinc_outgoing_K0_%s_%d", nuhel[ii], jj),
                                                         Form("%s flux corrected based on sadly nonexistent thin target data on K^{0} interaction, univ. #%i", nulabel[ii], jj));

            hthin_nua[ii][jj]             = initHist(Form("hthin_nua_%s_%d", nuhel[ii], jj),
                                                     Form("%s flux corrected based on thin target data on N+detector material, univ. #%i", nulabel[ii], jj));
            hthin_nua_inC_inPS[ii][jj]    = initHist(Form("hthin_nua_inC_inPS_%s_%d", nuhel[ii], jj),
                                                     Form("%s flux corrected based on thin target data on N+detector material, univ. #%i", nulabel[ii], jj));
            hthin_nua_inC_OOPS[ii][jj]    = initHist(Form("hthin_nua_inC_OOPS_%s_%d", nuhel[ii], jj),
                                                     Form("%s flux corrected based on thin target data on N+detector material, univ. #%i", nulabel[ii], jj));
            hthin_nua_outC_Ascale[ii][jj] = initHist(Form("hthin_nua_outC_Ascale_%s_%d", nuhel[ii], jj),
                                                     Form("%s flux corrected based on thin target data on N+detector material, univ. #%i", nulabel[ii], jj));
            hthin_nua_outC_OOPS[ii][jj]   = initHist(Form("hthin_nua_outC_OOPS_%s_%d", nuhel[ii], jj),
                                                     Form("%s flux corrected based on thin target data on N+detector material, univ. #%i", nulabel[ii], jj));
            hthin_nua_other[ii][jj]       = initHist(Form("hthin_nua_other_%s_%d", nuhel[ii], jj),
                                                     Form("%s flux corrected based on thin target data on N+detector material, univ. #%i", nulabel[ii], jj));

            // hpCQEL[ii][jj]    = new TH1D(Form("hpCQEL_%s_%d",    nuhel[ii],jj),Form("%s flux
            // affected by p+C QEL: N#rightarrowN+X at x_{F}>0.95 or x_{F} > p_{T}/(GeV/c) +
            // 0.5, univ. #%i", nulabel[ii], jj), NbinsE,emin,emax);
            // hpCQEL[ii][jj]->SetXTitle(xtitle);
            // hpCQEL[ii][jj]->SetYTitle(ytitle);

            // hthin_nuAlFe[ii][jj]    = new TH1D(Form("hthin_nuAlFe_%s_%d",
            // nuhel[ii],jj),Form("%s flux affected by N+Al and N+Fe interactions, univ. #%i",
            // nulabel[ii], jj), NbinsE,emin,emax); hthin_nuAlFe[ii][jj]->SetXTitle(xtitle);
            // hthin_nuAlFe[ii][jj]->SetYTitle(ytitle);
        }
    }

    // Loading ntuples:
    TChain* chain_evts = new TChain("dk2nuTree");
    TChain* chain_meta = new TChain("dkmetaTree");
    bsim::Dk2Nu* dk2nu = new bsim::Dk2Nu;
    bsim::DkMeta* dkmeta = new bsim::DkMeta;

    std::cout << " Adding ntuple at: " << inputFile << std::endl;

    chain_evts->Add(inputFile);
    chain_evts->SetBranchAddress("dk2nu", &dk2nu);
    int nentries = chain_evts->GetEntries();

    chain_meta->Add(inputFile);
    chain_meta->SetBranchAddress("dkmeta", &dkmeta);
    chain_meta->GetEntry(0); // all entries are the same

    std::string detname = "UserPosition";
    if (doing_precalculated_pos) {
        detname = (dkmeta->location)[idet].name;
    }
    std::cout << "=> Doing the analysis for: " << detname << std::endl;

    std::cout << "#POT = " << dkmeta->pots << std::endl;
    TH1D hpot("hpot", "", 1, 0.0, 1.0);
    hpot.SetYTitle("#POT");
    hpot.SetBinContent(1, dkmeta->pots);

    std::cout << "N of entries: " << nentries << std::endl;
    double fluxWGT = 0.0;
    double nuenergy = 0.0;

    for (int ii = 0; ii < nentries; ii++) {
        if (ii % 1000 == 0) {
            std::cout << ii / 1000 << " k / " << nentries / 1000 << " k evts" << std::endl;
        }

        chain_evts->GetEntry(ii);
        const int nuidx = idx_hel(dk2nu->decay.ntype);

        /**
         * Cut on parents
         */
        // std::cout << "Here0\n";
        makerew->calculateWeights(dk2nu, dkmeta);
        // std::cout << "Here0.5\n";
        if (doing_precalculated_pos) {
            fluxWGT = ((dk2nu->nuray)[idet].wgt) * (dk2nu->decay.nimpwt) / 3.1416;
            nuenergy = (dk2nu->nuray)[idet].E;
        } else {
            std::vector<double> vdet = {atof(cxxdet), atof(cyydet), atof(czzdet)};
            auto const nuweight = std::make_unique<NeutrinoFluxAuxiliar::NuWeight>(vdet);
            nuweight->calculate_weight(dk2nu);
            fluxWGT = (nuweight->wgt) * (dk2nu->decay.nimpwt) / 3.1416;
            nuenergy = nuweight->enu;
        }

        // std::cout << "fluxWGT = " << fluxWGT << std::endl;
        // std::cout << "nuenergy = " << nuenergy << std::endl;
        // std::cout << std::endl;

        if (nuidx < 0) {
            std::cout << "=> Wrong neutrino file" << std::endl;
        }
        const int parent_id = dk2nu->decay.ptype;
        switch (parent_id) {
        case PIP:
        case PIM:
            hnom[nuidx][1].Fill(nuenergy, fluxWGT);
            break;
        case KP:
        case KM:
            hnom[nuidx][2].Fill(nuenergy, fluxWGT);
            break;
        case K0L:
            hnom[nuidx][3].Fill(nuenergy, fluxWGT);
            break;
        case MUP:
        case MUM:
            hnom[nuidx][4].Fill(nuenergy, fluxWGT);
            break;
        default:
            std::cout << "Unexpected neutrino parent id " << parent_id << std::endl;
        }
        hnom[nuidx][0].Fill(nuenergy, fluxWGT); // sum of all parents
        hcv[nuidx].Fill(nuenergy, fluxWGT * makerew->GetCVWeight());


        auto const vwgt_mipp_pi = makerew->GetWeights("MIPPNumiPionYields");
        auto const vwgt_mipp_K = makerew->GetWeights("MIPPNumiKaonYields");
        auto const vwgt_abs = makerew->GetWeights("TotalAbsorption");
        auto const vwgt_att = makerew->GetWeights("TargetAttenuation");
        auto const vwgt_ttpCpi = makerew->GetWeights("ThinTargetpCPion");
        auto const vwgt_ttpCk = makerew->GetWeights("ThinTargetpCKaon");
        auto const vwgt_ttnCpi = makerew->GetWeights("ThinTargetnCPion");
        auto const vwgt_ttpCnu = makerew->GetWeights("ThinTargetpCNucleon");
        auto const vwgt_ttmesinc = makerew->GetWeights("ThinTargetMesonIncident");
        auto const vwgt_ttmesinc_incoming_pip = makerew->GetWeights("ThinTargetMesonIncident_IncomingPip");
        auto const vwgt_ttmesinc_incoming_pim = makerew->GetWeights("ThinTargetMesonIncident_IncomingPim");
        auto const vwgt_ttmesinc_incoming_Kp = makerew->GetWeights("ThinTargetMesonIncident_IncomingKp");
        auto const vwgt_ttmesinc_incoming_Km = makerew->GetWeights("ThinTargetMesonIncident_IncomingKm");
        auto const vwgt_ttmesinc_incoming_K0 = makerew->GetWeights("ThinTargetMesonIncident_IncomingK0");
        auto const vwgt_ttmesinc_outgoing_pip = makerew->GetWeights("ThinTargetMesonIncident_OutgoingPip");
        auto const vwgt_ttmesinc_outgoing_pim = makerew->GetWeights("ThinTargetMesonIncident_OutgoingPim");
        auto const vwgt_ttmesinc_outgoing_Kp = makerew->GetWeights("ThinTargetMesonIncident_OutgoingKp");
        auto const vwgt_ttmesinc_outgoing_Km = makerew->GetWeights("ThinTargetMesonIncident_OutgoingKm");
        auto const vwgt_ttmesinc_outgoing_K0 = makerew->GetWeights("ThinTargetMesonIncident_OutgoingK0");
        auto const vwgt_ttnua = makerew->GetWeights("ThinTargetnucleonA");
        auto const vwgt_ttnua_inC_inPS = makerew->GetWeights("ThinTargetnucleonAInCInPS");
        auto const vwgt_ttnua_inC_OOPS = makerew->GetWeights("ThinTargetnucleonAInCOOPS");
        auto const vwgt_ttnua_outC_Ascale = makerew->GetWeights("ThinTargetnucleonAOutCAscale");
        auto const vwgt_ttnua_outC_OOPS = makerew->GetWeights("ThinTargetnucleonAOutCOOPS");
        auto const vwgt_ttnua_other = makerew->GetWeights("ThinTargetnucleonAOther");
        // auto const vwgt_ttnuAlFe = makerew->GetWeights("ThinTargetnucleonAlFe");
        auto const vwgt_oth = makerew->GetWeights("Other");
        // auto const vwgt_pCQEL  = makerew->GetWeights("ThinTargetpCQEL");

        for (int jj = 0; jj < Nuniverses; jj++) {
            const double wgt_thin = vwgt_ttpCpi[jj]
                                  * vwgt_ttpCk[jj]
                                  * vwgt_ttnCpi[jj]
                                  * vwgt_ttpCnu[jj]
                                  * vwgt_ttmesinc[jj]
                                  * vwgt_ttnua_inC_inPS[jj]
                                  * vwgt_ttnua_inC_OOPS[jj]
                                  * vwgt_ttnua_outC_Ascale[jj]
                                  * vwgt_ttnua_outC_OOPS[jj]
                                  * vwgt_ttnua_other[jj];
                                  // * vwgt_ttnua[jj]*vwgt_ttnuAlFe[jj];
            const double wgt_mipp = vwgt_mipp_pi[jj] * vwgt_mipp_K[jj];
            const double wgt_att = vwgt_att[jj] * vwgt_abs[jj];
            hthin[nuidx][jj].Fill(nuenergy, fluxWGT * wgt_thin);
            hmipp[nuidx][jj].Fill(nuenergy, fluxWGT * wgt_mipp);
            hatt[nuidx][jj].Fill(nuenergy, fluxWGT * wgt_att);
            hothers[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_oth[jj]);
            htotal[nuidx][jj].Fill(nuenergy, fluxWGT * wgt_thin * wgt_mipp * wgt_att * vwgt_oth[jj]);
            hthin_pCpi[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttpCpi[jj]);
            hthin_pCk[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttpCk[jj]);
            hthin_nCpi[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnCpi[jj]);
            hthin_pCnu[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttpCnu[jj]);
            hthin_mesinc[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttmesinc[jj]);
            hthin_mesinc_incoming_pip[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttmesinc_incoming_pip[jj]);
            hthin_mesinc_incoming_pim[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttmesinc_incoming_pim[jj]);
            hthin_mesinc_incoming_Kp[nuidx][jj].Fill(nuenergy,  fluxWGT * vwgt_ttmesinc_incoming_Kp[jj]);
            hthin_mesinc_incoming_Km[nuidx][jj].Fill(nuenergy,  fluxWGT * vwgt_ttmesinc_incoming_Km[jj]);
            hthin_mesinc_incoming_K0[nuidx][jj].Fill(nuenergy,  fluxWGT * vwgt_ttmesinc_incoming_K0[jj]);
            hthin_mesinc_outgoing_pip[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttmesinc_outgoing_pip[jj]);
            hthin_mesinc_outgoing_pim[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttmesinc_outgoing_pim[jj]);
            hthin_mesinc_outgoing_Kp[nuidx][jj].Fill(nuenergy,  fluxWGT * vwgt_ttmesinc_outgoing_Kp[jj]);
            hthin_mesinc_outgoing_Km[nuidx][jj].Fill(nuenergy,  fluxWGT * vwgt_ttmesinc_outgoing_Km[jj]);
            hthin_mesinc_outgoing_K0[nuidx][jj].Fill(nuenergy,  fluxWGT * vwgt_ttmesinc_outgoing_K0[jj]);
            hthin_nua[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnua[jj]);
            hthin_nua_inC_inPS[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnua_inC_inPS[jj]);
            hthin_nua_inC_OOPS[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnua_inC_OOPS[jj]);
            hthin_nua_outC_Ascale[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnua_outC_Ascale[jj]);
            hthin_nua_outC_OOPS[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnua_outC_OOPS[jj]);
            hthin_nua_other[nuidx][jj].Fill(nuenergy, fluxWGT * vwgt_ttnua_other[jj]);
            // hthin_nuAlFe[nuidx][jj]->Fill(nuenergy, fluxWGT * vwgt_ttnuAlFe[jj]);

            // hpCQEL [nuidx][jj]->Fill(nuenergy, fluxWGT * vwgt_pCQEL[jj]);
            // std::cout << "Here5\n";
        }
    }

    std::cout << "storing general histos" << std::endl;
    fOut.cd();

    fOut.mkdir("nom");
    fOut.mkdir("nom/parent");
    for (int ii = 0; ii < Nnuhel; ii++) {
        fOut.mkdir(Form("%s_thintarget", nuhel[ii]));
        fOut.mkdir(Form("%s_mippnumi", nuhel[ii]));
        fOut.mkdir(Form("%s_attenuation", nuhel[ii]));
        fOut.mkdir(Form("%s_others", nuhel[ii]));
        fOut.mkdir(Form("%s_total", nuhel[ii]));

        fOut.mkdir(Form("%s_thintarget/pCpi", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/pCk", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nCpi", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/pCnu", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_incoming_pip", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_incoming_pim", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_incoming_Kp", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_incoming_Km", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_incoming_K0", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_outgoing_pip", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_outgoing_pim", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_outgoing_Kp", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_outgoing_Km", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/mesinc_outgoing_K0", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nua", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nua_inC_inPS", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nua_inC_OOPS", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nua_outC_Ascale", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nua_outC_OOPS", nuhel[ii]));
        fOut.mkdir(Form("%s_thintarget/nua_other", nuhel[ii]));
        // fOut.mkdir(Form("%s_thintarget/pCfwd",  nuhel[ii]));
        // fOut.mkdir(Form("%s_thintarget/nuAlFe", nuhel[ii]));

        // fOut.mkdir(Form("%s_pCQEL",  nuhel[ii]));
    }

    for (int ii = 0; ii < Nnuhel; ii++) {
        fOut.cd("nom");
        hcv[ii].Write();
        hnom[ii][0].Write();
        fOut.cd("nom/parent");
        for (int kk = 1; kk < Nparent; kk++) {
            hnom[ii][kk].Write();
        }
        for (int jj = 0; jj < Nuniverses; jj++) {
            fOut.cd(Form("%s_thintarget", nuhel[ii]));
            hthin[ii][jj].Write();
            fOut.cd(Form("%s_mippnumi", nuhel[ii]));
            hmipp[ii][jj].Write();
            fOut.cd(Form("%s_attenuation", nuhel[ii]));
            hatt[ii][jj].Write();
            fOut.cd(Form("%s_others", nuhel[ii]));
            hothers[ii][jj].Write();
            fOut.cd(Form("%s_total", nuhel[ii]));
            htotal[ii][jj].Write();

            fOut.cd(Form("%s_thintarget/pCpi", nuhel[ii]));
            hthin_pCpi[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/pCk", nuhel[ii]));
            hthin_pCk[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nCpi", nuhel[ii]));
            hthin_nCpi[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/pCnu", nuhel[ii]));
            hthin_pCnu[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc", nuhel[ii]));
            hthin_mesinc[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_incoming_pip" , nuhel[ii]));
            hthin_mesinc_incoming_pip[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_incoming_pim" , nuhel[ii]));
            hthin_mesinc_incoming_pim[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_incoming_Kp" , nuhel[ii]));
            hthin_mesinc_incoming_Kp[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_incoming_Km" , nuhel[ii]));
            hthin_mesinc_incoming_Km[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_incoming_K0" , nuhel[ii]));
            hthin_mesinc_incoming_K0[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_outgoing_pip" , nuhel[ii]));
            hthin_mesinc_outgoing_pip[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_outgoing_pim" , nuhel[ii]));
            hthin_mesinc_outgoing_pim[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_outgoing_Kp" , nuhel[ii]));
            hthin_mesinc_outgoing_Kp[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_outgoing_Km" , nuhel[ii]));
            hthin_mesinc_outgoing_Km[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/mesinc_outgoing_K0" , nuhel[ii]));
            hthin_mesinc_outgoing_K0[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nua", nuhel[ii]));
            hthin_nua[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nua_inC_inPS", nuhel[ii]));
            hthin_nua_inC_inPS[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nua_inC_OOPS", nuhel[ii]));
            hthin_nua_inC_OOPS[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nua_outC_Ascale", nuhel[ii]));
            hthin_nua_outC_Ascale[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nua_outC_OOPS", nuhel[ii]));
            hthin_nua_outC_OOPS[ii][jj].Write();
            fOut.cd(Form("%s_thintarget/nua_other", nuhel[ii]));
            hthin_nua_other[ii][jj].Write();
            // fOut.cd(Form("%s_thintarget/nuAlFe" , nuhel[ii]));
            // hthin_nuAlFe[ii][jj]->Write();

            // fOut.cd(Form("%s_pCQEL"  , nuhel[ii])); hpCQEL [ii][jj]->Write();
        }
    }
    fOut.cd();
    hpot.Write();

    fOut.Close();

    // Releasing memory:
    makerew->resetInstance();

    std::cout << "End of run()" << std::endl;
}

int idx_hel(const int pdgcode) {
    switch (pdgcode) {
    case 14:
        return 0;
    case -14:
        return 1;
    case 12:
        return 2;
    case -12:
        return 3;
    default:
        return -1;
    }
}

void usage() {
    std::cout << "This script calculates the flux at one position in the NuMI beamline: "
              << std::endl;
    std::cout << "Using a precalculated detector positions:" << std::endl;
    std::cout << "  bin/doReweight_dk2nu [inputFile] [outputFile] [optionsFile] [idet]"
              << std::endl;
    std::cout << "Using a user input position:" << std::endl;
    std::cout
        << "  bin/doReweight_dk2nu [inputFile] [outputFile] [optionsFile] [xpos] [ypos] [zpos]"
        << std::endl;
    std::cout << "  " << std::endl;
    std::cout << "Inputs  " << std::endl;
    std::cout << "[inputFile] : g4numi ntuple in dk2nu/dkmeta format (v6 minerva branch is "
                 "recommended)"
              << std::endl;
    std::cout << "[outputFile] : user definied output file name." << std::endl;
    std::cout << "[optionsFile] :xml file with the ppfx input parameters (look at "
                 "${PPFX_DIR}/script/input_default.xml)"
              << std::endl;
    std::cout << "[idet] : index of the precalculated detector (look at the location.name in "
                 "the dkmeta tree of the g4numi ntuple)"
              << std::endl;
    std::cout << "[xpos], [ypos], [zpos] : position (cm) respect to the MC NuMI coordinate "
                 "system to calculate the flux"
              << std::endl;
    std::cout << "  " << std::endl;
}
//////////////

#ifndef __CINT__
int main(int argc, const char* argv[]) {
    if (argc == 5) {
        doReweight_dk2nu(argv[1], argv[2], argv[3], argv[4], "none", "none");
    } else if (argc == 7) {
        doReweight_dk2nu(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
    } else {
        std::cout << "Error: Invalid number of arguments (" << (argc - 1) << ")" << std::endl;
        usage();
        exit(1);
    }

    return 0;
}
#endif
