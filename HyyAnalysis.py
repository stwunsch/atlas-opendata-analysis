import ROOT
ROOT.gROOT.SetBatch(True)
import argparse
import os


def get_data_samples(path):
    samples = ROOT.std.vector("string")()
    for tag in ["A", "B", "C", "D"]:
        samples.push_back(os.path.join(path, "GamGam/Data/data_{}.GamGam.root".format(tag)))
    return samples


def get_ggH125_samples(path):
    samples = ROOT.std.vector("string")()
    samples.push_back(os.path.join(path, "GamGam/MC/mc_343981.ggH125_gamgam.GamGam.root"))
    return samples


def get_VBF125_samples(path):
    samples = ROOT.std.vector("string")()
    samples.push_back(os.path.join(path, "GamGam/MC/mc_345041.VBFH125_gamgam.GamGam.root"))
    return samples


def main(path):
    # Run analysis using multi-threading
    ROOT.ROOT.EnableImplicitMT()

    # Set up the dataframes with the datasets corresponding to the respective process
    df = {}
    df["data"] = ROOT.RDataFrame("mini", get_data_samples(path))
    df["ggH"] = ROOT.RDataFrame("mini", get_ggH125_samples(path))
    df["VBF"] = ROOT.RDataFrame("mini", get_VBF125_samples(path))
    processes = list(df.keys())

    # Apply scale factors and MC weight for simulated events and a weight of 1 for the data
    for p in ["ggH", "VBF"]:
        df[p] = df[p].Define("weight",
                "scaleFactor_PHOTON * scaleFactor_PhotonTRIGGER * scaleFactor_PILEUP * mcWeight");
    df["data"] = df["data"].Define("weight", "1.0")

    # Make selection on the datasets
    for p in processes:
        # Apply preselection cut on photon trigger
        df[p] = df[p].Filter("trigP == true", "Preselection cut on photon trigger")

        # Find two good muons with pt > 25 GeV and not in the transition region between barrel and encap
        df[p] = df[p].Define("goodphotons",
                             "(photon_pt > 25000) && (abs(photon_eta) < 2.37) && ((abs(photon_eta) < 1.37) || (abs(photon_eta) > 1.52))")\
                     .Filter("Sum(goodphotons) == 2",
                             "Exactly two good photons with pt > 25 GeV and in good detector region")

        # Take only isolated photons
        df[p] = df[p].Filter("Sum(photon_ptcone30[goodphotons] / photon_pt[goodphotons] < 0.065) == 2",
                             "Events with two isolated photons (ptcone)")\
                     .Filter("Sum(photon_etcone20[goodphotons] / photon_pt[goodphotons] < 0.065) == 2",
                             "Events with two isolated photons (etcone)")

    # Make fourvectors for the two good photons and compute invariant mass for further analysis
    ROOT.gInterpreter.Declare(
    """
    using Vec_t = const ROOT::VecOps::RVec<float>;
    using Fourvec_t = ROOT::Math::PtEtaPhiEVector;
    std::vector<Fourvec_t> MakeFourvectors(Vec_t& pt, Vec_t& eta, Vec_t& phi, Vec_t& e) {
        std::vector<Fourvec_t> fourvecs;
        fourvecs.emplace_back(Fourvec_t(pt[0], eta[0], phi[0], e[0]));
        fourvecs.emplace_back(Fourvec_t(pt[1], eta[1], phi[1], e[1]));
        return fourvecs;
    }
    """)

    hists = {}
    for p in processes:
        # Make four vectors and compute invariant mass
        df[p] = df[p].Define("p4s", "MakeFourvectors(photon_pt[goodphotons], photon_eta[goodphotons], photon_phi[goodphotons], photon_E[goodphotons])")
        df[p] = df[p].Define("m_yy", "(p4s[0] + p4s[1]).mass() / 1000.0")

        # Make additional kinematic cuts
        df[p] = df[p].Filter("p4s[0].pt() / 1000.0 * m_yy > 0.35", "Kinematic cut on leading photon")\
                     .Filter("p4s[1].pt() / 1000.0 * m_yy > 0.25", "Kinematic cut on subleading photon")\
                     .Filter("(m_yy > 105) && (m_yy < 160)", "Select mass window between 105 GeV and 160 GeV")

        # Book histogram of the invariant mass with this selection
        hists[p] = df[p].Histo1D(
                ROOT.ROOT.RDF.TH1DModel(p, "Diphoton invariant mass; m_{#gamma#gamma} [GeV];Events / bin", 30, 105, 160),
                "m_yy", "weight")

        # Book histogram in unconverted central category
        df[p] = df[p].Filter("Sum((abs(photon_eta[goodphotons]) < 0.75) && (photon_convType[goodphotons] == 0)) == 2",
                             "Unconverted central category")
        hists[p + "_cat"] = df[p].Histo1D(
                ROOT.ROOT.RDF.TH1DModel(p + "_cat", "Diphoton invariant mass; m_{#gamma#gamma} [GeV];Events / bin", 30, 105, 160),
                "m_yy", "weight")

    # Book cut-flow reports
    reports = {}
    for p in processes:
        reports[p] = df[p].Report()

    # Write histograms to file
    f = ROOT.TFile("HyyAnalysis.root", "RECREATE")
    for h in hists:
        hists[h].Write()
    f.Close()

    # Print cut-flow reports
    for p in processes:
        print(">>> Cut-flow report for {}".format(p))
        reports[p].Print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to base directory with datasets")
    args = parser.parse_args()
    main(args.path)
