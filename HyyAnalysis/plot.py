import ROOT
ROOT.gROOT.SetBatch(True)
import argparse


def main(path):
    # Create canvas with pads for main plot and data/MC ratio
    c = ROOT.TCanvas("c", "", 700, 750)

    pad0 = ROOT.TPad("pad0", "", 0, 0.29, 1, 1, 0, 0, 0)
    pad0.SetTickx(False)
    pad0.SetTicky(False)
    pad0.SetTopMargin(0.05)
    pad0.SetBottomMargin(0)
    pad0.SetLeftMargin(0.14)
    pad0.SetRightMargin(0.05)
    pad0.SetFrameBorderMode(0)
    pad0.SetTopMargin(0.06)
    pad0.Draw()

    pad1 = ROOT.TPad("pad1", "", 0, 0, 1, 0.29, 0, 0, 0)
    pad1.SetTickx(False)
    pad1.SetTicky(False)
    pad1.SetTopMargin(0.0)
    pad1.SetBottomMargin(0.5)
    pad1.SetLeftMargin(0.14)
    pad1.SetRightMargin(0.05)
    pad1.SetFrameBorderMode(0)
    pad1.Draw()

    # Load histograms from file
    f = ROOT.TFile(path, "READ")
    ggh = f.Get("ggH")
    vbf = f.Get("VBF")
    data = f.Get("data")

    # Scale simulated events with luminosity * cross-section / sum of weights
    # and merge to single Higgs signal
    lumi = 10064.0
    ggh.Scale(lumi * 0.102 / ggh.Integral())
    vbf.Scale(lumi * 0.008518764 / vbf.Integral())
    higgs = ggh
    higgs.Add(vbf)

    # Save plot to file
    c.SaveAs("HyyAnalysis.png")
    c.SaveAs("HyyAnalysis.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to ROOT file with histograms")
    args = parser.parse_args()
    main(args.path)
