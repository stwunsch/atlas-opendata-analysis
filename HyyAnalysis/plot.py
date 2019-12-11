import ROOT
ROOT.gROOT.SetBatch(True)
import argparse


def main(path):
    # Set styles
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetEndErrorSize(0)
    ROOT.gStyle.SetErrorX(0)

    # Create canvas with pads for main plot and data/MC ratio
    c = ROOT.TCanvas("c", "", 700, 750)

    upper_pad= ROOT.TPad("upper_pad", "", 0, 0.29, 1, 1, 0, 0, 0)
    lower_pad = ROOT.TPad("lower_pad", "", 0, 0, 1, 0.29, 0, 0, 0)
    for p in [upper_pad, lower_pad]:
        p.SetLeftMargin(0.14)
        p.SetRightMargin(0.05)

    upper_pad.Draw()
    lower_pad.Draw()

    # Load histograms from file
    f = ROOT.TFile(path, "READ")
    ggh = f.Get("ggH")
    vbf = f.Get("VBF")
    data = f.Get("data")

    # Fit signal + background model to data
    upper_pad.cd()
    fit = ROOT.TF1("fit", "([0]+[1]*x+[2]*x^2+[3]*x^3)+[4]*exp(-0.5*((x-[5])/[6])^2)", 105, 160)
    fit.FixParameter(5, 125.0)
    fit.FixParameter(4, 119.1)
    fit.FixParameter(6, 2.39)
    data.Fit("fit", "", "same:e", 105, 160)
    fit.SetLineColor(2)
    fit.SetLineStyle(1)
    fit.SetLineWidth(2)
    fit.Draw("SAME")

    # Draw background
    bkg = ROOT.TF1("bgd", "([0]+[1]*x+[2]*x^2+[3]*x^3)", 105, 160)
    for i in range(4):
        bkg.SetParameter(i, fit.GetParameter(i))
    bkg.SetLineColor(4)
    bkg.SetLineStyle(2)
    bkg.SetLineWidth(2)
    bkg.Draw("SAME")

    # Draw data
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E1 SAME")
    data.SetMinimum(0)
    data.SetMaximum(16e3)

    # Scale simulated events with luminosity * cross-section / sum of weights
    # and merge to single Higgs signal
    lumi = 10064.0
    ggh.Scale(lumi * 0.102 / ggh.Integral())
    vbf.Scale(lumi * 0.008518764 / vbf.Integral())
    higgs = ggh
    higgs.Add(vbf)
    higgs.Draw("HIST SAME")

    # Add legend
    legend = ROOT.TLegend(0.60, 0.55, 0.89, 0.85)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.05)
    legend.SetTextAlign(32)
    legend.AddEntry(data, "Data" ,"lep")
    legend.AddEntry(bkg, "Background", "l")
    legend.AddEntry(fit, "Signal + Bkg.", "l")
    legend.AddEntry(higgs, "Signal", "l")
    legend.Draw("SAME")

    # Add ATLAS label
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(72)
    text.SetTextSize(0.05)
    text.DrawLatex(0.18, 0.84, "ATLAS")

    text.SetTextFont(42)
    text.DrawLatex(0.18 + 0.13, 0.84, "Open Data")

    text.SetTextSize(0.04)
    text.DrawLatex(0.18, 0.78, "#sqrt{s} = 13 TeV, 10 fb^{-1}")

    # Plot ratio
    # FIXME
    lower_pad.cd()
    ratio = data.Clone()
    ratio.Add(bkg, -1.0)
    ratio.Draw("HIST SAME")
    ratio.SetMinimum(0)
    ratio.SetMaximum(500)
    higgs.Draw("HIST SAME")

    # Save plot to file
    c.SaveAs("HyyAnalysis.png")
    c.SaveAs("HyyAnalysis.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Path to ROOT file with histograms")
    args = parser.parse_args()
    main(args.path)
