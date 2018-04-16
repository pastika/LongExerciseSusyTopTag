#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"

#include <memory>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>

//This is a helper function which will keep the plot from overlapping with the legend
void smartMax(const TH1 * const h, const TLegend* const l, const TPad* const p, double& gmin, double& gmax, double& gpThreshMax, const bool error)
{
    const bool isLog = p->GetLogy();
    double min = 9e99;
    double max = -9e99;
    double pThreshMax = -9e99;
    int threshold = static_cast<int>(h->GetNbinsX()*(l->GetX1() - p->GetLeftMargin())/((1 - p->GetRightMargin()) - p->GetLeftMargin()));

    for(int i = 1; i <= h->GetNbinsX(); ++i)
    {
        double bin = 0.0;
        if(error) bin = h->GetBinContent(i) + h->GetBinError(i);
        else      bin = h->GetBinContent(i);
        if(bin > max) max = bin;
        else if(bin > 1e-10 && bin < min) min = bin;
        if(i >= threshold && bin > pThreshMax) pThreshMax = bin;
    }

    gpThreshMax = std::max(gpThreshMax, pThreshMax);
    gmax = std::max(gmax, max);
    gmin = std::min(gmin, min);
}

//Class to hold TH1* with various helper functions 
class histInfo
{
public:
    std::string legEntry, histFile, histName, drawOptions;
    int color, rebin;
    std::shared_ptr<TH1> h;

    //helper function to get histogram from file and configure its optional settings
    void retrieveHistogram()
    {
        //Open the file for this histogram
        TFile *f = TFile::Open(histFile.c_str());

        //check that the file was opened successfully
        if(!f)
        {
            printf("File \"%s\" could not be opened!!!\n", histFile.c_str());
            h = nullptr;
            return;
        }

        //get the histogram from the file
        h.reset(static_cast<TH1*>(f->Get(histName.c_str())));

        //with the histogram retrieved, close the file
        f->Close();
        delete f;

        //check that the histogram was retireved from the file successfully
        if(!h)
        {
            printf("Histogram \"%s\" could not be found in file \"%s\"!!!\n", histName.c_str(), histFile.c_str());
            return;
        }

        //set the histogram color
        h->SetLineColor(color);
        h->SetLineWidth(3);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(20);

        // rebin the histogram if desired
        if(rebin > 0) h->Rebin(rebin);
    }

    //helper function for axes
    void setupAxes(double xOffset, double yOffset, double xTitle, double yTitle, double xLabel, double yLabel)
    {
        h->SetStats(0);
        h->SetTitle(0);
        h->GetXaxis()->SetTitleOffset(xOffset);
        h->GetYaxis()->SetTitleOffset(yOffset);
        h->GetXaxis()->SetTitleSize(xTitle);
        h->GetYaxis()->SetTitleSize(yTitle);
        h->GetXaxis()->SetLabelSize(xLabel);
        h->GetYaxis()->SetLabelSize(yLabel);
        if(h->GetXaxis()->GetNdivisions() % 100 > 5) h->GetXaxis()->SetNdivisions(6, 5, 0);
    }

    //helper function for pads
    void setupPad(double left, double right, double top, double bottom)
    {
        gPad->SetLeftMargin(left);
        gPad->SetRightMargin(right);
        gPad->SetTopMargin(top);
        gPad->SetBottomMargin(bottom);
        gPad->SetTicks(1,1);
    }

    //wrapper to draw histogram
    void draw(const std::string& additionalOptions = "", bool noSame = false) const
    {
        h->Draw(((noSame?"":"same " + drawOptions + " " + additionalOptions)).c_str());
    }

    void setFillColor(int newColor = -1)
    {
        if(newColor >= 0) h->SetFillColor(newColor);
        else              h->SetFillColor(color);
    }

    histInfo(const std::string& legEntry, const std::string& histFile, const std::string& drawOptions, const int color) : legEntry(legEntry), histFile(histFile), histName(""), drawOptions(drawOptions), color(color), rebin(-1), h(nullptr)
    {
    }

    histInfo(TH1* h) : legEntry(h->GetName()), histFile(""), histName(h->GetName()), drawOptions(""), color(0), rebin(0), h(h)
    {
    }

    ~histInfo()
    {
    }
};

class Plotter
{
private:
    //entry for data
    histInfo data_;
    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries_;
    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries_;
    
public:
    Plotter(histInfo&  data, std::vector<histInfo>&  bgEntries, std::vector<histInfo>&  sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}
    Plotter(histInfo&& data, std::vector<histInfo>&& bgEntries, std::vector<histInfo>&& sigEntries) : data_(data), bgEntries_(bgEntries), sigEntries_(sigEntries) {}

    void plot(const std::string& histName, const std::string& xAxisLabel, const std::string& yAxisLabel = "Events", const bool isLogY = false, const double xmin = 999.9, const double xmax = -999.9, int rebin = -1, double lumi = 36100)
    {
//        printf("Begin plotting histogram %s\n", histName.c_str());

        //This is a magic incantation to disassociate opened histograms from their files so the files can be closed
        TH1::AddDirectory(false);

        //create the canvas for the plot
        TCanvas *c = new TCanvas("c1", "c1", 800, 800);
        //switch to the canvas to ensure it is the active object
        c->cd();

        // Upper plot will be in pad1: TPad(x1, y1, x2, y2)
        TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
        //pad1->SetBottomMargin(0); // Upper and lower plot are joined
        pad1->SetGridy();         // Horizontal grid
        pad1->Draw();             // Draw the upper pad: pad1
        pad1->cd();               // pad1 becomes the current pad

        //Create TLegend: TLegend(x1, y1, x2, y2)
        TLegend *leg = new TLegend(0.50, 0.56, 0.89, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        //get maximum from histos and fill TLegend
        double min = 0.0;
        double max = 0.0;
        double lmax = 0.0;

        //Create the THStack for the background ... warning, THStacks are terrible and must be filled "backwards"
        THStack *bgStack = new THStack();
        //Make seperate histogram from sum of BG histograms because I don't know how to make a THStack give me this 
        TH1* hbgSum = nullptr;

        for(int iBG = bgEntries_.size() - 1; iBG >= 0; --iBG)
        {
            //Get new histogram
            bgEntries_[iBG].histName = histName;
            bgEntries_[iBG].rebin = rebin;
            bgEntries_[iBG].retrieveHistogram();

            bgStack->Add(bgEntries_[iBG].h.get(), bgEntries_[iBG].drawOptions.c_str());
            if(!hbgSum) hbgSum = static_cast<TH1*>(bgEntries_[iBG].h->Clone());
            else        hbgSum->Add(bgEntries_[iBG].h.get());
        }

        //data
        //get new histogram from file
        data_.histName = histName;
        data_.rebin = rebin;
        data_.retrieveHistogram();

        char label[256];
        sprintf(label, "%s (%0.1e)", data_.legEntry.c_str(), data_.h->Integral(0, data_.h->GetNbinsX() + 1));
        leg->AddEntry(data_.h.get(), label, data_.drawOptions.c_str());
        smartMax(data_.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, true);

        //background
        for(auto& entry : bgEntries_)
        {
            //set fill color so BG will have solid fill
            entry.setFillColor();

            sprintf(label, "%s (%0.1e)", entry.legEntry.c_str(), entry.h->Integral(0, entry.h->GetNbinsX() + 1));

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), label, "F");
        }
        smartMax(hbgSum, leg, static_cast<TPad*>(gPad), min, max, lmax, false);

        //signal 
        for(auto& entry : sigEntries_)
        {
            //get new histogram
            entry.histName = histName;
            entry.rebin = rebin;
            entry.retrieveHistogram();

            //add histograms to TLegend
            leg->AddEntry(entry.h.get(), entry.legEntry.c_str(), "L");
            smartMax(entry.h.get(), leg, static_cast<TPad*>(gPad), min, max, lmax, false);
        }

        //Set Canvas margin (gPad is root magic to access the current pad, in this case canvas pad1)
        //gPad->SetLeftMargin(0.12);
        //gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0.08);
        //gPad->SetBottomMargin(0.0);

        //create a dummy histogram to act as the axes
        histInfo dummy(new TH1D("dummy", "dummy", 1000, data_.h->GetBinLowEdge(1), data_.h->GetBinLowEdge(data_.h->GetNbinsX()) + data_.h->GetBinWidth(data_.h->GetNbinsX())));
        // set pad margins: setupPad(left, right, top, bottom)
        dummy.setupPad(0.12, 0.06, 0.08, 0.0);
        dummy.setupAxes(1.1, 1.0, 0.06, 0.06, 0.05, 0.05);
        dummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        //dummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        dummy.h->GetXaxis()->SetTickLength(0.03);
        dummy.h->GetYaxis()->SetTickLength(0.03);

        //Set the y-range of the histogram
        if(isLogY)
        {
            double locMin = std::min(0.2, std::max(0.2, 0.05 * min));
            double legSpan = (log10(3*max) - log10(locMin)) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            double legMin = legSpan + log10(locMin);
            if(log10(lmax) > legMin)
            {
                double scale = (log10(lmax) - log10(locMin)) / (legMin - log10(locMin));
                max = pow(max/locMin, scale)*locMin;
            }
            dummy.h->GetYaxis()->SetRangeUser(locMin, 10*max);
        }
        else
        {
            double locMin = 0.0;
            double legMin = (1.2*max - locMin) * (leg->GetY1() - gPad->GetBottomMargin()) / ((1 - gPad->GetTopMargin()) - gPad->GetBottomMargin());
            if(lmax > legMin) max *= (lmax - locMin)/(legMin - locMin);
            dummy.h->GetYaxis()->SetRangeUser(0.0, max*1.3);
        }
        //set x-axis range
        if(xmin < xmax) dummy.h->GetXaxis()->SetRangeUser(xmin, xmax);

        //draw dummy axes
        dummy.draw();

        //Switch to logY if desired
        gPad->SetLogy(isLogY);

        //plot background stack
        bgStack->Draw("same");

        //plot signal histograms
        for(const auto& entry : sigEntries_)
        {
            entry.draw();
        }

        //plot data histogram
        data_.draw();

        //plot legend
        leg->Draw("same");

        //Draw dummy hist again to get axes on top of histograms
        dummy.draw("AXIS");

        //Draw CMS and lumi lables
        char lumistamp[128];
        sprintf(lumistamp, "%.1f fb^{-1} (13 TeV)", lumi / 1000.0);

        TLatex mark;
        mark.SetNDC(true);

        //Draw CMS mark
        mark.SetTextAlign(11);
        mark.SetTextSize(0.050);
        mark.SetTextFont(61);
        mark.DrawLatex(gPad->GetLeftMargin(), 1 - (gPad->GetTopMargin() - 0.017), "CMS"); // #scale[0.8]{#it{Preliminary}}");
        mark.SetTextSize(0.040);
        mark.SetTextFont(52);
        mark.DrawLatex(gPad->GetLeftMargin() + 0.11, 1 - (gPad->GetTopMargin() - 0.017), "Preliminary");

        //Draw lumistamp
        mark.SetTextFont(42);
        mark.SetTextAlign(31);
        mark.DrawLatex(1 - gPad->GetRightMargin(), 1 - (gPad->GetTopMargin() - 0.017), lumistamp);

        // lower plot will be in pad2
        c->cd();          // Go back to the main canvas before defining pad2
        TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.0, 1, 0.3);
        //pad2->SetTopMargin(0);
        //pad2->SetBottomMargin(0.2);
        pad2->SetGridy(); // Horizontal grid
        pad2->Draw();
        pad2->cd();       // pad2 becomes the current pad        

        //histInfo dummy(new TH1D("dummy", "dummy", 1000, data_.h->GetBinLowEdge(1), data_.h->GetBinLowEdge(data_.h->GetNbinsX()) + data_.h->GetBinWidth(data_.h->GetNbinsX())));
        //TH1* ratio = (TH1*)data_.h->Clone();

        //Set Canvas margin (gPad is root magic to access the current pad, in this case pad2)
        //gPad->SetLeftMargin(0.12);
        //gPad->SetRightMargin(0.06);
        //gPad->SetTopMargin(0);
        //gPad->SetBottomMargin(0.40);

        //make ratio dummy
        histInfo ratioDummy(new TH1D("rdummy", "rdummy", 1000, data_.h->GetBinLowEdge(1), data_.h->GetBinLowEdge(data_.h->GetNbinsX()) + data_.h->GetBinWidth(data_.h->GetNbinsX())));
        ratioDummy.h->GetXaxis()->SetTitle(xAxisLabel.c_str());
        //ratioDummy.h->GetYaxis()->SetTitle(yAxisLabel.c_str());
        ratioDummy.h->GetYaxis()->SetTitle("Data / BG");
        ratioDummy.h->GetXaxis()->SetTickLength(0.1);
        ratioDummy.h->GetYaxis()->SetTickLength(0.045);
        ratioDummy.setupAxes(1.2, 0.4, 0.15, 0.15, 0.13, 0.13);
        ratioDummy.h->GetYaxis()->SetNdivisions(6, 5, 0);
        ratioDummy.h->GetXaxis()->SetRangeUser(xmin, xmax);
        ratioDummy.h->GetYaxis()->SetRangeUser(0.5, 1.5);
        ratioDummy.h->SetStats(0);
        //ratioDummy.h->SetMinimum(0.5);
        //ratioDummy.h->SetMaximum(1.5);

        //Make ratio histogram for data / background.
        histInfo ratio((TH1*)data_.h->Clone());

        // set pad margins: setupPad(left, right, top, bottom)
        ratio.setupPad(0.12, 0.06, 0.0, 0.40);
        
        ratio.drawOptions = "ep";
        ratio.color = kBlack;

        //ratio.h->SetLineColor(kBlack);
        //ratio.h->Sumw2();
        //ratio.h->SetStats(0);
        ratio.h->Divide(hbgSum);
        ratio.h->SetMarkerStyle(21);

        ratioDummy.draw();
        ratio.draw("same");

        //save new plot to file
        std::string name = histName;
        size_t pos = name.find("/");
        do
        {
            name[pos] = '_';
        }
        while( (pos = name.find("/", pos + 1)) != std::string::npos);
        c->Print((name + ".png").c_str());

        //clean up dynamic memory
        delete c;
        delete leg;
        delete bgStack;
        delete hbgSum;

//        printf("Finish plotting histogram %s\n", histName.c_str());
    }
};

int main()
{
    //entry for data
    //this uses the initializer syntax to initialize the histInfo object
    //               leg entry  root file                 draw options  draw color
    histInfo dataPhoton = {"Data",    "../TT_Data_SinglePhoton-2018-3-26_noWgt_v2.root", "PEX0",       kBlack};
    histInfo dataMuon   = {"Data",    "../TT_Data_SingleMuon-2018-3-26_noWgt_v2.root",   "PEX0",       kBlack};
    histInfo dataMET    = {"Data",    "../TT_Data_MET-2018-3-26_noWgt_v2.root",          "PEX0",       kBlack};
    histInfo dataJetHT  = {"Data",    "../TT_Data_JetHT-2018-3-26_noWgt_v2.root",        "PEX0",       kBlack};

    //vector summarizing background histograms to include in the plot
    std::vector<histInfo> bgEntries = {
        {"QCD",                "../TT_QCD-2018-3-26_noWgt_v2.root",               "hist", kOrange},
        //{"t#bar{t}",           "../TT_TTbarSingleLep-2018-3-26_noWgt_v2.root",    "hist", kRed},
        {"t#bar{t}",           "../TT_TTbar-2018-3-26_noWgt_v2.root",             "hist", kRed},
        {"G+Jets",             "../TT_GJets-2018-3-26_noWgt_v2.root",             "hist", kGreen + 2},
        {"Z#rightarrowll",     "../TT_DYJetsToLL-2018-3-26_noWgt_v2.root",        "hist", kBlue},
        {"Z#rightarrow#nu#nu", "../TT_ZJetsToNuNu-2018-3-26_noWgt_v2.root",       "hist", kBlue + 2},
        {"W+Jets",             "../TT_WJetsToLNu-2018-3-26_noWgt_v2.root",        "hist", kGray},
        {"TTG",                "../TT_TTG-2018-3-26_noWgt_v2.root",               "hist", kYellow + 3},
        {"TTZ",                "../TT_TTZ-2018-3-26_noWgt_v2.root",               "hist", kMagenta + 2},
        {"diboson",            "../TT_Diboson-2018-3-26_noWgt_v2.root",           "hist", kPink - 2}
    };

    //vector summarizing signal histograms to include in the plot
    std::vector<histInfo> sigEntries = {
//        {"T2tt (1000, 1)", "myhistos/Signal_fastsim_T2tt_mStop-1000.root", "hist", kGreen + 2},
    };

    //make plotter object with the required sources for histograms specified
    Plotter pltMET(dataMET, bgEntries, sigEntries);
    Plotter pltPhoton(dataPhoton, bgEntries, sigEntries);
    Plotter pltMuon(dataMuon, bgEntries, sigEntries);
    Plotter pltJetHT(dataJetHT, bgEntries, sigEntries);

    std::vector<std::pair<std::string, Plotter*>> controlRegions = {
        {"ttbar", &pltMET},
        {"ttbarNob", &pltMET},
        {"photon", &pltPhoton},
        {"dilepton", &pltMuon},
        {"ttbarLep", &pltMuon},
        {"QCD", &pltJetHT},
        {"QCDb", &pltJetHT},
    };

    for(auto& cr : controlRegions)
    {

        cr.second->plot(cr.first + "/HT", "H_{T} [GeV]", "Events", true, 0, 2000, 5);
        cr.second->plot(cr.first + "/MET", "MET [GeV]", "Events", true, 0, 1000, 5);
        cr.second->plot(cr.first + "/nJets", "N_{j}", "Events", true);
        cr.second->plot(cr.first + "/nBJets", "N_{b}", "Events", true, -0.5, 9.5);
        cr.second->plot(cr.first + "/nTops", "N_{t}", "Events", true);
        cr.second->plot(cr.first + "/fakerateHT2", "H_{T} [GeV]", "Events", true, 0, 2000, 5);
        cr.second->plot(cr.first + "/fakerateNj2", "N_{j}");
        cr.second->plot(cr.first + "/fakerateNb2", "N_{b}", "Events", true, -0.5, 9.5);
        cr.second->plot(cr.first + "/randomTopPt", "rand top p_{T} [GeV]", "Events", false, -1, -1, 5);
        cr.second->plot(cr.first + "/randomTopCandPt", "rand top p_{T} [GeV]", "Events", false, -1, -1, 5);
        cr.second->plot(cr.first + "/nVertices", "NPV");
    }
}
