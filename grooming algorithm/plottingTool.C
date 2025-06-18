#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TLegend.h>
#include <TF1.h>



void plottingTool(std::string fileName){
    double beta_values[6] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    TFile *input = new TFile(fileName.c_str());
    TCanvas *c1 = new TCanvas();
    int meow[5] = {158, 200, 315, 501, 1000};
    

    
    //quark vs gluon rG histograms
    for (int pT = 0; pT < 4; pT ++) {
        double meanRatio = 0;
        auto legend2 = new TLegend(.86,.8,.98,.94, "Quarks V Gluons", "brNDC");
        legend2->SetHeader("Quarks V Gluons","C");
        legend2->SetTextSize(.023);
        auto legend3 = new TLegend(.85,.70,.985,.78);
        

        double maxY2 = 0;
        TH1D *maximumChecker = (TH1D*)input->Get(Form("quark_rG_Histogram_%d", pT));

        if (maximumChecker->GetMaximum() > maxY2) {
            maxY2 = maximumChecker->GetMaximum();
        }
        maximumChecker = (TH1D*)input->Get(Form("gluon_rG_Histogram_%d", pT));
        if (maximumChecker->GetMaximum() > maxY2) {
            maxY2 = maximumChecker->GetMaximum();
        }
            
        //create gluon histogram distribution
        TH1D *special_h = (TH1D*)input->Get(Form("quark_rG_Histogram_%d", pT));
            
        special_h->SetMaximum(maxY2 * 1.1);
        special_h->GetXaxis()->SetTitle("rG");
        special_h->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{dN_{jet}}{dr_{g}}");
        special_h->GetXaxis()->SetTitleColor(4);
        special_h->GetYaxis()->SetTitleColor(4);
        special_h->GetXaxis()->SetTitleSize(0.045);
        special_h->GetYaxis()->SetTitleSize(.035);
        special_h->SetTitle(Form("Quark vs. Gluon rG Histogram w/ pT between: %d - %d", meow[pT], meow[pT+1]));
        special_h->SetLineColor(2);
        special_h->SetStats(false);
        legend2->AddEntry(special_h, "Quarks", "l");
        meanRatio = special_h->GetMean();

        special_h->Draw("HIST");
        

        //create gluon histogram distribution
        special_h = (TH1D*)input->Get(Form("gluon_rG_Histogram_%d", pT));
        special_h->SetLineWidth(2);
        special_h->SetLineColor(3);
        legend2->AddEntry(special_h, "Gluons", "l");
        meanRatio = meanRatio / (special_h->GetMean());
       
        special_h->Draw("HIST SAME");
        

        //draw legend 3
        legend3->SetTextSize(.02);
        legend3->SetFillColor(30);
        legend3->SetHeader(Form("q-jet/g-jet mean:\n%.3f",meanRatio),"C");
        legend3->Draw();


        //draw plot
        legend2->Draw();
        c1->GetPad(0)->SetLogx();
        c1->Print(Form("doubleRefinedPlots/quark_v_gluon_pT_range-%d.pdf", pT));


    }

    c1->Clear();


    //make quark vs gluon count charts
    for (int pT = 0; pT < 4; pT++) {
        TH1D *countH = (TH1D*)input->Get(Form("count_%d", pT));
        countH->SetLineColor(4);
        countH->SetStats(0);
        countH->SetLineColor(2);
        countH->SetFillColor(2);
        countH->SetTitle(Form("Quark VS. Gluon Counts (pT range: %d - %d)", meow[pT], meow[pT+1]));
        countH->Draw();
        c1->GetPad(0)->SetLogx(0);
        c1->Print(Form("doubleRefinedPlots/quark_v_gluon_counts-%d.pdf", pT));

    }

    c1->Clear();

    //make rG plots
    for (int pT = 0; pT < 4; pT++) {
        TH1D *qH = (TH1D*)input->Get(Form("quark_rG_Histogram_%d", pT));

        for (int z = 0; z < 5; z++) {
            bool first = true;
            bool styled = false;
            int color = 2;

            auto legend = new TLegend(.86,.74,.98,.94, "Beta Values", "brNDC");
            legend->SetHeader("Beta Values", "C");
            legend->SetTextSize(.025);

            double maxY = 0;
            for(int i = 0; i < 6; i++) {
                TH1D *h = (TH1D*)input->Get(Form("rG_Histogram_%d_%d_%d", pT, i, z));
                if (h->GetMaximum() > maxY) {
                    maxY = h->GetMaximum();
                }
            }

            for (int b = 0; b < 6; b++) {
                TH1D *h = (TH1D*)input->Get(Form("rG_Histogram_%d_%d_%d", pT, b, z));
                h->SetMaximum(maxY * 1.1);

                if(!styled) {
                    h->GetXaxis()->SetTitle("rG");
                    h->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{dN_{jet}}{dr_{g}}");
                    h->GetXaxis()->SetTitleColor(4);
                    h->GetYaxis()->SetTitleColor(4);
                    h->GetXaxis()->SetTitleSize(0.045);
                    h->GetYaxis()->SetTitleSize(.035);
                    h->SetLineWidth(2);
                    styled = true;
                }
                if(first){
                    h->SetLineColor(1);
                    h->SetStats(false);
                    h->Draw("HIST");
                    legend->AddEntry(h, Form("%.1f", beta_values[b]), "l");
                    first = false;
                }else{
                    h->SetLineWidth(2);
                    h->SetLineColor(color);
                    legend->AddEntry(h, Form("%.1f", beta_values[b]), "l");
                    color++;
                    h->Draw("HIST SAME");
                }
                
                
            }

            legend->Draw();
            c1->GetPad(0)->SetLogx();
            c1->Print(Form("doubleRefinedPlots/jet_girths_%d_%d.pdf", pT, z));
        }
    }

    c1->Clear();

    //last but not least, create the three jet vs jet 2D histograms of constituents

    //(1)
    TH2D *jetConstituents_OG = (TH2D*)input->Get("h_og");
    //jetConstituents_OG->GetZaxis()->SetTitle("pT");
    jetConstituents_OG->SetStats(false);
    jetConstituents_OG->Draw("colz");
    c1->GetPad(0)->SetLogx(0);
    c1->Print("doubleRefinedPlots/xOGjet.pdf");

    c1->Clear();

    //(2)
    TH2D *jetConstituents_z01 = (TH2D*)input->Get("h_z01");
    //jetConstituents_z01->GetZaxis()->SetTitle("pT");
    jetConstituents_z01->SetStats(false);
    jetConstituents_z01->Draw("colz");
    c1->Print("doubleRefinedPlots/xZ01jet.pdf");

    c1->Clear();

    //(3)
    TH2D *jetConstituents_z03 = (TH2D*)input->Get("h_z03");
    //jetConstituents_z03->GetZaxis()->SetTitle("pT");
    jetConstituents_z03->SetStats(false);
    jetConstituents_z03->Draw("colz");
    c1->Print("doubleRefinedPlots/xZ03jet.pdf");

    input->Close();


}