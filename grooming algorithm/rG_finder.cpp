// A simple Pythia example, just to demonstrate we can produce jets and do useful things with them.

#include <Pythia8/Pythia.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/Selector.hh> 
#include <fastjet/contrib/SoftDrop.hh>

#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TH2D.h>

using namespace Pythia8;
using namespace fastjet;
using namespace std;

int main()
{
    // Default values for collision parameters
    // Recall all energies are in GeV
    int com_energy = 5020;
    int nevents = 10000; // Default value
    int minpthat = 170;
    int min_jet_pt = 50;
    pair<int, int> pT_ranges[4] = {{158, 200}, {200, 315}, {315,501}, {501,1000}}; //GeV ranges
    double jet_radius = 0.4;
    double max_eta = 2.1;

    //quark + gluon beta = 0.0
    //quark + gluon zcut = .2


    // Softdrop parameters - 
    // beta: {0-1} (.2 increments)
    // z_cut: {.1-.5} (.1 increments)
    double beta[6] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0}; // beta values
    double z_cut[5] = {0.1, 0.2, 0.3, 0.4, 0.5}; // z_cut values

    cout << "Chose number of events." << endl;
    cin >> nevents;

    Pythia pythia;
    
    // Event generator settings
    pythia.readString("PhaseSpace:pTHatMin = " + to_string(minpthat));
    pythia.readString("Beams:eCM = " + to_string(com_energy));
    pythia.readString("HardQCD:all = on");
    pythia.readString("Tune:pp = 21");
    
    // Initialize with these settings
    if (!pythia.init()) return 1;

    JetDefinition fjdef(antikt_algorithm, jet_radius); // antikt algorithm


    //softdrop vector in order of {0, 0.1}, {0, 0.2}, ..., {1, 0.5}
    contrib::SoftDrop* softdrops[6][5];
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 5; j++) {
            softdrops[i][j] = new contrib::SoftDrop(beta[i], z_cut[j], jet_radius);
        }
    }
    
    //contrib::SoftDrop softdrop(beta, z_cut, jet_radius);



    //defining Histograms
    int Nbins = 80;
    double xmin = log10(.003);
    double xmax = log10(.4);
    double logwidth = (xmax - xmin)/Nbins;
    vector<double> bins;
    for(int i = 0; i <= Nbins; i++){
        bins.push_back(pow(10, xmin + i * logwidth) );
    }

    TH1D* rG_Histogram[4][6][5];
    TH1D* quark_rG_Histogram[4];
    TH1D* gluon_rG_Histogram[4];
    TH1D* count_Histogram[4];

    TH2D *h_og;
    TH2D *h_z01; 
    TH2D *h_z03;



    for (int pT = 0; pT < 4; ++pT) {
        for (int b = 0; b < 6; ++b) {
            for (int z = 0; z < 5; ++z) {
                rG_Histogram[pT][b][z] = new TH1D(Form("rG_Histogram_%d_%d_%d", pT, b, z), 
                                                    Form("rG Histogram w/ pT between %d and %d. [zcut: %.1f]", pT_ranges[pT].first,pT_ranges[pT].second, z_cut[z]), 
                                                    Nbins, &bins[0]);
            }
        }
    }

    for (int pT = 0; pT < 4; ++pT) {
        //initialize pT histograms
        quark_rG_Histogram[pT] = new TH1D(Form("quark_rG_Histogram_%d", pT), "", Nbins, &bins[0]);
        gluon_rG_Histogram[pT] = new TH1D(Form("gluon_rG_Histogram_%d", pT), "", Nbins, &bins[0]);
        //initialize count histograms
        count_Histogram[pT] = new TH1D(Form("count_%d", pT), "", 2, 0, 2);
        count_Histogram[pT]->GetXaxis()->SetBinLabel(1, "Quark");
        count_Histogram[pT]->GetXaxis()->SetBinLabel(2, "Gluon");
        count_Histogram[pT]->SetMinimum(0);
    }

    bool create2DHistograms = false;
    
    for (int iEvt = 0; iEvt < nevents; iEvt++)
    {
        // Generate event
        if (!pythia.next()) continue;

        
        // Get particle info of final particles
        vector<PseudoJet> particles;
        vector<PseudoJet> quarks;
        vector<PseudoJet> gluons;
        for (auto i : pythia.event)
        {
            if (abs(i.status()) == 23 && abs(i.id()) == 21)
            {
                PseudoJet j(i.px(), i.py(), i.pz(), i.e());
                gluons.push_back(j);
            }
            if (abs(i.status()) == 23 && abs(i.id()) < 7)
            {
                if (i.idAbs() == 7) cout << "Top Quark Scattering" << endl;
                PseudoJet j(i.px(), i.py(), i.pz(), i.e());
                quarks.push_back(j);
                
            }
            if (i.isFinal() && abs(i.eta()) < max_eta)
            {
                PseudoJet j(i.px(), i.py(), i.pz(), i.e());
                j.set_user_index(i.id());
                particles.push_back(j);
                
            }
        }

        // Find jets using antikT algorithm
        ClusterSequence foundJets(particles, fjdef);
        vector<PseudoJet> truthJets = foundJets.inclusive_jets();

        // Apply pT cut
        for (auto it = truthJets.begin(); it != truthJets.end();)
        {
            if (it->pt() < min_jet_pt) {
                it = truthJets.erase(it);
            }
            else {
                it++;
            }
        }


        // Groom jets with SoftDrop and find rg
        Selector sel = SelectorCircle(.4);

        for (auto jet : truthJets){
            
            //fill three histograms
            if(!create2DHistograms) {
                PseudoJet groomed_01 = (*softdrops[0][0])(jet);
                auto groomed_info_01 = groomed_01.structure_of<contrib::SoftDrop>();
                double rG_01 = groomed_info_01.has_substructure () ? groomed_info_01.delta_R() : -1;

                PseudoJet groomed_03 = (*softdrops[0][2])(jet);
                auto groomed_info_03 = groomed_03.structure_of<contrib::SoftDrop>();
                double rG_03 = groomed_info_03.has_substructure () ? groomed_info_03.delta_R() : -1;

                if(rG_01 > 0 && rG_03 > 0 && rG_03 < (.5 * rG_01)) {
                    double minEta_og = 1e9;
                    double maxEta_og = -1e9;
                    double minEta_z01 = 1e9;
                    double maxEta_z01 = -1e9;
                    double minEta_z03 = 1e9;
                    double maxEta_z03 = -1e9;

                    double minPhi_og = 1e9;
                    double maxPhi_og = -1e9;
                    double minPhi_z01 = 1e9;
                    double maxPhi_z01 = -1e9;
                    double minPhi_z03 = 1e9;
                    double maxPhi_z03 = -1e9;

                    cout << "rG condition met!";
                    create2DHistograms = true;

                    for(auto &p : jet.constituents()) {
                        if(p.eta() < minEta_og) minEta_og = p.eta();
                        if(p.eta() > maxEta_og) maxEta_og = p.eta();

                        if(p.phi() < minPhi_og) minPhi_og = p.phi();
                        if(p.phi() > maxPhi_og) maxPhi_og = p.phi();

                        h_og = new TH2D("h_og", "Original Cluster;#eta;#phi", 25, minEta_og, maxEta_og, 25, minPhi_og, maxPhi_og);
                        h_og->GetXaxis()->SetRangeUser(minEta_og, maxEta_og);
                        h_og->GetYaxis()->SetRangeUser(minPhi_og, maxPhi_og);
                    }
                    for(auto& p : jet.constituents()) {
                        h_og->Fill(p.eta(), p.phi(), jet.pt());
                    }


                    for(auto &p : groomed_info_01.constituents(groomed_01)) {
                        if(p.eta() < minEta_z01) minEta_z01 = p.eta();
                        if(p.eta() > maxEta_z01) maxEta_z01 = p.eta();

                        if(p.phi() < minPhi_z01) minPhi_z01 = p.phi();
                        if(p.phi() > maxPhi_z01) maxPhi_z01 = p.phi();

                        h_z01 = new TH2D("h_z01", "Z-Cut: .1;#eta;#phi", 25, minEta_z01, maxEta_z01, 25, minPhi_z01, maxPhi_z01);
                        h_z01->GetXaxis()->SetRangeUser(minEta_z01, maxEta_z01);
                        h_z01->GetYaxis()->SetRangeUser(minPhi_z01, maxPhi_z01);
                    }
                    for(auto &p : groomed_info_01.constituents(groomed_01)) {
                        h_z01->Fill(p.eta(), p.phi(), p.pt());
                    }

                    for(auto &p : groomed_info_03.constituents(groomed_03)) {
                        if(p.eta() < minEta_z03) minEta_z03 = p.eta();
                        if(p.eta() > maxEta_z03) maxEta_z03 = p.eta();

                        if(p.phi() < minPhi_z03) minPhi_z03 = p.phi();
                        if(p.phi() > maxPhi_z03) maxPhi_z03 = p.phi();

                        h_z03 = new TH2D("h_z03", "Z-Cut: .3;#eta;#phi", 25, minEta_z03, maxEta_z03, 25, minPhi_z03, maxPhi_z03);
                        h_z03->GetXaxis()->SetRangeUser(minEta_z03, maxEta_z03);
                        h_z03->GetYaxis()->SetRangeUser(minPhi_z03, maxPhi_z03);
                    }
                    for(auto &p : groomed_info_03.constituents(groomed_03)) {
                        h_z03->Fill(p.eta(), p.phi(), p.pt());
                    }

                }
            }







            for(int i = 0; i < 6; i++) {        //beta values
                for(int j = 0; j < 5; j++) {    //zcut values
                    PseudoJet groomed_jet = (*softdrops[i][j])(jet);
                    auto groomed_info = groomed_jet.structure_of<contrib::SoftDrop>();
                    double rg = groomed_info.delta_R(); // delta R of the last iteration (definition of rg)
                    bool failed = !groomed_info.has_substructure(); // Never passed the soft drop condition
            
                    if(!failed){
                        if (jet.pt() >= pT_ranges[0].first && jet.pt() < pT_ranges[0].second) {
                                    rG_Histogram[0][i][j]->Fill(rg);
                            }
                        else if (jet.pt() >= pT_ranges[1].first && jet.pt() < pT_ranges[1].second) {
                                    rG_Histogram[1][i][j]->Fill(rg);
                            }
                        else if (jet.pt() >= pT_ranges[2].first && jet.pt() < pT_ranges[2].second) {
                                    rG_Histogram[2][i][j]->Fill(rg);
                                }     
                        else if (jet.pt() >= pT_ranges[3].first && jet.pt() < pT_ranges[3].second) {
                                    rG_Histogram[3][i][j]->Fill(rg);
                                }
                        

                        if(i == 0 && j == 1 ) {
                            // If the jet is close enough to the original quark or gluon, call it a quark or gluon jet.
                            sel.set_reference(jet);
                            vector<PseudoJet> nearby_quarks = sel(quarks);
                            vector<PseudoJet> nearby_gluons = sel(gluons);
                            

                            // Print their size
                            //cout << "nearby_quarks.size() = " << nearby_quarks.size() << endl;
                            //cout << "nearby_gluons.size() = " << nearby_gluons.size() << endl;
                        
                            if (nearby_quarks.size() == 0 && nearby_gluons.size() == 0) {
                                //cout << "Could not match parton to jet!" << endl;
                            }
                            else if (nearby_quarks.size() > 0 && nearby_gluons.size() == 0) {
                                // Fill quark jets histogram
                                if (jet.pt() >= pT_ranges[0].first && jet.pt() < pT_ranges[0].second) {
                                        quark_rG_Histogram[0]->Fill(rg);
                                        count_Histogram[0]->Fill(.5);
                                        
                                }
                                else if (jet.pt() >= pT_ranges[1].first && jet.pt() < pT_ranges[1].second) {
                                            quark_rG_Histogram[1]->Fill(rg);
                                            count_Histogram[1]->Fill(.5);
                                           
                                    }
                                else if (jet.pt() >= pT_ranges[2].first && jet.pt() < pT_ranges[2].second) {
                                            quark_rG_Histogram[2]->Fill(rg);
                                            count_Histogram[2]->Fill(.5);
                                         
                                        }     
                                else if (jet.pt() >= pT_ranges[3].first && jet.pt() < pT_ranges[3].second) {
                                            quark_rG_Histogram[3]->Fill(rg);
                                            count_Histogram[3]->Fill(.5);
                                            
                                        }
                            }
                            else if (nearby_quarks.size() == 0 && nearby_gluons.size() > 0) {
                                // Fill gluon jets histogram
                                if (jet.pt() >= pT_ranges[0].first && jet.pt() < pT_ranges[0].second) {
                                        gluon_rG_Histogram[0]->Fill(rg);
                                        count_Histogram[0]->Fill(1.5);
                              
                                }
                                else if (jet.pt() >= pT_ranges[1].first && jet.pt() < pT_ranges[1].second) {
                                        gluon_rG_Histogram[1]->Fill(rg);
                                        count_Histogram[1]->Fill(1.5);
                                        
                                    }
                                else if (jet.pt() >= pT_ranges[2].first && jet.pt() < pT_ranges[2].second) {
                                        gluon_rG_Histogram[2]->Fill(rg);
                                        count_Histogram[2]->Fill(1.5);
                                    
                                    }     
                                else if (jet.pt() >= pT_ranges[3].first && jet.pt() < pT_ranges[3].second) {
                                        gluon_rG_Histogram[3]->Fill(rg);
                                        count_Histogram[3]->Fill(1.5);
                                       
                                    }

                            }
                            else {
                                //cout << "Multiple partons matched! Can not unambiguously assign quark or gluon to jet." << endl;
                            }

                            
                        }


                    }



                }
            }
        }


    }


    TFile *f = TFile::Open("jet_girths.root", "RECREATE");


    // Normalize and write all histograms in the 3D array
    for (int pT = 0; pT < 4; ++pT) {
        TH1D* quark_hist = quark_rG_Histogram[pT];
        int totalQuarkJets = quark_hist->GetEntries();
        for (int b = 1; b <= quark_hist->GetNbinsX(); ++b) {
            double bin_content = quark_hist->GetBinContent(b);
            double bin_width = quark_hist->GetBinWidth(b);
            double normalizedBinValue = (totalQuarkJets > 0) ? bin_content / (bin_width * totalQuarkJets) : 0;
            quark_hist->SetBinContent(b, normalizedBinValue);
        }
        quark_hist->Write();


        TH1D* gluon_hist = gluon_rG_Histogram[pT];
        int totalGluonJets = gluon_hist->GetEntries();
        for (int b = 1; b <= gluon_hist->GetNbinsX(); ++b) {
            double bin_content = gluon_hist->GetBinContent(b);
            double bin_width = gluon_hist->GetBinWidth(b);
            double normalizedBinValue = (totalGluonJets > 0) ? bin_content / (bin_width * totalGluonJets) : 0;
            gluon_hist->SetBinContent(b, normalizedBinValue);
        }
        gluon_hist->Write();


        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 5; ++j) {
                TH1D* hist = rG_Histogram[pT][i][j];
                int totalJets = hist->GetEntries();
                for (int b = 1; b <= hist->GetNbinsX(); ++b) {
                    double bin_content = hist->GetBinContent(b);
                    double bin_width = hist->GetBinWidth(b);
                    double normalizedBinValue = (totalJets > 0) ? bin_content / (bin_width * totalJets) : 0;
                    hist->SetBinContent(b, normalizedBinValue);
                }
                hist->Write();
            }
        }

        count_Histogram[pT]->Write();

    }


    h_og->Write();
    h_z01->Write();
    h_z03->Write();



    f->Close();

    pythia.stat();
    return 0;
}