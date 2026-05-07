//ROOT
#include <TFile.h> 
#include <TGraph.h> 
#include <TTree.h> 
#include <TCanvas.h>
#include <TParameter.h> 
#include <TError.h> 
#include <TString.h> 
#include <TH2D.h> 
#include <TAxis.h> 
#include <TLegend.h> 

//stdlib
#include <string>
#include <vector> 

#include "compute_heat_capacity.C"

int draw_heat_capacities()
{
    std::string input_cooling = "data/out_cooling_n-%i_beta-19.83-66.67_12623676.root"; 

    std::string input_heating = "data/out_heating_n-%i_beta-12.24-24.00_12623703.root"; 

    int n_min=4; 
    int n_max=7; 

    std::vector<TGraph*> graphs_C, graphs_E; 

    double g2_min{6./66.7}, g2_max{6./12.24}; 

    
    double frame_x[] = { g2_min, g2_max }; 
    double frame_y[] = { 0., 0. }; 

    auto canvC = new TCanvas; 
    auto legend = new TLegend; 
    canvC->cd(); 
    auto graph_frameC = new TGraph(2, frame_x, frame_y); 
    graph_frameC->SetMinimum(0.);
    graph_frameC->SetMaximum(2.5e-3);  
    graph_frameC->SetTitle("Heat Capacity;g_{0}^{2};dE/dg_{0}^{2}"); 
    graph_frameC->Draw(); 

    for (int n=n_min; n<=n_max; n++) {
        std::string path_cooling = Form(input_cooling.c_str(), n); 
        std::string path_heating = Form(input_heating.c_str(), n); 
        
        TGraph* gcC, *ghC; 
        TGraph* gcE, *ghE; 

        std::vector<double> pts_cg2, pts_hg2, pts_cC, pts_hC, pts_E; 

        compute_heat_capacity(path_cooling, "", pts_cg2, pts_cC, pts_E);         
        compute_heat_capacity(path_heating, "", pts_hg2, pts_hC, pts_E); 

        gcE = new TGraph(pts_cg2.size(), pts_cg2.data(), pts_E.data()); 
        gcE->SetTitle("Energy;g_{0}^{2};dE/dg_{0}^{2}"); 
        
        gcC = new TGraph(pts_cg2.size(), pts_cg2.data(), pts_cC.data()); 
        gcC->SetLineColor(n-n_min+1); 
        gcC->SetLineStyle(kDashed); 
        gcC->Draw("SAME"); 
        legend->AddEntry(gcC, Form("%i-cool",n)); 
        
        ghC = new TGraph(pts_hg2.size(), pts_hg2.data(), pts_hC.data()); 
        ghC->SetLineColor(n-n_min+1); 
        ghC->SetLineStyle(kSolid); 
        ghC->Draw("SAME"); 
        legend->AddEntry(ghC, Form("%i-heat",n)); 
        
        graphs_C.push_back(gcC); 
        graphs_E.push_back(gcE); 
    }

    legend->Draw(); 
    canvC->SaveAs("Heat_capacity.png"); 

    return 0; 
}