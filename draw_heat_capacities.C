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
#include <cmath> 

#include "compute_heat_capacity.C"

namespace {
    double vector_min(const std::vector<double>& v) {
        double min=v.front(); 
        for (double x : v) { min = std::min(x,min); }
        return min; 
    };
    
    double vector_max(const std::vector<double>& v) {
        double max=v.front(); 
        for (double x : v) { max = std::max(x,max); }
        return max; 
    };
}

int draw_heat_capacities()
{
    std::string input_cooling = "data/out_cooling_n-%i_beta-19.83-66.67_12623676.root"; 

    std::string input_heating = "data/out_heating_n-%i_beta-12.24-24.00_12623703.root"; 

    int n_min=4; 
    int n_max=7; 

    double g2_min{6./66.7}, g2_max{6./12.24}; 
    
    double frame_x[] = { g2_min, g2_max }; 
    double frame_y[] = { 0., 0. }; 

    auto canvC = new TCanvas; 
    auto canvE = new TCanvas; 

    canvC->Divide(1, n_max-n_min+1, 0.,0.01);
    canvE->Divide(1, n_max-n_min+1, 0.,0.01); 

    for (int n=n_min; n<=n_max; n++) {

        
        std::string path_cooling = Form(input_cooling.c_str(), n); 
        std::string path_heating = Form(input_heating.c_str(), n); 
        
        TGraph* gcC, *ghC; 
        TGraph* gcE, *ghE; 

        std::vector<double> pts_cg2, pts_hg2, pts_cC, pts_hC, pts_cE, pts_hE; 

        compute_heat_capacity(path_cooling, "", pts_cg2, pts_cC, pts_cE);         
        compute_heat_capacity(path_heating, "", pts_hg2, pts_hC, pts_hE); 
    
        //heat capacity ------------------------------------------------------------------------------
        {
            auto graph_frameC = new TGraph(2, frame_x, frame_y); 
            graph_frameC->SetMinimum( 0. );
            graph_frameC->SetMaximum( 1.1*std::max(vector_max(pts_cC), vector_max(pts_hC)) ); 
            graph_frameC->SetTitle("Heat Capacity;g_{0}^{2};dE/dg_{0}^{2}"); 

            graph_frameC->GetXaxis()->SetNdivisions(25);
            graph_frameC->GetYaxis()->SetNdivisions(0);
        
            canvC->cd(n-n_min+1); 
            graph_frameC->Draw(); 

            gcC = new TGraph(pts_cg2.size(), pts_cg2.data(), pts_cC.data()); 
            //gcC->SetLineColor(n-n_min+1); 
            gcC->SetLineStyle(kDashed); 
            gcC->SetLineWidth(3); 
            gcC->SetLineColor(kBlue);
            gcC->Draw("SAME"); 
            
            ghC = new TGraph(pts_hg2.size(), pts_hg2.data(), pts_hC.data()); 
            //ghC->SetLineColor(n-n_min+1); 
            ghC->SetLineStyle(kSolid); 
            ghC->SetLineWidth(3); 
            ghC->SetLineColor(kRed);
            ghC->Draw("SAME"); 
            
            if (n==n_min) {
                auto legend = new TLegend; 
                legend->AddEntry(gcC, "cooling"); 
                legend->AddEntry(ghC, "heating");
                legend->Draw();  
            }
        }

        //energy ------------------------------------------------------------------------------
        {
            auto graph_frameC = new TGraph(2, frame_x, frame_y); 
            graph_frameC->SetMinimum( 0. );
            graph_frameC->SetMaximum( 1.1*std::max(vector_max(pts_cE), vector_max(pts_hE)) ); 
            graph_frameC->SetTitle("Energy;g_{0}^{2};Energy"); 

            graph_frameC->GetXaxis()->SetNdivisions(25);
            graph_frameC->GetYaxis()->SetNdivisions(10);
        
            canvE->cd(n-n_min+1); 
            graph_frameC->Draw(); 

            gcC = new TGraph(pts_cg2.size(), pts_cg2.data(), pts_cE.data()); 
            //gcC->SetLineColor(n-n_min+1); 
            gcC->SetLineStyle(kDashed); 
            gcC->SetLineWidth(3); 
            gcC->SetLineColor(kBlue);
            gcC->Draw("SAME"); 
            
            ghC = new TGraph(pts_hg2.size(), pts_hg2.data(), pts_hE.data()); 
            //ghC->SetLineColor(n-n_min+1); 
            ghC->SetLineStyle(kSolid); 
            ghC->SetLineWidth(3); 
            ghC->SetLineColor(kRed);
            ghC->Draw("SAME"); 
            
            if (n==n_min) {
                auto legend = new TLegend; 
                legend->AddEntry(gcC, "cooling"); 
                legend->AddEntry(ghC, "heating");
                legend->Draw();  
            }
        }
    }

    //legend->Draw(); 
    canvC->SaveAs("Heat_capacity.png"); 

    return 0; 
}