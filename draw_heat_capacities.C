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
#include <TF1.h> 
#include <TFitResult.h>
#include <TFitResultPtr.h> 

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

    void fit_gaus_to_vectors( 
        const std::vector<double>& pts_x,
        const std::vector<double>& pts_y,   
        double radius, 
        double &center, 
        double &sigma, 
        bool do_draw=false
    ) 
    {
        //first guess for the center of the fit. we assume here that the gaus we wanna fit is the tallest peak in the hist. 
        auto g = new TGraph(pts_x.size(), pts_x.data(), pts_y.data()); 

        double y_max=pts_y.front();
        double y_max_x=pts_x.front(); 
        for (int i=0; i<pts_x.size(); i++) {
            if (y_max < pts_y[i]) {
                y_max_x = pts_x[i];
                y_max = pts_y[i]; 
            }
        }

        auto x_ax = g->GetXaxis(); 
        center = y_max_x; 
        
        //first guess for sigma 
        sigma = radius/8.; 

        double x_min = std::max( center-radius, x_ax->GetXmin() );    
        double x_max = std::min( center+radius, x_ax->GetXmax() );

        //the guess for the 'base' will be based on the bins on either extreme end of the fit
        double base = 0.5*(
            g->Eval(center + radius) + g->Eval(center - radius)
        ); 

        //amplitude over background
        double amplitude = y_max - base; 
        
        //our function to fit the histogram with   
        auto gaus_fit = new TF1( "gaus_fit", [](const double* X, const double *par){
            double arg = (X[0]-par[1])/par[2]; 
            return par[0]*std::exp( -arg*arg/2. ) + par[3]; 
        }, x_min, x_max,  4); 
    
        //set the parameters 
        gaus_fit->SetParameter( 0, amplitude ); 
        gaus_fit->SetParameter( 1, center ); 
        gaus_fit->SetParameter( 2, sigma );  
        gaus_fit->SetParameter( 3, base ); 

        std::printf("amplitude: %f center: %f sigma %f base %f\n", amplitude,center,sigma,base); 
        
        auto fit_ptr = g->Fit("gaus_fit", (do_draw ? "S L R Q" : "S L N R Q")); 
        
        if (!fit_ptr.Get() || fit_ptr->IsValid()==false) {
            //throw std::logic_error("<fit_gaus_to_hist>: fit failed."); 
            Error(__func__, "Fit failed"); 
            return; 
        }

        center = fit_ptr->Parameter(1); 
        sigma  = std::fabs(fit_ptr->Parameter(2)); 
        
        return; 
    }
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

    //check the finite temperature scaling 
    std::vector<double> pts_N, pts_tc_cool, pts_tc_hot; 

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

            double center, sigma; 
            
            gcC = new TGraph(pts_cg2.size(), pts_cg2.data(), pts_cC.data()); 
            //gcC->SetLineColor(n-n_min+1); 
            gcC->SetLineStyle(kDashed); 
            gcC->SetLineWidth(3); 
            gcC->SetLineColor(kBlue);
            gcC->Draw("SAME"); 
            fit_gaus_to_vectors(pts_cg2, pts_cC, 0.025, center, sigma, true); 
            pts_tc_cool.push_back(center); 

            ghC = new TGraph(pts_hg2.size(), pts_hg2.data(), pts_hC.data()); 
            //ghC->SetLineColor(n-n_min+1); 
            ghC->SetLineStyle(kSolid); 
            ghC->SetLineWidth(3); 
            ghC->SetLineColor(kRed);
            ghC->Draw("SAME"); 
            fit_gaus_to_vectors(pts_hg2, pts_hC, 0.025, center, sigma, true); 
            pts_tc_hot.push_back(center); 

            pts_N.push_back( 4.*std::pow(n,4) );
            
            
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

    new TCanvas;
    auto legend = new TLegend; 
    
    auto g_cool = new TGraph(pts_N.size(), pts_N.data(), pts_tc_cool.data());
    g_cool->SetMaximum(0.5);
    g_cool->SetMinimum(0.1);
    g_cool->SetTitle("Phase transition temperature;N = 4*L^{4};T_{c}"); 
    g_cool->SetLineStyle(kDashed); 
    g_cool->SetLineWidth(3); 
    g_cool->SetLineColor(kBlue);
    g_cool->SetMarkerStyle(kOpenCircle);
    g_cool->SetMarkerColor(kBlue);  
    g_cool->SetMarkerSize(1);
    g_cool->Draw(); 

    auto scaling_fit = new TF1("fit_temp", [](const double *x, const double* par){
        return par[0]*std::pow( x[0]/par[1], -4. ) + par[2];
    }, pts_N.front(), pts_N.back(), 3);

    scaling_fit->SetParameter( 0, 0. );
    scaling_fit->SetParameter( 1, 1.e3 );
    scaling_fit->SetParameter( 2, pts_tc_cool.back() );

    g_cool->Fit("fit_temp", "R N Q"); 
    scaling_fit->SetLineColor(kBlue); 
    //scaling_fit->DrawCopy("SAME"); 

    auto g_hot = new TGraph(pts_N.size(), pts_N.data(), pts_tc_hot.data());
    g_hot->SetMaximum(0.5);
    g_hot->SetMinimum(0.1);
    g_hot->SetTitle("Phase transition temperature;N = 4*L^{4};T_{c}"); 
    g_hot->SetLineStyle(kSolid); 
    g_hot->SetLineWidth(3); 
    g_hot->SetLineColor(kRed);
    g_hot->SetMarkerStyle(kOpenCircle);
    g_hot->SetMarkerColor(kRed);  
    g_hot->SetMarkerSize(1);
    g_hot->Draw("SAME, PL"); 

    scaling_fit->SetParameter( 2, pts_tc_hot.back() );

    g_hot->Fit("fit_temp", "R N Q"); 
    scaling_fit->SetLineColor(kRed); 
    //scaling_fit->DrawCopy("SAME"); 

    return 0; 
}