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
//stdlib
#include <string>
#include <cmath>
#include <type_traits> 
#include <memory> 
#include <stdexcept> 
#include <vector> 

template<typename T> T GetParameter(const std::string& name, TFile* file)
{   
    //check to make sure this type is default-constructable
    static_assert(std::is_default_constructible<T>{}, "type 'T' must be default constructable");

    if (!file || !file->IsOpen()) {
        Error(__func__, "TFile is not open / valid"); 
        return T{};
    }

    TParameter<T>* param = file->Get<TParameter<T>>(name.c_str());
    if (!param) {
        Error(__func__, "TParameter could not be fetched"); 
        return T{}; 
    }
    return param->GetVal();
}

int compute_heat_capacity(std::string path_infile, std::string path_outfile)
{
    TFile* file; 
    file = new TFile(path_infile.c_str(), "READ");  
    
    /*try {
        
    } catch (const std::logic_error& e) {
        Error(__func__, "Something went wrong trying to define TFile at path '%s'\nwhat(): %s", path_infile.c_str(), e.what()); 
        return -1; 
    }*/

    std::cout << "file: " << file << std::endl; 

    if (!file || !file->IsOpen()) {
        Error(__func__, "TFile with following path could not be opened: '%s'", path_infile.c_str()); 
        return -1; 
    }
    
    //get the parameters we need
    const int steps_per_bin = 5; 

    const int n_steps = GetParameter<int>("n_steps", file); 

    const double beta_min = GetParameter<double>("beta_min", file); 
    const double beta_max = GetParameter<double>("beta_max", file); 

    const int n_bins = 100; 

    auto tree = file->Get<TTree>("data_tree"); 

    const double energy_max = 1.;
    const int energy_res = 500; 

    const double g0_min = std::sqrt(6./beta_max);
    const double g0_max = std::min( std::sqrt(6./beta_min), 1. );
    

    TH2D* hist = new TH2D(
        "hist", "data",
        n_bins, g0_min, g0_max, 
        energy_res,0.,energy_max
    );
    
    tree->Draw("energy:sqrt(6./beta)>>hist", "", "goff");  
    
    auto x_ax = hist->GetXaxis(); 
    auto y_ax = hist->GetYaxis(); 

    std::vector<double> pts_beta, pts_C; 

    pts_beta.reserve(n_bins); 
    pts_C   .reserve(n_bins); 

    for (int bx=1; bx<=x_ax->GetNbins(); bx++) {

        double g = x_ax->GetBinCenter(bx); 

        double sum_EE =0.;
        double sum_E =0.; 
        double sum_N =0.; 
        
        pts_beta.push_back( g );

        for (int by=1; by<y_ax->GetNbins(); by++) {

            //the stats of the bin
            double N = hist->GetBinContent(bx,by);

            //the 'energy' of this bin
            double E = y_ax->GetBinCenter(by);   

            sum_E   += E*N; 
            sum_EE  += E*E*N;
            sum_N   += N; 
        }  
        sum_E   *= 1./sum_N; 
        sum_EE  *= 1./sum_N; 

        pts_C.push_back( ( sum_EE - sum_E*sum_E )/(g*g) ); 
    }

    auto graph = new TGraph(n_bins, pts_beta.data(), pts_C.data()); 
    auto c = new TCanvas; 
    graph->SetTitle("Heat capacity;#g_{0};dE / dg_{0}"); 
    graph->Draw(); 
    c->SaveAs(path_outfile.c_str());
    
    return 0; 
}

