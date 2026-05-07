#ifndef compute_heat_capacity_C 
#define compute_heat_capacity_C 

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

int compute_heat_capacity(std::string path_infile, std::string path_outfile, 
    std::vector<double>& pts_g2,
    std::vector<double>& pts_C,
    std::vector<double>& pts_E)
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

    const double g0_min = 6./beta_max; //std::sqrt(6./beta_max);
    const double g0_max = std::min( 6./beta_min, 1. ); //std::min( std::sqrt(6./beta_min), 1. );
    
    TH2D* hist = new TH2D(
        "hist", "data",
        n_bins, g0_min, g0_max, 
        energy_res,0.,energy_max
    );
    
    tree->Draw("energy:6./beta>>hist", "", "goff");  
    
    auto x_ax = hist->GetXaxis(); 
    auto y_ax = hist->GetYaxis(); 

    pts_g2.clear(); pts_g2.reserve(n_bins); 
    pts_C.clear(); pts_C.reserve(n_bins); 
    pts_E.clear(); pts_E.reserve(n_bins); 

    for (int bx=1; bx<=x_ax->GetNbins(); bx++) {

        double g2 = x_ax->GetBinCenter(bx); 

        double sum_EE =0.;
        double sum_E  =0.; 
        double sum_N  =0.; 
        
        pts_g2.push_back( g2 );

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

        pts_C.push_back( ( sum_EE - sum_E*sum_E )/(g2*g2) ); 
        pts_E.push_back( sum_E ); 
    }

    file->Close(); 
    return 0; 
}

#endif 

