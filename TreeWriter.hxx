#ifndef TreeWriter_HXX
#define TreeWriter_HXX

#include <TFile.h>
#include <TTree.h>  
#include <vector> 
#include <string> 
#include <fstream> 
#include <memory> 
#include <mutex>

/// writes line of data to TFile  
class TreeWriter {
public: 
    TreeWriter(std::string tree_name, std::string path_outfile); 

    //don't allow copies of this object, so we can avoid accidentally writing to the same TFile 
    TreeWriter(const TreeWriter&) = delete; 

    void WriteLine(double beta, double energy, double p_accept, double theta, double frob_distsance, double time); 

    void CloseFile(); 

private: 

    //branches
    double fBeta, fEnergy, fP_accept, fTheta, fFrob_distance, fTime; 

    std::mutex fWriteMutex; 
    std::unique_ptr<TFile> fFile; 
    TTree* fTree; 
};

#endif