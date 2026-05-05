#include "TreeWriter.hxx"
#include <stdexcept> 
#include <iostream> 
#include <TString.h> 

//_____________________________________________________________________________________________
TreeWriter::TreeWriter(std::string tree_name, std::string path_outfile)
{
    try {
        fFile = std::unique_ptr<TFile>(new TFile(path_outfile.data(), "RECREATE"));
    } catch (const std::exception& e) {
        throw std::logic_error(Form("<TreeWriter::TreeWriter> something went wrong trying to open file '%s'.\nwhat(): %s",path_outfile.c_str(),e.what())); 
        return; 
    }

    //now, try to define the output branches
    fTree = new TTree(tree_name.c_str(), "data");
    
    fTree->Branch("beta",               &fBeta,             "beta/D");
    fTree->Branch("energy",             &fEnergy,           "energy/D");
    fTree->Branch("p_accept",           &fP_accept,         "p_accept/D");
    fTree->Branch("theta",              &fTheta,            "theta/D");
    fTree->Branch("frob_distance",      &fFrob_distance,    "frob_distance/D");
    fTree->Branch("time_per_update",    &fTime,             "time_per_update/D");
}
//_____________________________________________________________________________________________
void TreeWriter::CloseFile()
{
    if (fFile != nullptr && fFile->IsOpen()) {
        if (fTree) fTree->Write(); 
        fFile->Close(); 
        fTree = nullptr; 
    }
}
//_____________________________________________________________________________________________
void TreeWriter::WriteLine(double beta, double energy, double p_accept, double theta, double frob_distsance, double time)
{    
    if (!fTree) {
        throw std::logic_error("<TreeWriter::WriteLine> tree is null"); 
        return;
    }
    fWriteMutex.lock();

    fBeta=beta;
    fEnergy=energy;
    fP_accept=p_accept;
    fTheta=theta;
    fFrob_distance=frob_distsance;
    fTime=time;

    //write this line in the tree
    fTree->Fill(); 
    
    fWriteMutex.unlock(); 
} 
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
//_____________________________________________________________________________________________
