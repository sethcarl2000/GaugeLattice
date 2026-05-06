#include "TreeWriter.hxx"
#include <stdexcept> 
#include <iostream> 
#include <TString.h> 

//_____________________________________________________________________________________________
TreeWriter::TreeWriter(std::string tree_name, std::string path_outfile, long int n_entries)
{
    try {
        fFile = std::unique_ptr<TFile>(new TFile(path_outfile.data(), "RECREATE"));
    } catch (const std::exception& e) {
        throw std::logic_error(Form("<TreeWriter::TreeWriter> something went wrong trying to open file '%s'.\nwhat(): %s",path_outfile.c_str(),e.what())); 
        return; 
    }

    if (n_entries > 0) fData.reserve(n_entries); 
}
//_____________________________________________________________________________________________
void TreeWriter::CloseFile()
{   
    //write all elements to vector
    if (fFile == nullptr || !fFile->IsOpen()) {
        throw std::logic_error("<TreeWriter::CloseFile>: File is not open!"); 
        return; 
    }

    fTree = new TTree("data_tree", "data tree"); 

    Data my_data; 

    fTree->Branch("beta",           &my_data.beta,          "beta/D"); 
    fTree->Branch("energy",         &my_data.energy,        "energy/D"); 
    fTree->Branch("p_accept",       &my_data.p_accept,      "p_accept/D"); 
    fTree->Branch("theta",          &my_data.theta,         "theta/D"); 
    fTree->Branch("frob_distance",  &my_data.frob_distance, "frob_distance/D"); 
    fTree->Branch("time",           &my_data.time,          "time/D"); 

    for (const auto& evt : fData) {
        //copy this event into data
        my_data = evt; 

        //save this event
        fTree->Fill(); 
    }

    //save this tree to the file 
    fTree->Write(); 
    
    //close & save the file 
    fFile->Close(); 
    fTree = nullptr; 
}
//_____________________________________________________________________________________________
void TreeWriter::WriteLines(const std::vector<TreeWriter::Data>& data_in)
{
    fWriteMutex.lock();
    fData.insert( fData.end(), data_in.begin(), data_in.end() ); 
    fWriteMutex.unlock(); 
}
//_____________________________________________________________________________________________
void TreeWriter::WriteLine(double beta, double energy, double p_accept, double theta, double frob_distsance, double time)
{   
    fWriteMutex.lock();
    fData.emplace_back(beta, energy, p_accept, theta, frob_distsance, time); 
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
