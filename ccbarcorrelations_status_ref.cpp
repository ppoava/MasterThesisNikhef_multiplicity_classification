// ccbarcorrelations_status.cpp                                                                                                       
// Script that used PYTHIA to simulate collision events according to the settings given in the pythiasettings_Hard_Low_cc.cmnd file.   

// The variables saved are necessary to create azimuthal angular correlation plots for beauty hadrons, which are defined as follows:   
// pT = transverse momentum                                                                                                            
// eta = pseudorapidity                                                                                                                
// phi = azimuthal angle                                                                                                               
// ID = PDG ID of particle                                                                                                             
// mother = index of mother particle                                                                                                 
// motherID = PDG ID of mother particle                                                                                             

// Other variables are saved to be used to create a few histograms that provide quick checks on the quality of the simulated output

//C++ libraries we are using
#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>
#include <vector>
// PYTHIA library
#include "Pythia8/Pythia.h"
// ROOT libraries 
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TTree.h"
// include for mac compilation (getpid)
#include <unistd.h>

#define PI 3.14159265
using namespace std;
using namespace Pythia8;

// Simple function that checks if a particle has any charm content
bool IsCharm( Int_t particlepdg){
	Int_t pdg = abs(particlepdg);
	pdg /= 10; // Last digit does not have to do with quark content
	if(pdg % 10 == 4) return true; // 3rd quark
	pdg /= 10;
	if(pdg % 10 == 4) return true; // 2nd quark
	pdg /= 10;
	if(pdg % 10 == 4) return true; // 1st quark
	return false;
	}

// Returns delta phi in range -pi/2 3pi/2
Double_t DeltaPhi(Double_t phi1, Double_t phi2) {
	return fmod(phi1-phi2+2.5*PI,2*PI)-0.5*PI;
	}

int main(int argc, char** argv) {

	if(argc != 4) { 		
	        cout<<"Error in the number of arguments provided"<<endl;
		cout<<"Expected input: ./ccbarcorrelations_status output_name random_number1 random_number2"<<endl;
		cout<<"Terminating program"<<endl;
		return 0;
	}

	// For larger simulations it is essential to know the duration of the execution of this script
	// Start time here
	auto start = chrono::high_resolution_clock::now();
	
	// Create output file
	TFile* output = new TFile(argv[1],"RECREATE");
	if(!output->IsOpen()){
		cout<<"Error: File "<<argv[1]<<" already exists terminating program!"<<endl;
		return 1;
	}
	
	// Define TTree that will contain event information
	TTree *tree = new TTree("tree","ccbar correlations");
	
	// Define variables that will be generated
	Double_t pT, eta, y, phi, pTTrigger, pTAssociate, charge, DeltaPhiDD, status, mother, motherID;
	Int_t id, idCharm, MULTIPLICITY, nMPIs, sph, nJets, nEvents, charmness;
	
	// Define vectors that contain event-level information on the particles
	vector<Int_t> vID;
	vector<Double_t> vPt;
	vector<Double_t> vEta; // pseudorapidity
	vector<Double_t> vY; // rapidity
	vector<Double_t> vPhi;
	vector<Double_t> vCharge;
	vector<Double_t> vStatus;
	vector<Double_t> vMother1;
	vector<Double_t> vMotherID;
	
	// Setting up tree branches
	// The histograms are used to provide quick checks on the output
	tree->Branch("ID",&vID);
	tree->Branch("PT",&vPt);
	tree->Branch("ETA",&vEta);
	tree->Branch("Y",&vY);
	tree->Branch("PHI",&vPhi);
	tree->Branch("CHARGE",&vCharge);
	tree->Branch("STATUS",&vStatus);
	tree->Branch("MOTHER",&vMother1);
	tree->Branch("MOTHERID",&vMotherID);
	tree->Branch("MULTIPLICITY",&MULTIPLICITY,"x/I");
	tree->Branch("nMPIs", &nMPIs, "nMPIs/I"); // make caps
	// tree->Branch("SPH", &sph, "sph/D");
	tree->Branch("nJets", &nJets, "nJets/I"); // make caps
	TH1D* hMULTIPLICITY = new TH1D("hMULTIPLICITY","Multiplicity",301,-0.5,300.5);
	TH1D *hNMPIs = new TH1D("hNMPIs", "number of MPIs in event", 50, 0, 50);
	// TH1D *hSph = new TH1D("hSph", "sphericity of event", 100, 0, 1);
	TH1D *hNJets = new TH1D("hNJets", "number of jets in event", 20, 0, 20);
	TH1D* hidCharm = new TH1D("hidCharm","PDG Codes for Charm hadrons",12000,-6000,6000);
	TH1D* hPtTrigger = new TH1D("hPtTrigger","p_{T} for triger D^{+} ",50,0,10);
	TH1D* hPtAssociate = new TH1D("hPtAssociate", "p_{T} for associate D^{-}",50,0,10);
	TH1D* hDeltaPhiDD = new TH1D("hDeltaPhiDD","D^{+}D^{-} correlations",100,-PI/2,3*PI/2);
	TH1D* hCharmPart = new TH1D("hCharmPart", "Charm Particles Per Event",200,-0.5,200.5);
	
	// Kinematics constraints
	const Double_t pTmin = 0.15; // minimum pT
	const Double_t etamax = 4.; // maximum pseudorapidity (absolute value)
	
	// Get PYTHIA
	Pythia pythia;

	// Get event analysis tools
	// Used for sphericity and slowjet (to analyse jet structure of events)
	// PROBLEM: sphericity returns 0 always? Seems to be related to not being able to apply kinematic cuts to selection?
	// Paramters for Sphericity:
	// int power_sphericity = 2;
	// int select_sphericity = 3;		   // only final-state charged particles
	// Parameters for SlowJet:
	int power_slowJet = -1;		   // anti-kt
	double R = 0.4;		   // jet radius
	double pTjetMin = 5.0; // minimum pT
	double etaMax = 4.0;   // max pseudorapidity
	int select_slowJet = 3;		   // only final-state charged particles
	int massSet = 2;	   // invariant mass scheme
	SlowJetHook *sjHookPtr = nullptr;
	bool useFJcore = true;
	bool useStandardR = true;
	// Initialise
	// Sphericity sphericity(power_sphericity, select_sphericity);
	SlowJet slowJet(power_slowJet, R, pTjetMin, etaMax, select_slowJet, massSet, sjHookPtr, useFJcore, useStandardR);
	
	// Simulation settings from pythiasettings_Hard_Low_cc.cmnd
	// The settings used are documented in that file
	pythia.readFile("pythiasettings_Hard_Low_cc.cmnd");
	nEvents = pythia.mode("Main:numberOfEvents");
	
	// Create a random seed so that the outcome will be truly random
        Int_t proccessid = getpid();
        int seedMod1 = std::stoi(argv[2]);
        int seedMod2 = std::stoi(argv[3]);
        string seedstr = "Random:seed = "+std::to_string((time(0)+proccessid+seedMod1+seedMod2)%900000000);
        pythia.readString("Random:setSeed = on");
        pythia.readString(seedstr);

	// Initializing simulation
	pythia.init(); 
	cout<<"Generating "<<nEvents<<" events!"<<endl;
	
	// Event loop
	for(int iEvent = 0; iEvent<nEvents; iEvent++){
		if(!pythia.next()) continue;
		if (!slowJet.analyze(pythia.event)) {
			std::cout << "ERROR: Failed to do slowJet in event" << std::endl;
		}
		int nPart = pythia.event.size(); // Number of particles produced in this event
		nMPIs = pythia.info.nMPI();
		// sph = sphericity.sphericity();
		// std::cout << "sphericity in event = " << sph << std::endl;
		nJets = slowJet.sizeJet();
		MULTIPLICITY = 0; // Initialiazing for multiplicity plot
		charmness = 0; // Intialiazing for charm production plot
		
		// Initializing vectors for event-level information
		vID.clear();
		vPt.clear();
		vEta.clear();
		vY.clear();
		vPhi.clear();
		vCharge.clear();
		vStatus.clear();
		vMother1.clear();
		vMotherID.clear();

		// Particle loop
		for(int iPart = 0; iPart<nPart; iPart++){
			const Particle &particle = pythia.event[iPart];
			if(!particle.isFinal()) continue; // Skip if the particle is not at its final state
			
			id = particle.id();
			pT = particle.pT();
			eta = particle.eta();
			y = particle.y();
			phi = particle.phi();
			charge = particle.charge();
			status = static_cast<Double_t> (particle.status());
			mother = static_cast<Double_t> (particle.mother1());
			motherID = static_cast<Double_t> (pythia.event[particle.mother1()].id());
			
			// Kinematics check
			if(pT < pTmin || abs(eta) > etamax ) continue;
			
			if(abs(id) == 11 || abs(id) == 13 || abs(id) == 211 || abs(id) == 321 || abs(id) == 2212) { // Needs to be a primary particle                                                                                                                             
                if(81 <= status && status >= 89) { // Needs to be prompt                                                                                                                                                                                                
                    MULTIPLICITY++;
                }
            }
			
			if(!IsCharm(id)) continue;
				idCharm = id;
				hidCharm->Fill((Double_t) id);
				charmness++;
				// Filling vectors
				vID.push_back(id);
				vPt.push_back(pT);
				vEta.push_back(eta);
				vY.push_back(y);
				vPhi.push_back(phi);
				vCharge.push_back(charge);
				vStatus.push_back(status);
				vMother1.push_back(mother);
				vMotherID.push_back(motherID);
				
				// Creating D+D- correlation plots as check
				if(id == 421){ // D+ meson triger
				
					for(int jPart = 0; jPart<nPart; jPart++){
						const Particle particleB = pythia.event[jPart];

						if(jPart == iPart) continue; // So we do not correlate a particle with itself.
						Int_t associate_id = particleB.id();
						Double_t associate_pT = particleB.pT();
						Double_t associate_eta = particleB.eta();

						if(associate_pT < pTmin  || abs(associate_eta) > etamax) continue;
						if(associate_id == -421){ // D- meson associate
							pTTrigger = pT;
							DeltaPhiDD = DeltaPhi(phi,particleB.phi());
							hPtTrigger->Fill(pT);
							hPtAssociate->Fill(associate_pT);
							hDeltaPhiDD->Fill(DeltaPhiDD);					
						}	
						else continue;	
					} // End of B meson triger particle loop.
				} // End of B meson triger
		} // 1st particle loop

		// Don't consider events that don't pass kinematic cuts
		if (MULTIPLICITY == 0)
		{
			continue;
		}

		hMULTIPLICITY->Fill((Double_t) MULTIPLICITY);
		hCharmPart->Fill((Double_t) charmness);

		// In order not to fill tree with empty vectors.
		if(vID.empty() || vPt.empty() || vEta.empty() || vY.empty() || vPhi.empty() || vCharge.empty() || vStatus.empty() || vMother1.empty() || vMotherID.empty() ) continue; 
		tree->Fill();
	} // End of event loop
	
	// Write output and close it
	output->Write();
	cout<<"File has been created and its name is: "<<output->GetName()<<endl;
	output->Close();
	
	// Return duration of the simulation
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start);
        cout << "This script took " << duration.count() << " minutes to run." << endl;
	
	return 0;
}




