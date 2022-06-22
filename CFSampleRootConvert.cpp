
#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TVector3.h>

using namespace std;


vector<string> split(string main, char delim = ' ');


void CFSampleRootConvert(string in_file_path, string out_file_path) {
	// PDG Database
	TDatabasePDG *db = new TDatabasePDG();
	TParticlePDG *p_info = new TParticlePDG();

	// Constants
	const int proton_pid = 2212;
	float p_min = 0.1;  // STAR minimum momentum acceptance
	float pt_min = 0.01;  // Issues arise in eta if pt is zero
	float qvec_pt_min = 0.2;
	float qvec_pt_max = 2.0;
	float ref_eta_max = 0.5;
	float ref2_eta_min = 0.5;
	float ref2_eta_max = 1.0;
	float ref3_eta_max = 1.0;
	float mass_qa_percent = 1.0;  // % difference in mass to output warning

	float eta_max = 1.0;  // Track selection max eta
	float pt_max_cut = 2.2;  // Track selection max pt
	float pt_min_cut = 0.3;  // Track selection min pt

	int buffer_size = 5000000;
	int split_level = 1;

	// Input file variables
	int pid;  // track particle id
	float px, py, pz, p0, mass;  // track variables

	// Output tree variables
	int event=0, refmult, refmult2, refmult3;  // event variables
	float qx, qy;  // event variables
	vector<int> pid_vec;  // track variables
	vector<float> px_vec;  // Apparently need to be declared on separate lines
	vector<float> py_vec;
	vector<float> pz_vec;  // track variables


	// Open file and define tree
	TFile f(out_file_path.data(), "RECREATE");
	TTree tr("tree", "CooperFryeSampler Data");

	// Define event branches:-------------------------------------------
	tr.Branch("event",    &event,     "event/I");
	tr.Branch("refmult",  &refmult,   "refmult/I");
	tr.Branch("refmult2", &refmult2,  "refmult2/I");
	tr.Branch("refmult3", &refmult3,  "refmult3/I");
	tr.Branch("qx",       &qx,        "qx/F");
	tr.Branch("qy",       &qy,        "qy/F");

	// particle branches:
	tr.Branch("pid",      &pid_vec, buffer_size, split_level);
	tr.Branch("px",       &px_vec,  buffer_size, split_level);
	tr.Branch("py",       &py_vec,  buffer_size, split_level);
	tr.Branch("pz",       &pz_vec,  buffer_size, split_level);


	cout<<"making .root file from .dat file..."<<endl;

	ifstream infile;
	infile.open(in_file_path);
	string line;

	while(getline(infile, line)) {
		vector<string> line_vec = split(line, ' ');
		if (line_vec[0] == "Event") {  // New event started
			event += 1;
			if (event > 1) { tr.Fill(); }  // Previous event should be recorded

			// Clear variables for next event
			qx = 0.; qy = 0.;
			refmult = 0; refmult2 = 0; refmult3 = 0;
			pid_vec.clear(); px_vec.clear(); py_vec.clear(); pz_vec.clear();

			getline(infile, line);  // Skip header line for event
		} else {  // Particle line, record particle
			pid = stoi(line_vec[0]);
			p0 = stof(line_vec[1]);
			px = stof(line_vec[2]);
			py = stof(line_vec[3]);
			pz = stof(line_vec[4]);
			mass = sqrt(p0*p0 - px*px - py*py - pz*pz);

			p_info = db->GetParticle((int)pid);
			if(!p_info) { cout << "pid: " << pid << " not in TDatabasePDG" << endl; continue; }
			else if(p_info->Mass() * 1+mass_qa_percent/100 < mass || p_info->Mass() * 1-mass_qa_percent/100 > mass) {
				cout << "pid: " << pid << " with mass " << mass << "  expected mass: " << p_info->Mass() << endl;
			}
			if(fabs((int)p_info->Charge()) == 0) continue;

			TVector3 p_mom(px, py, pz);
			if(p_mom.Mag() < p_min) continue; // Questionable

			double pt = p_mom.Perp();
			if(pt < pt_min) continue;  // Avoid bad pseudorapidity warning, should be out of eta cut anyway.

			double eta = p_mom.PseudoRapidity();

			double phi = p_mom.Phi();

			// Refmult and event plane logic

			// ref
			if(fabs(eta) < ref_eta_max) { refmult++; }

			// ref2
			if(fabs(eta) > ref2_eta_min && fabs(eta) < ref2_eta_max && fabs((int)p_info->Charge()) == 3) refmult2++;

			// ref3
			if(fabs(eta) < ref3_eta_max && fabs(pid) != proton_pid && fabs((int)p_info->Charge()) == 3) refmult3++;

			// event plane q-vector
			if(pt > qvec_pt_min && pt < qvec_pt_max && fabs(eta) < ref3_eta_max) {
				qx += cos(2*phi); qy += sin(2*phi);
			}

			// Record particle track if proton within STAR acceptance. Only record protons to save space. --> Now record anti-protons as well
			if (fabs(pid) == proton_pid && eta <= eta_max && pt >= pt_min_cut && pt <= pt_max_cut) {
				px_vec.push_back(px);
				py_vec.push_back(py);
				pz_vec.push_back(pz);
				pid_vec.push_back(pid);
			}
		}
	}// event loop end

	tr.Fill();  // Record last event


	f.cd();
	tr.Write();
	f.Close();
	cout<<" Tree written succesfully "<<endl;
}


// Emulation of Python split function. Split string into vector of strings on delim.
vector<string> split(string main, char delim) {
	if(main.size() == 0) { return {}; }

	vector<string> split_strings {""};
	for(char x:main) {
		if(x == delim) {
			if(split_strings.back() != "") {
				split_strings.push_back("");
			}
		} else {
			split_strings.back() += x;
		}
	}
	return(split_strings);
}
