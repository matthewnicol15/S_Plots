#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TBenchmark.h>
#include <vector>
#include <TCanvas.h>
#include <chrono>
#include <sstream>
#include <string>

// Macro name
void Strangeness_1_S_Weight(){


  auto start = std::chrono::high_resolution_clock::now();

  // Read file with information on vectors
  // gROOT->ProcessLine(".L ~/work/Macros/Q_Weight/Loader.C+");
  gROOT->ProcessLine(".L /media/mn688/Elements1/PhD/Macros/Loader.C+");

  // Read input root file and assign it to 'f'
  TFile *fin = new TFile("/media/mn688/Elements1/PhD/Trees/Dibaryon/RGA/Strangeness_1/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_260521_01.root");
  // Read TTree within root file and assign it to 't1'
  TTree *t1 = (TTree*)fin->Get("RGA_Skim4_Tree_260521_01");

  // Recording the calculated mass of kaons
  Double_t mass_kp;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Functions to fit S weight plots
  // auto* func4=new TF1("func4","[0]*exp(-pow(x-[1],2)/(2*pow([2],2))) + [3]*exp(-pow(x-[4],2)/(2*pow(3.2*[2],2))) + pol1(5)",0.36,0.65);
  auto* func1=new TF1("func1","gaus(0) + pol4(3)",0.36,0.65);
  auto* func5=new TF1("func5","gaus(0)",0.36,0.65);
  auto* func6=new TF1("func6","pol4(0)",0.36,0.65);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  // Creating components to read from TTree
  // Set any vectors to 0 before reading from the TTree
  // Event information
  TLorentzVector *readbeam=NULL;  // Information on the beam
  TLorentzVector *readtarget=NULL; // Information on the target
  Double_t start_time; // Event start time
  Int_t readrunno; // Run number
  Int_t readeventno; // Event number
  Int_t readtriggerno; // Trigger number
  // Number of given particle or charge track in each event
  Int_t readchargetracks; // Number of positive or negative charge tracks
  Int_t readothertracks; // Number of particles excluding p, e^-, pions or kaons
  Int_t region; // which region the particles go in (FT, FD, CD)


  // Particle information
  vector<TLorentzVector> *v_p4=0;   // 4-vectors for the detected particles
  vector<TLorentzVector> *v_vertex=0;   // Vertex information for particles
  vector<double> *v_path=0;   // Measured path of particles
  vector<double> *v_time=0;   // Measured time of flight of particles
  vector<double> *v_beta=0;   // Beta measured from FTOF
  vector<double> *v_energy=0;   // Energy of particle
  vector<double> *v_charge=0;   // Charge of particle from Drift Chambers
  vector<double> *v_PID=0;   // PID from FTOF
  vector<double> *v_chi2PID=0;   // Chi^2 of the PID
  vector<Int_t> *v_region=0; // region particle goes in

  // Setting the branch addresses to read from original TTree
  t1->SetBranchAddress("p4",&v_p4);
  t1->SetBranchAddress("vertex",&v_vertex);
  t1->SetBranchAddress("path",&v_path);
  t1->SetBranchAddress("beta",&v_beta);
  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);
  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&v_energy);
  t1->SetBranchAddress("charge",&v_charge);
  t1->SetBranchAddress("PID",&v_PID);
  t1->SetBranchAddress("chi2PID",&v_chi2PID);
  t1->SetBranchAddress("region",&v_region);
  t1->SetBranchAddress("time",&v_time);
  t1->SetBranchAddress("eventno",&readeventno);
  t1->SetBranchAddress("runno",&readrunno);
  t1->SetBranchAddress("triggerno",&readtriggerno);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Getting particle database to use for masses
  auto db=TDatabasePDG::Instance();

  // Create histograms here

  auto* hmass_kp=new TH1F("hmass_kp","K^{+} mass;M(K^{+});Counts",100,0.2,0.8);
  auto* hmiss_1=new TH1F("hmiss_1","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_1_sig=new TH1F("hmiss_1_sig","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);
  auto* hmiss_1_back=new TH1F("hmiss_1_back","MM(e' K^{+});MM(e' K^{+}) [GeV];Counts",800,-1,3);

  auto* h_delta_beta_kp=new TH2F("h_delta_beta_kp","",200,0,11,200,-1,1);
  auto* h_delta_beta_pr=new TH2F("h_delta_beta_pr","",200,0,11,200,-1,1);
  auto* h_delta_beta_pim=new TH2F("h_delta_beta_pim","",200,0,11,200,-1,1);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Create all the various vectors for the different particles

  // Create vectors of TLorentzVectors to store information of
  // all particles of a given type (important when you have more than 1
  // of any particle in the same event)
  vector<TLorentzVector> v_el;  // e^-
  vector<TLorentzVector> v_pim; // pi^-
  vector<TLorentzVector> v_pr; // protons
  vector<TLorentzVector> v_kp; // K^+
  vector<TLorentzVector> v_othertracks; // Any other particles are assigned to this

  // TLorentzVectors for individual particles
  TLorentzVector el;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;


  // Vertex information for different particle types
  vector<TLorentzVector> v_vertex_el; // e^-
  vector<TLorentzVector> v_vertex_pim; // pi^-
  vector<TLorentzVector> v_vertex_pr; // protons
  vector<TLorentzVector> v_vertex_kp; // K^+
  TLorentzVector vertex_el;
  TLorentzVector vertex_pim;
  TLorentzVector vertex_pr;
  TLorentzVector vertex_kp;

  // These are used to define the missing masses later
  TLorentzVector missall;
  TLorentzVector miss1;
  TLorentzVector lambda;

  // After information is read from the TTree, particles are identified using
  // PID and assigned to that particle type
  // e^-
  vector<Double_t> v_beta_tof_el;  // Beta measured
  Double_t beta_tof_el;
  vector<Double_t> v_P_el;  // Momentum measured
  Double_t P_el;
  vector<Double_t>v_path_el; // Path measured
  Double_t path_el;
  vector<Double_t>v_TOF_el;  // TOF measured
  Double_t TOF_el;
  vector<Double_t> v_beta_calc_el;  // Beta calculated
  Double_t beta_calc_el;
  vector<Double_t> v_delta_beta_el;  // Beta calculated - beta measured
  Double_t delta_beta_el;
  vector<Double_t> v_vertex_time_el;  // Vertex time calculated
  Double_t vertex_time_el;
  vector<Double_t> v_region_el; // region hit
  Double_t region_el;

  // pi^-
  vector<Double_t> v_beta_tof_pim;  // Beta measured
  Double_t beta_tof_pim;
  vector<Double_t> v_P_pim;  // Momentum measured
  Double_t P_pim;
  vector<Double_t>v_path_pim; // Path measured
  Double_t path_pim;
  vector<Double_t>v_TOF_pim;  // TOF measured
  Double_t TOF_pim;
  vector<Double_t> v_beta_calc_pim;  // Beta calculated
  Double_t beta_calc_pim;
  vector<Double_t> v_delta_beta_pim;  // Beta calculated - beta measured
  Double_t delta_beta_pim;
  vector<Double_t> v_vertex_time_pim;  // Vertex time calculated
  Double_t vertex_time_pim;
  vector<Double_t> v_region_pim; // region hit
  Double_t region_pim;

  // protons
  vector<Double_t> v_beta_tof_pr;  // Beta measured
  Double_t beta_tof_pr;
  vector<Double_t> v_P_pr;  // Momentum measured
  Double_t P_pr;
  vector<Double_t>v_path_pr; // Path measured
  Double_t path_pr;
  vector<Double_t>v_TOF_pr;  // TOF measured
  Double_t TOF_pr;
  vector<Double_t> v_beta_calc_pr;  // Beta calculated
  Double_t beta_calc_pr;
  vector<Double_t> v_delta_beta_pr;  // Beta calculated - beta measured
  Double_t delta_beta_pr;
  vector<Double_t> v_vertex_time_pr;  // Vertex time calculated
  Double_t vertex_time_pr;
  vector<Double_t> v_region_pr; // region hit
  Double_t region_pr;

  // K^+
  vector<Double_t> v_beta_tof_kp;  // Beta measured
  Double_t beta_tof_kp, beta_tof_kp_other;
  vector<Double_t> v_P_kp;  // Momentum measured
  Double_t P_kp;
  vector<Double_t>v_path_kp; // Path measured
  Double_t path_kp;
  vector<Double_t>v_TOF_kp;  // TOF measured
  Double_t TOF_kp;
  vector<Double_t> v_beta_calc_kp;  // Beta calculated
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;  // Beta calculated - beta measured
  Double_t delta_beta_kp;
  vector<Double_t> v_vertex_time_kp;  // Vertex time calculated
  Double_t vertex_time_kp;
  vector<Double_t> v_region_kp; // region hit
  Double_t region_kp;
  vector<Double_t> v_mass_kp; // region hit


  Double_t c=30;  // Speed of light used for calculating vertex time


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop over the events

  // Reads the total number of entries in the TTree
  // Long64_t nentries = t1->GetEntries();
  // You can just run over a set number of events for fast analysis
  Long64_t nentries = 100000;
  cout<<t1->GetEntries()<<endl;
  // This is used to print out the percentage of events completed so far
  Int_t Percentage = nentries/100;
  // This loops over all the entries in the TTree
  for(Long64_t i = 0; i < nentries; i++){
    // Printing the event number every 100 events
    if((i/10)%10==0 && i%10==0){
      cout<<i<<endl;
      // auto finish = std::chrono::high_resolution_clock::now();
      // std::chrono::duration<double> elapsed = finish - start;
      // std::cout << i<<" Elapsed time: " << elapsed.count();
    }

    // Get entry 'i' from the TTree
    t1->GetEntry(i);


    // All the vectors must be cleared at the start of each event entry
    // e^-
    v_el.clear();
    v_beta_tof_el.clear();
    v_P_el.clear();
    v_path_el.clear();
    v_TOF_el.clear();
    v_beta_calc_el.clear();
    v_delta_beta_el.clear();
    v_vertex_time_el.clear();
    v_vertex_el.clear();
    v_region_el.clear();

    // pi^-
    v_pim.clear();
    v_beta_tof_pim.clear();
    v_P_pim.clear();
    v_path_pim.clear();
    v_TOF_pim.clear();
    v_beta_calc_pim.clear();
    v_delta_beta_pim.clear();
    v_vertex_pim.clear();
    v_vertex_time_pim.clear();
    v_vertex_pim.clear();
    v_region_pim.clear();

    // protons
    v_pr.clear();
    v_beta_tof_pr.clear();
    v_P_pr.clear();
    v_path_pr.clear();
    v_TOF_pr.clear();
    v_beta_calc_pr.clear();
    v_delta_beta_pr.clear();
    v_vertex_pr.clear();
    v_vertex_time_pr.clear();
    v_vertex_pr.clear();
    v_region_pr.clear();

    // K^+
    v_kp.clear();
    v_beta_tof_kp.clear();
    v_P_kp.clear();
    v_path_kp.clear();
    v_TOF_kp.clear();
    v_beta_calc_kp.clear();
    v_delta_beta_kp.clear();
    v_vertex_kp.clear();
    v_vertex_time_kp.clear();
    v_vertex_kp.clear();
    v_region_kp.clear();
    v_mass_kp.clear();


    // This reads the number of particles in the current entry
    Int_t Nparticles = v_p4->size();

    // This loops over all the particles in the current entry
    for(Int_t j=0; j<Nparticles; j++){

      // Filling histogram to show the different regions hit
      // hregion->Fill(v_region->at(j));

      // Checking PID and assigning particles
      // e^-
      if(v_PID->at(j)==11){
        // Setting the 4-vector and assigning mass from PDG
        el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(11)->Mass());
        TOF_el = v_time->at(j); // Measured time
        path_el = v_path->at(j); // Measured path
        beta_tof_el = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_el = el.Rho()/(sqrt((pow(el.Rho(),2))+(pow(el.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_el = beta_calc_el-beta_tof_el;
        vertex_time_el = TOF_el - path_el / (beta_tof_el*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_el.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_el);
        region_el = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_el.push_back(el);
        v_beta_tof_el.push_back(beta_tof_el);
        v_P_el.push_back(P_el);
        v_path_el.push_back(path_el);
        v_TOF_el.push_back(TOF_el);
        v_beta_calc_el.push_back(beta_calc_el);
        v_delta_beta_el.push_back(delta_beta_el);
        v_vertex_time_el.push_back(vertex_time_el);
        v_vertex_el.push_back(vertex_el);
        v_region_el.push_back(region_el);
      }

      // pi^-
      else if(v_PID->at(j)==-211){
        // Setting the 4-vector and assigning mass from PDG
        pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(-211)->Mass());
        TOF_pim = v_time->at(j); // Measured time
        path_pim = v_path->at(j); // Measured path
        beta_tof_pim = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pim = pim.Rho()/(sqrt((pow(pim.Rho(),2))+(pow(pim.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pim = beta_calc_pim-beta_tof_pim;
        vertex_time_pim = TOF_pim - path_pim / (beta_tof_pim*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pim.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pim);
        region_pim = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pim.push_back(pim);
        v_beta_tof_pim.push_back(beta_tof_pim);
        v_P_pim.push_back(P_pim);
        v_path_pim.push_back(path_pim);
        v_TOF_pim.push_back(TOF_pim);
        v_beta_calc_pim.push_back(beta_calc_pim);
        v_delta_beta_pim.push_back(delta_beta_pim);
        v_vertex_time_pim.push_back(vertex_time_pim);
        v_vertex_pim.push_back(vertex_pim);
        v_region_pim.push_back(region_pim);
      }

      // protons
      else if(v_PID->at(j)==2212){
        // Setting the 4-vector and assigning mass from PDG
        pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(2212)->Mass());
        TOF_pr = v_time->at(j); // Measured time
        path_pr = v_path->at(j); // Measured path
        beta_tof_pr = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_pr = pr.Rho()/(sqrt((pow(pr.Rho(),2))+(pow(pr.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_pr = beta_calc_pr-beta_tof_pr;
        vertex_time_pr = TOF_pr - path_pr / (beta_tof_pr*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_pr.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_pr);
        region_pr = v_region->at(j);

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_pr.push_back(pr);
        v_beta_tof_pr.push_back(beta_tof_pr);
        v_P_pr.push_back(P_pr);
        v_path_pr.push_back(path_pr);
        v_TOF_pr.push_back(TOF_pr);
        v_beta_calc_pr.push_back(beta_calc_pr);
        v_delta_beta_pr.push_back(delta_beta_pr);
        v_vertex_time_pr.push_back(vertex_time_pr);
        v_vertex_pr.push_back(vertex_pr);
        v_region_pr.push_back(region_pr);
      }

      // K^+
      else if(v_PID->at(j)==321){
        // Setting the 4-vector and assigning mass from PDG
        kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),db->GetParticle(321)->Mass());

        TOF_kp = v_time->at(j); // Measured time
        path_kp = v_path->at(j); // Measured path
        beta_tof_kp = v_beta->at(j); // Measured beta from FTOF
        // Calculating beta from momentum and mass
        beta_calc_kp = kp.Rho()/(sqrt((pow(kp.Rho(),2))+(pow(kp.M(),2))));
        // Difference between calculated and measured beta
        delta_beta_kp = beta_calc_kp-beta_tof_kp;
        vertex_time_kp = TOF_kp - path_kp / (beta_tof_kp*c); // Calculate vertex time
        // Setting the vertex information now vertex time has been calculated
        vertex_kp.SetXYZT(v_vertex->at(j).X(), v_vertex->at(j).Y(), v_vertex->at(j).Z(), vertex_time_kp);
        region_kp = v_region->at(j);
        mass_kp = sqrt((pow(v_p4->at(j).Rho(),2) / (pow(beta_tof_kp,2))) - pow(v_p4->at(j).Rho(),2));

        // Pushing back all that iformation into the vectors
        // Again this is done so you can store information on multiple particles
        // of the same type in one place
        v_kp.push_back(kp);
        v_beta_tof_kp.push_back(beta_tof_kp);
        v_P_kp.push_back(P_kp);
        v_path_kp.push_back(path_kp);
        v_TOF_kp.push_back(TOF_kp);
        v_beta_calc_kp.push_back(beta_calc_kp);
        v_delta_beta_kp.push_back(delta_beta_kp);
        v_vertex_time_kp.push_back(vertex_time_kp);
        v_vertex_kp.push_back(vertex_kp);
        v_region_kp.push_back(region_kp);
        v_mass_kp.push_back(mass_kp);
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Here you can apply conditions on the events you want to analyse
    if(v_kp.size()!=1 || v_el.size()!=1 || v_pr.size()!=1 || v_pim.size()!=1) continue;
    {

      // Define reconstructed particles
      // Missing mass of e' K^{+}, looking for lambda ground state
      miss1 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0);
      // Missing mass of e' K^{+} p, looking for pi^{-}
      // miss2 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0);
      // Missing mass of e', looking for phi meson background
      // miss3 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_pr.at(0);
      // Invariant mass of p pi^{-}, looking for lambda ground state
      lambda = v_pr.at(0) + v_pim.at(0);
      // Missing mass of all detected particles, should have peak at 0
      missall = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - v_el.at(0) - v_kp.at(0) - v_pr.at(0) - v_pim.at(0);


      // Filling invariant and missing mass histograms
      hmass_kp->Fill(mass_kp);
      hmiss_1->Fill(miss1.M());
      // hmiss_2->Fill(miss2.M2());
      // hmiss_3->Fill(miss3.M());
      // hinv_lambda->Fill(lambda.M());
      // hmiss_mass_all->Fill(missall.M2());
      // hmiss_momentum_all->Fill(missall.Rho());

      // Checking the delta beta for each particle
      h_delta_beta_kp->Fill(v_kp.at(0).Rho(),v_delta_beta_kp.at(0));
      h_delta_beta_pr->Fill(v_pr.at(0).Rho(),v_delta_beta_pr.at(0));
      h_delta_beta_pim->Fill(v_pim.at(0).Rho(),v_delta_beta_pim.at(0));


    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Setting fit parameters and fitting

  // Setting parameteres for total function
  func1->SetParameter(0,0.8 * hmass_kp->GetMaximum());
  func1->SetParLimits(0,0.5 * hmass_kp->GetMaximum(), 1.1 * hmass_kp->GetMaximum());
  func1->SetParameter(1,0.493);
  func1->SetParameter(2,0.08);
  // func1->SetParLimits(4,800,10000);
  // func1->SetParLimits(4,800,10000);



  // Fitting and drawing functions
  hmass_kp->Fit("func1","RB");
  hmass_kp->Draw();

  // Fixing parameteres for signal function
  func5->FixParameter(0,func1->GetParameter(0));
  func5->FixParameter(1,func1->GetParameter(1));
  func5->FixParameter(2,func1->GetParameter(2));

  // Fixing parameters for background function
  func6->FixParameter(0,func1->GetParameter(3));
  func6->FixParameter(1,func1->GetParameter(4));
  func6->FixParameter(2,func1->GetParameter(5));
  func6->FixParameter(3,func1->GetParameter(6));
  func6->FixParameter(4,func1->GetParameter(7));

  // Drawing signal and background functions
  func5->SetLineColor(kBlue);
  func6->SetLineColor(kGreen);
  func5->Draw("same");
  func6->Draw("same");


  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count();
}
