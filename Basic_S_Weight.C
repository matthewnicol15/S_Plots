{

  //////////////////////////////////////////////////////////////////////////////
  ////Define variables for naming and limits ///////////
  //////////////////////////////////////////////////////////////////////////////

  // Information for canvas and histogram name
  // ostringstream Data;
  // ostringstream Quantity;
  // ostringstream Date;
  // ostringstream Version;
  // ostringstream Output_File_Name;

  // Setting the strings for canvas name
  // Data<<"RGA_Spring2019_Inbending_dst_Tree_04";
  // Quantity<<"Total";
  // Date<<"19112021";
  // Version<<"01";

  // Output_File_Name<<"/media/mn688/Elements1/PhD/Analysis_Output/S_Weight_Strangeness_Analysis_"<<Data.str().c_str()<<"_"<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<"_"<<Version.str().c_str()<<".root";


  //////////////////////////////////////////////////////////////////////////////
  //// Getting input file and histograms    ////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Input file
  TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/Strangeness_Analysis_RGA_Spring2019_Inbending_dst_Tree_Total_22112021_01.root");

  TFile *f2 = new TFile("/shared/storage/physhad/JLab/mn688/Analysis_Output/Strangeness_RGA_SPRING_2019_Inbending_eFD_Kp_201021_04_part2_Total_221121_03.root");

  // Getting the multidimensional histogram plots
  TH2D *hmiss_1_a__S1_kp_1 = (TH2D*)f1->Get("hmiss_1_a__S1_kp_1");
  // Get the number of x bins
  Int_t maxbin = hmiss_1_a__S1_kp_1->GetNbinsY();
  TH3F *hmiss_s2_a__S2_kp_1__S2_kp_2 = (TH3F*)f1->Get("hmiss_s2_a__S2_kp_1__S2_kp_2");
  Int_t zbinmax = hmiss_s2_a__S2_kp_1__S2_kp_2->GetNbinsZ();

  hmiss_s2_a__S2_kp_1__S2_kp_2->Project3D("y")->Draw();

  // Get sideband subtracted result
  TH1F *hmiss_1_a_sig = (TH1F*)f2->Get("hmiss_1_a_sig");
  TH1F *hmiss_1_a_result = (TH1F*)f2->Get("hmiss_1_a_result");

  // Make array for all the 1D x projections
  TH1D *h_projectionx[100];

  //////////////////////////////////////////////////////////////////////////////
  //// Creating and fitting functions    ///////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Define functions for fitting kaon calculated mass
  // Functions for strangeness 1 - kaon 1
  TF1 *func1 = new TF1("func1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2 = new TF1("func2","gaus(0)",0.36,0.7);
  TF1 *func3 = new TF1("func3","pol3(0)",0.36,0.7);
  TF1 *func4 = new TF1("func4","gaus(0)",0.36,0.7);
  TF1 *func5 = new TF1("func5","gaus(0) + gaus(3)",0.36,0.7);

  // Functions for strangeness 2 - kaon 1
  TF1 *func1_s2_kp1 = new TF1("func1_s2_kp1","gaus(0) + pol3(3) + gaus(7)",0.36,0.7);
  TF1 *func2_s2_kp1 = new TF1("func2_s2_kp1","gaus(0)",0.36,0.7);
  TF1 *func3_s2_kp1 = new TF1("func3_s2_kp1","pol3(0)",0.36,0.7);
  TF1 *func4_s2_kp1 = new TF1("func4_s2_kp1","gaus(0)",0.36,0.7);
  TF1 *func5_s2_kp1 = new TF1("func5_s2_kp1","gaus(0) + gaus(3)",0.36,0.7);


  // Setting parameters before fitting
  // Strangeness 1 - kaon 1
  func1->SetParameters(hmass_S1_kp_1->GetMaximum(),0.493,0.02); // Amplitude, mean, sigma for firs gauss
  func1->SetParameter(7,hmass_S1_kp_1->GetMaximum() / 2); // amplitude for second gauss
  func1->SetParameter(8,0.493); // mean for second gauss
  func1->SetParameter(9,0.02); // sigma for second gauss
  // Setting parameter limits before fitting
  func1->SetParLimits(0,hmass_S1_kp_1->GetMaximum() / 3,hmass_S1_kp_1->GetMaximum()); // amplitude for first gauss
  func1->SetParLimits(1,0.480,0.505); // mean for first gauss
  func1->SetParLimits(2,0.005,0.03); // sigma for first gauss
  func1->SetParLimits(7,hmass_S1_kp_1->GetMaximum() / 3,hmass_S1_kp_1->GetMaximum()); // amplitude for second gauss
  func1->SetParLimits(8,0.480,0.505); // mean for second gauss
  func1->SetParLimits(9,0.005,0.05); // sigma for second gauss

  // Strangeness 2 - kaon 1
  func1_s2_kp1->SetParameters(hmiss_s2_a__S2_kp_1__S2_kp_2_y->GetMaximum() / 2,0.493,0.02);
  func1_s2_kp1->SetParameter(7,hmiss_s2_a__S2_kp_1__S2_kp_2_y->GetMaximum() / 2);
  func1_s2_kp1->SetParameter(8,0.493);
  func1_s2_kp1->SetParameter(9,0.02);
  // Setting parameter limits before fitting
  func1_s2_kp1->SetParLimits(0,hmiss_s2_a__S2_kp_1__S2_kp_2_y->GetMaximum() / 3,hmiss_s2_a__S2_kp_1__S2_kp_2_y->GetMaximum());
  func1_s2_kp1->SetParLimits(1,0.480,0.505);
  func1_s2_kp1->SetParLimits(2,0.005,0.03);
  func1_s2_kp1->SetParLimits(7,hmiss_s2_a__S2_kp_1__S2_kp_2_y->GetMaximum() / 3,hmiss_s2_a__S2_kp_1__S2_kp_2_y->GetMaximum());
  func1_s2_kp1->SetParLimits(8,0.480,0.505);
  func1_s2_kp1->SetParLimits(9,0.005,0.05);

  // Fitting functions and getting parameters
  // Strangeness 1 - kaon 1
  hmass_S1_kp_1->Fit("func1","RB");
  func2->FixParameter(0, func1->GetParameter(0));
  func2->FixParameter(1, func1->GetParameter(1));
  func2->FixParameter(2, func1->GetParameter(2));
  func3->FixParameter(0, func1->GetParameter(3));
  func3->FixParameter(1, func1->GetParameter(4));
  func3->FixParameter(2, func1->GetParameter(5));
  func3->FixParameter(3, func1->GetParameter(6));
  func4->FixParameter(0, func1->GetParameter(7));
  func4->FixParameter(1, func1->GetParameter(8));
  func4->FixParameter(2, func1->GetParameter(9));
  func5->FixParameter(0, func1->GetParameter(0));
  func5->FixParameter(1, func1->GetParameter(1));
  func5->FixParameter(2, func1->GetParameter(2));
  func5->FixParameter(3, func1->GetParameter(7));
  func5->FixParameter(4, func1->GetParameter(8));
  func5->FixParameter(5, func1->GetParameter(9));

  // Strangeness 2 - kaon 1
  hmiss_s2_a__S2_kp_1__S2_kp_2_y->Fit("func1_s2_kp1","RB");
  func2_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(0));
  func2_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(1));
  func2_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(2));
  func3_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(3));
  func3_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(4));
  func3_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(5));
  func3_s2_kp1->FixParameter(3, func1_s2_kp1->GetParameter(6));
  func4_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(7));
  func4_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(8));
  func4_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(9));
  func5_s2_kp1->FixParameter(0, func1_s2_kp1->GetParameter(0));
  func5_s2_kp1->FixParameter(1, func1_s2_kp1->GetParameter(1));
  func5_s2_kp1->FixParameter(2, func1_s2_kp1->GetParameter(2));
  func5_s2_kp1->FixParameter(3, func1_s2_kp1->GetParameter(7));
  func5_s2_kp1->FixParameter(4, func1_s2_kp1->GetParameter(8));
  func5_s2_kp1->FixParameter(5, func1_s2_kp1->GetParameter(9));


  //////////////////////////////////////////////////////////////////////////////
  //// Calculating scaling factor    //////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  hmiss_1_a__S1_kp_1->ProjectionY();
  Int_t Lower_end = hmiss_1_a__S1_kp_1_py->FindBin(func1->GetParameter(1) - 2 * func1->GetParameter(2));
  Int_t Upper_end = hmiss_1_a__S1_kp_1_py->FindBin(func1->GetParameter(1) + 2 * func1->GetParameter(2));

  hmiss_s2_a__S2_kp_1__S2_kp_2->Project3D("y");
  Int_t Peak_Lower_end_s2 = hmiss_s2_a__S2_kp_1__S2_kp_2->FindBin(func1_s2_kp1->GetParameter(1) - 2 * func1_s2_kp1->GetParameter(2));
  Int_t Peak_Upper_end_s2 = hmiss_s2_a__S2_kp_1__S2_kp_2->FindBin(func1_s2_kp1->GetParameter(1) + 2 * func1_s2_kp1->GetParameter(2));
  Int_t Back_Left_Lower_end_s2 = hmiss_s2_a__S2_kp_1__S2_kp_2->FindBin(func1_s2_kp1->GetParameter(1) - 7 * func1_s2_kp1->GetParameter(2));
  Int_t Back_Left_Upper_end_s2 = hmiss_s2_a__S2_kp_1__S2_kp_2->FindBin(func1_s2_kp1->GetParameter(1) - 5 * func1_s2_kp1->GetParameter(2));
  Int_t Back_Right_Lower_end_s2 = hmiss_s2_a__S2_kp_1__S2_kp_2->FindBin(func1_s2_kp1->GetParameter(1) + 5 * func1_s2_kp1->GetParameter(2));
  Int_t Back_Right_Upper_end_s2 = hmiss_s2_a__S2_kp_1__S2_kp_2->FindBin(func1_s2_kp1->GetParameter(1) + 7 * func1_s2_kp1->GetParameter(2));

 Double_t Scaling_Factor;

  // Loop over the x bins to determine scaling factors
  for(Int_t bin = 1; bin < 100; bin++){

  cout<<bin<<endl;
    // Check the values are above zero
    if(func5->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) > 0 && func3->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) > 0){

      // Determine the scaling factor looking at the signal and background functions
      Scaling_Factor = func5->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) / (func3->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) + func5->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)));

      // Get the x projection for the current bin
      h_projectionx[bin] = (TH1D*)hmiss_1_a__S1_kp_1->ProjectionX("",bin,bin)->Clone("projection_x");
      h_projectionx[bin]->Scale(Scaling_Factor);
    }
    else{
      h_projectionx[bin]->Scale(0);
    }
  }
  for(Int_t bin = 2; bin < 100; bin++){
    cout<<"bin"<<endl;
    h_projectionx[Lower_end]->Add(h_projectionx[bin]);
  }

  //////////////////////////////////////////////////////////////////////////////
  //// Creating canvases and drawing plots    //////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // canvas for projections
  // auto *c1 = new TCanvas("c1","original and backgroud subtracted",800,800);
  // c1->cd();
  // hmiss_1_a__S1_kp_1->ProjectionX()->Draw();
  // // h_projectionx[Lower_end]->Scale(1.3);
  // h_projectionx[Lower_end]->SetLineColor(kRed);
  // // hmiss_1_a_result->Scale(1.3);
  // hmiss_1_a_result->SetLineColor(kGreen);
  // h_projectionx[Lower_end]->Draw("hist,same");
  // hmiss_1_a_result->Draw("hist,same");
  //
  //
  // auto *c2 = new TCanvas("c2","sideband",800,800);
  // c2->cd();
  // hmiss_1_a_sig->Draw("hist");
  // hmiss_1_a_result->Draw("hist,same");
  // h_projectionx[Lower_end]->Draw("hist,same");
  //
  // auto *c3 = new TCanvas("c3","pion peak",800,800);
  // c3->cd();
  // // hmiss_1_a__S1_kp_1->ProjectionY()->Draw();
  // hmiss_1_a__S1_kp_1->ProjectionX();
  // // hmiss_1_a__S1_kp_1->ProjectionY("",hmiss_1_a__S1_kp_1_px->FindBin(0.8),hmiss_1_a__S1_kp_1_px->FindBin(0.99))->Draw();
  // hmiss_1_a__S1_kp_1->ProjectionY()->Draw();
  // func1->Draw("same");
  // func2->Draw("same");
  // func3->Draw("same");
  // func4->Draw("same");
  // func5->Draw("same");
  //
  auto *c4 = new TCanvas("c4","strangeness 2 peak",800,800);
  c4->cd();
  hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX("_px",Peak_Lower_end_s2, Peak_Upper_end_s2,0,zbinmax)->Draw();
  hmiss_s2_a__S2_kp_1__S2_kp_2_px->Draw();
  TH1D *S2_peak = hmiss_s2_a__S2_kp_1__S2_kp_2_px->Clone("S2_peak");

  auto *c5 = new TCanvas("c5","strangeness 2 background",800,800);
  c5->cd();
  hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX("_px",Back_Left_Lower_end_s2, Back_Left_Upper_end_s2,0,zbinmax);
  TH1D *S2_background_left = hmiss_s2_a__S2_kp_1__S2_kp_2_px->Clone("S2_back_left");
  hmiss_s2_a__S2_kp_1__S2_kp_2->ProjectionX("_px",Back_Right_Lower_end_s2, Back_Right_Upper_end_s2,0,zbinmax);
  TH1D *S2_background_Right = hmiss_s2_a__S2_kp_1__S2_kp_2_px->Clone("S2_back_Right");
  TH1D *S2_background_Total = (TH1D*)S2_background_left->Clone("S2_background_Total");
  S2_background_Total->Add(S2_background_Right);
  S2_background_Total->Draw();

  TH1D *S2_Result = (TH1D*)S2_peak->Clone("S2_Result");
  S2_Result->Add(S2_background_Total,-1);

  auto *c6 = new TCanvas("c6","strangeness 2 after sideband subtraction",800,800);
  c6->cd();
  // S2_peak->Draw();
  S2_Result->SetLineColor(kRed);
  S2_Result->Draw();
}
