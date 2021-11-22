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
  //// Getting input file, histograms and functions    /////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Input file
  TFile *f1=new TFile("/media/mn688/Elements1/PhD/Analysis_Output/RGA_Spring2019_Inbending_dst_Tree_04_Total_19112021_02.root");

  TFile *f2 = new TFile("/shared/storage/physhad/JLab/mn688/Analysis_Output/Strangeness_RGA_SPRING_2019_Inbending_eFD_Kp_201021_04_part2_Total_221121_03.root");

  // Getting the kaon mass fit functions from part 1
  // Strangeness 1 - kaon 1
  TF1 *func1 = (TF1*)f1->Get("func1");  // Total function
  TF1 *func2 = (TF1*)f1->Get("func2"); // Signal function, gauss 1
  TF1 *func3 = (TF1*)f1->Get("func3"); // Background function, pol3
  TF1 *func4 = (TF1*)f1->Get("func4"); // Signal function, gauss 2
  TF1 *func5 = (TF1*)f1->Get("func5"); // Total signal function, 2 gauss


  // Getting the multidimensional histogram plots
  TH2D *hmiss_1_a__S1_kp_1 = (TH2D*)f1->Get("hmiss_1_a__S1_kp_1");
  // Get the number of x bins
  Int_t maxbin = hmiss_1_a__S1_kp_1->GetNbinsY();


  // Get sideband subtracted result
  TH1F *hmiss_1_a_sig = (TH1F*)f2->Get("hmiss_1_a_sig");
  TH1F *hmiss_1_a_result = (TH1F*)f2->Get("hmiss_1_a_result");

  // Make array for all the 1D x projections
  TH1D *h_projectionx[100];

  //////////////////////////////////////////////////////////////////////////////
  //// Calculating scaling factor    /////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  hmiss_1_a__S1_kp_1->ProjectionY()->Draw();
  // func1->Draw("same");
  // func2->Draw("same");
  // func3->Draw("same");
  // func4->Draw("same");
  // func5->Draw("same");

 Double_t Scaling_Factor;

  // Loop over the x bins to determine scaling factors
  for(Int_t bin = 1; bin < maxbin; bin++){

  cout<<bin<<endl;
    // Check the values are above zero
    if(func5->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) > 0 && func3->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) > 0){

      // Determine the scaling factor looking at the signal and background functions
      Scaling_Factor = func2->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) / (func3->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)) + func5->Eval(hmiss_1_a__S1_kp_1_py->GetBinCenter(bin)));

      // Get the x projection for the current bin
      h_projectionx[bin] = (TH1D*)hmiss_1_a__S1_kp_1->ProjectionX("",bin,bin)->Clone("projection_x");
      h_projectionx[bin]->Scale(Scaling_Factor);
    }
    else{
      h_projectionx[bin]->Scale(0);
    }
  }
  for(Int_t bin = 2; bin < maxbin; bin++){
    cout<<"bin"<<endl;
    h_projectionx[1]->Add(h_projectionx[bin]);
  }



  // canvas for projections
  auto *c1 = new TCanvas("c1","original and backgroud subtracted",800,800);
  c1->cd();
  hmiss_1_a__S1_kp_1->ProjectionX()->Draw();
  // h_projectionx[1]->Scale(2);
  h_projectionx[1]->SetLineColor(kRed);
  hmiss_1_a_result->SetLineColor(kGreen);
  h_projectionx[1]->Draw("hist,same");
  hmiss_1_a_result->Draw("hist,same");


  auto *c2 = new TCanvas("c2","sideband",800,800);
  c2->cd();
  hmiss_1_a_sig->Draw("hist");
  hmiss_1_a_result->Draw("hist,same");






}
