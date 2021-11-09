#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooDataHist.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include <iomanip>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;

// see below for implementation
void AddModel(RooWorkspace *);
void AddData(RooWorkspace *);
void DoSPlot(RooWorkspace *);
void MakePlots(RooWorkspace *);

void S_Plot_Strangeness_1()
{

   // Create a new workspace to manage the project.
   RooWorkspace *wspace = new RooWorkspace("myWS");

   // add the signal and background models to the workspace.
   // Inside this function you will find a description of our model.
   AddModel(wspace);

   // add some toy data to the workspace
   AddData(wspace);

   // inspect the workspace if you wish
   //  wspace->Print();

   // do sPlot.
   // This will make a new dataset with sWeights added for every event.
   DoSPlot(wspace);

   // Make some plots showing the discriminating variable and
   // the control variable after unfolding.
   MakePlots(wspace);

   // cleanup
   delete wspace;
}

//____________________________________
void AddModel(RooWorkspace *ws)
{

   // Make models for signal (Higgs) and background (Z+jets and QCD)
   // In real life, this part requires an intelligent modeling
   // of signal and background -- this is only an example.

   // set range of observable
   Double_t lowRange = 0.37, highRange = 0.65;

   // make a RooRealVar for the observables
   RooRealVar *mass_kp=new RooRealVar("mass_kp", "M_{kaon}", lowRange, highRange); // Kaon mass
   RooRealVar *missing_mass=new RooRealVar("missing_mass", "missing_mass", 0., 5.); // Missing mass

   // --------------------------------------
   // make 2-d model for Z including the invariant mass
   // distribution  and an missing_mass distribution which we want to
   // unfold from QCD.
   std::cout << "make z model" << std::endl;
   // mass model for Z
   RooRealVar mkaon("mkaon", "kaon Mass", 0.493677, 0.48, 0.51);
   RooRealVar sigmakaon("sigmakaon", "Width of Gaussian", 0.001, 0, 0.1, "GeV");
   RooGaussian mkaonModel("mkaonModel", "kaon mass signal Model", *mass_kp, mkaon, sigmakaon); // Model for kaon mass signal
   // we know Z mass
   mkaon.setConstant();
   // we leave the width of the Z free during the fit in this example.

   // missing_mass model for Z.  Only used to generate toy MC.
   // the exponential is of the form exp(c*x).  If we want
   // the missing_mass to decay an e-fold every R GeV, we use
   // c = -1/R.
   RooRealVar mlambda("mlambda", "lambda Mass", 1.116, 1.0, 1.3);
   RooRealVar sigmalambda("sigmalambda", "width of Gaussian", 0.01, 0, 0.2);
   RooGaussian missing_massModel("missing_massModel", "Missing Mass Signal Model", *missing_mass, mlambda, sigmalambda); // Model for kaon mass signal
   mlambda.setConstant();



   // make the combined Z model
   RooProdPdf signal_Model("signal_Model", "2-d model for Z", RooArgSet(mkaonModel, missing_massModel)); // Signal 1 v 2

   // --------------------------------------
   // make QCD model

   std::cout << "make qcd model" << std::endl;
   // mass model for QCD.
   // the exponential is of the form exp(c*x).  If we want
   // the mass to decay an e-fold every R GeV, we use
   // c = -1/R.
   // We can leave this parameter free during the fit.
   RooRealVar mkaonbackPar1("mkaonbackPar1", "Polynomial 1st parameter", 130000);
   RooRealVar mkaonbackPar2("mkaonbackPar2", "Polynomial 2nd parameter", -1000000);
   RooRealVar mkaonbackPar3("mkaonbackPar3", "Polynomial 3rd parameter", 3500000);
   RooRealVar mkaonbackPar4("mkaonbackPar4", "Polynomial 4th parameter", -4750000);
   RooRealVar mkaonbackPar5("mkaonbackPar5", "Polynomial 5th parameter", 2400000);
   RooPolynomial mkaonbackModel("mkaonbackModel", "Missing Mass model", mkaon, RooArgList(mkaonbackPar1, mkaonbackPar2, mkaonbackPar3, mkaonbackPar4, mkaonbackPar5),4); // signal for missing_mass

   // missing_mass model for QCD.  Only used to generate toy MC
   // the exponential is of the form exp(c*x).  If we wantmass_kp
   // the missing_mass to decay an e-fold every R GeV, we use
   // c = -1/R.
   RooRealVar missing_massbackPar1("missing_massbackPar1", "Polynomial 1st parameter", 1700);
   RooRealVar missing_massbackPar2("missing_massbackPar2", "Polynomial 2nd parameter", -4000);
   RooRealVar missing_massbackPar3("missing_massbackPar3", "Polynomial 3rd parameter", 2900);
   RooPolynomial missing_massbackModel("missing_massbackModel", "missing mass back model", *missing_mass, RooArgList(missing_massbackPar1, missing_massbackPar2, missing_massbackPar3),2); // backgroud for missing_mass

   // make the 2-d model
   RooProdPdf background_Model("background_Model", "2-d model for QCD", RooArgSet(mkaonbackModel, missing_massbackModel));  // Backgrond 1 v 2

   // --------------------------------------
   // combined model

   // These variables represent the number of Z or QCD events
   // They will be fitted.
   RooRealVar zYield("zYield", "fitted yield for Z", 500000, 0., 1500000);
   RooRealVar qcdYield("qcdYield", "fitted yield for QCD", 300000, 0., 1500000);

   // now make the combined model
   std::cout << "make full model" << std::endl;
   RooAddPdf model("model", "z+qcd background models", RooArgList(signal_Model, background_Model), RooArgList(zYield, qcdYield));

   // interesting for debugging and visualizing the model
   model.graphVizTree("fullModel.dot");

   std::cout << "import model" << std::endl;

   ws->import(model);
}

//____________________________________
void AddData(RooWorkspace *ws)
{

  // Grab data file
  TFile *f1 = new TFile("/media/mn688/Elements1/PhD/Trees/Dibaryon/RGA/Strangeness_1/RGA_Fall2018_Inbending_skim4_Exclusive_Tree_Out_271021_01.root");
  TTree *t2 = (TTree*)f1->Get("t2");

  // Add a toy dataset
   // how many events do we want?
   // Int_t nEvents = 1000;

   // get what we need out of the workspace to make toy data
   RooAbsPdf *model = ws->pdf("model");
   RooRealVar *mass_kp = ws->var("mass_kp");
   RooRealVar *missing_mass = ws->var("missing_mass");
   // RooArgSet observables = new RooArgSet(mass_kp, missing_mass);
   // observables->Add(ws->var("mass_kp"));
   // observables->Add(ws->var("missing_mass"));
   // ws->defineSet("observables", *observables);
   // RooRealVar* x = new RooRealVar("x","x",0.2,0.8);
   // make the toy data
   std::cout << "make data set and import to workspace" << std::endl;
   // t2->SetBranchAddress("mass_kp",&mass_kp);
   // t2->SetBranchAddress("missing_mass",&missing_mass);
   RooDataSet data("data","data",RooArgSet(*mass_kp,*missing_mass),Import(*t2));
   // RooDataSet *data = model->generate(RooArgSet(*mass_kp, *missing_mass), nEvents);

   // import data into workspace
   ws->import(data, Rename("data"));
}

//____________________________________
void DoSPlot(RooWorkspace *ws)
{
   std::cout << "Calculate sWeights" << std::endl;

   // get what we need out of the workspace to do the fit
   RooAbsPdf *model = ws->pdf("model");
   RooRealVar *zYield = ws->var("zYield");
   RooRealVar *qcdYield = ws->var("qcdYield");
   RooDataSet *data = (RooDataSet *)ws->data("data");

   // fit the model to the data.
   model->fitTo(*data, Extended());

   // The sPlot technique requires that we fix the parameters
   // of the model that are not yields after doing the fit.
   //
   // This *could* be done with the lines below, however this is taken care of
   // by the RooStats::SPlot constructor (or more precisely the AddSWeight
   // method).
   //
   // RooRealVar* sigmaZ = ws->var("sigmaZ");
   // RooRealVar* qcdMassDecayConst = ws->var("qcdMassDecayConst");
   // sigmaZ->setConstant();
   // qcdMassDecayConst->setConstant();

   RooMsgService::instance().setSilentMode(true);

   std::cout << "\n\n------------------------------------------\nThe dataset before creating sWeights:\n";
   data->Print();

   RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

   // Now we use the SPlot class to add SWeights to our data set
   // based on our model and our yield variables
   RooStats::SPlot sData("sData", "An SPlot", *data, model, RooArgList(*zYield, *qcdYield));

   std::cout << "\n\nThe dataset after creating sWeights:\n";
   data->Print();

   // Check that our weights have the desired properties

   std::cout << "\n\n------------------------------------------\n\nCheck SWeights:" << std::endl;

   std::cout << std::endl
             << "Yield of Z is\t" << zYield->getVal() << ".  From sWeights it is "
             << sData.GetYieldFromSWeight("zYield") << std::endl;

   std::cout << "Yield of QCD is\t" << qcdYield->getVal() << ".  From sWeights it is "
             << sData.GetYieldFromSWeight("qcdYield") << std::endl
             << std::endl;

   for (Int_t i = 0; i < 10; i++) {
      std::cout << "z Weight for event " << i << std::right << std::setw(12) << sData.GetSWeight(i, "zYield") << "  qcd Weight"
                << std::setw(12) << sData.GetSWeight(i, "qcdYield") << "  Total Weight" << std::setw(12) << sData.GetSumOfEventSWeight(i)
                << std::endl;
   }

   std::cout << std::endl;

   // import this new dataset with sWeights
   std::cout << "import new dataset with sWeights" << std::endl;
   ws->import(*data, Rename("dataWithSWeights"));

   RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
}

void MakePlots(RooWorkspace *ws)
{

   // Here we make plots of the discriminating variable (mass_kp) after the fit
   // and of the control variable (missing_mass) after unfolding with sPlot.
   std::cout << "make plots" << std::endl;

   // make our canvas
   TCanvas *cdata = new TCanvas("sPlot", "sPlot demo", 400, 600);
   cdata->Divide(1, 3);

   // get what we need out of the workspace
   RooAbsPdf *model = ws->pdf("model");
   RooAbsPdf *signal_Model = ws->pdf("signal_Model");
   RooAbsPdf *background_Model = ws->pdf("background_Model");

   RooRealVar *missing_mass = ws->var("missing_mass");
   RooRealVar *mass_kp = ws->var("mass_kp");

   // note, we get the dataset with sWeights
   RooDataSet *data = (RooDataSet *)ws->data("dataWithSWeights");

   // this shouldn't be necessary, need to fix something with workspace
   // do this to set parameters back to their fitted values.
//   model->fitTo(*data, Extended());

   // plot mass_kp for data with full model and individual components overlaid
   //  TCanvas* cdata = new TCanvas();
   cdata->cd(1);
   RooPlot *frame = mass_kp->frame();
   data->plotOn(frame);
   model->plotOn(frame, Name("FullModel"));
   model->plotOn(frame, Components(*signal_Model), LineStyle(kDashed), LineColor(kRed), Name("Signal_Model"));
   model->plotOn(frame, Components(*background_Model), LineStyle(kDashed), LineColor(kGreen), Name("Background_Model"));

   TLegend leg(0.11, 0.5, 0.5, 0.8);
   leg.AddEntry(frame->findObject("FullModel"), "Full model", "L");
   leg.AddEntry(frame->findObject("Signal_Model"), "Z model", "L");
   leg.AddEntry(frame->findObject("Background_Model"), "QCD model", "L");
   leg.SetBorderSize(0);
   leg.SetFillStyle(0);

   frame->SetTitle("Fit of model to discriminating variable");
   frame->Draw();
   leg.DrawClone();

   // Now use the sWeights to show missing_mass distribution for Z and QCD.
   // The SPlot class can make this easier, but here we demonstrate in more
   // detail how the sWeights are used.  The SPlot class should make this
   // very easy and needs some more development.

   // Plot missing_mass for Z component.
   // Do this by plotting all events weighted by the sWeight for the Z component.
   // The SPlot class adds a new variable that has the name of the corresponding
   // yield + "_sw".
   cdata->cd(2);

   // create weighted data set
   RooDataSet *dataw_z = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "zYield_sw");

   RooPlot *frame2 = missing_mass->frame();
   // Since the data are weighted, we use SumW2 to compute the errors.
   dataw_z->plotOn(frame2, DataError(RooAbsData::SumW2));

   frame2->SetTitle("missing_mass distribution with s weights to project out Z");
   frame2->Draw();

   // Plot missing_mass for QCD component.
   // Eg. plot all events weighted by the sWeight for the QCD component.
   // The SPlot class adds a new variable that has the name of the corresponding
   // yield + "_sw".
   cdata->cd(3);
   RooDataSet *dataw_qcd = new RooDataSet(data->GetName(), data->GetTitle(), data, *data->get(), 0, "qcdYield_sw");
   RooPlot *frame3 = missing_mass->frame();
   dataw_qcd->plotOn(frame3, DataError(RooAbsData::SumW2));

   frame3->SetTitle("missing_mass distribution with s weights to project out QCD");
   frame3->Draw();

   //  cdata->SaveAs("SPlot.gif");
}
