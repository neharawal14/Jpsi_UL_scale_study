// this file is just to analyze those produced histograms by the DY smaples
//
#include<iostream>
using namespace std;
#include <TROOT.h>  
#include <TChain.h>
#include <TFile.h> 
// Header file for the classes stored in the TTree if any.                                                                                                                          
#include "string" 
#include "vector" 
#include "RooRealVar.h"   
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h" 
#include "RooHistPdf.h" 
#include "RooPolynomial.h"
#include "RooAbsArg.h"   
#include "RooPlot.h"    
//#include "TRatioPlot.h" 
#include "RooAddPdf.h" 
#include "RooFitResult.h"
#include "TAxis.h"      
#include "TH1.h"                                                                                                                                                                    
using namespace RooFit ;

// base class having all the functions to analyze the histograms
// The functions defined in base class are called by derived class and used in that way
class class_reading{
	public:
   TFile *file = TFile::Open("Jpsi_ntuple_smear1percent_15KeV_1MeV_test.root"); 
   TString saving_path = "/cmsuf/data/store/user/t2/users/neha.rawal/UL/UL_samples/analyze_bins/pt_analysis/analysis_folders/1_percent_files/results/test_results/";    
   float pt_list_mean[9] = {1.25,2.5,3.5,4.5,5.5,6.5,7.5,9.5,15.5};                                                                                                                
   float pt_list[10] = {0.5,2,3,4,5,6,7,8,11,20};                                                                                                                                  
   TString bin_number_pt[9] = {"1st","2nd","3rd","4th","5th","6th","7th", "8th", "9th"};     

	 // hisotgram declaration to read the histograms from the root file
    TH1F * histogram_gen_Zmass; 
    TH1F * histogram_gen_Zmass_reconstruct; 
    TH1F * histogram_gen_positive_mu ;
    TH1F * histogram_id_positive_mu; 
    TH1F * histogram_gen_negative_mu; 
    TH1F * histogram_id_negative_mu; 
		TH1F * histogram_eta_gen_positive_mu; 
		TH1F * histogram_eta_gen_negative_mu; 
		TH1F * histogram_phi_gen_positive_mu; 
		TH1F * histogram_phi_gen_negative_mu; 
    TH1F * histogram_gen_positive_mu_smear; 
    TH1F * histogram_id_positive_mu_smear; 
    TH1F * histogram_gen_negative_mu_smear; 
    TH1F * histogram_id_negative_mu_smear; 
    TH1F * histogram_reco_Zmass; 
    TH1F * histogram_pt_gen_Z; 
    TH1F * histogram_eta_gen_Z; 
    TH1F * histogram_phi_gen_Z; 

     // this histograms are defined in a loop since its an array  
     TH1F *histogram_Zmass_positive_pt[9];
     TH1F *histogram_Zmass_negative_pt[9];
     TH1F *histogram_Zmass_positive_pt_smear[9];
     TH1F *histogram_Zmass_negative_pt_smear[9];
   
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
    
    TH1F *histogram_pT_positive_bin[9] ; 
    TH1F *histogram_pT_negative_bin[9] ; 
    TH1F *histogram_pT_positive_bin_smear[9] ; 
    TH1F *histogram_pT_negative_bin_smear[9] ; 
    TH1F *histogram_pT_negative_corresponding_positive_bin[9] ; 
    TH1F *histogram_pT_positive_corresponding_negative_bin[9] ; 
    TH1F *histogram_pT_negative_corresponding_positive_bin_smear[9] ; 
    TH1F *histogram_pT_positive_corresponding_negative_bin_smear[9] ; 

		// delcared variables to analyze the mean and sigma of the distribution
   float mean_positive_pt[9], mean_negative_pt[9], sigma_positive_pt[9], sigma_negative_pt[9];
   float mean_positive_pt_error[9], mean_negative_pt_error[9], sigma_positive_pt_error[9], sigma_negative_pt_error[9];
   float mean_positive_reco_pt[9], mean_negative_reco_pt[9], sigma_positive_reco_pt[9],sigma_negative_reco_pt[9];
   float mean_positive_reco_pt_error[9], mean_negative_reco_pt_error[9], sigma_positive_reco_pt_error[9],sigma_negative_reco_pt_error[9];
   float diff_Zmass_positive_gen_reco[9];
   float diff_Zmass_negative_gen_reco[9];
	 float diff_Zmass_positive_gen_reco_error[9];
   float diff_Zmass_negative_gen_reco_error[9];
 
	 // constructor to open the root file
	 void initializing();

   // functions used to fit the histograms with BW, DSCB and gaussian distributions	 
   std::pair<float,float> plotting_fitting_mass_BW(TH1F* hist_fit, TString saving_name, TString title_name); 
   std::pair<float,float> plotting_fitting_mass_BW_smearing(TH1F* hist_fit, TString saving_name, TString title_name);

	 // for gaussian fit
   std::pair<float,float> plotting_fitting_mass_gauss(TH1F* hist_fit, TString saving_name, TString title_name); 
   std::pair<float,float> plotting_fitting_mass_gauss_smearing(TH1F* hist_fit, TString saving_name, TString title_name);  

	 // for DSCB fit before smearing
	 std::pair<float,float> plotting_fitting_mass_DSCB(TH1F* hist_fit, TString saving_name, TString title_name);
   std::pair<float,float> plotting_fitting_mass_DSCB_bin(TH1F* hist_fit, TString saving_name, TString title_name); 

  	// for DSCB fit after smearing 
    //  std::pair<float,float> plotting_fitting_mass_DSCB_smearing(TH1F* hist_fit, TString saving_name, TString title_name);  
    //  std::pair<float,float> plotting_fitting_mass_DSCB_bin_smearing(TH1F* hist_fit, TString saving_name, TString title_name);  

	 void plotting_mass_distribution(TH1F *hist_inclusive, TH1F *hist_bin, float range_min, float range_max, TString title_name, TString saving_name, TString title_name_axis);        
  	// to plot histogram for pT distribution and mass distribution individually 
   void plotting_hist(TH1F* hist_draw, TString title, TString saving_name, TString title_name_axis);
	 	// to count total number of entries in each pT bin and also above 200 GeV and print them in a text file myfile_n_entries                                                          
  	// these 3 functions are for sanity check
    // void plot_entries();   
    // void total_entries();   
    // void evaluating_mean_pT(); 
    // plotting mean and sigma distribution 		

    // to save all the histograms in pdf plots	 
    void saving_histogram_pT();

		// graphs to plot mean and sigma, and difference due to smearing, and  also print them in  text files
	 	void graph_mean_pT_combine_smear();
    void graph_mean_pT_combine();    
    void graph_diff_Zmass_gen_reco(); 
};

// declaring other class and its variables - the one which is actually used -=> inherited class from class reading
class derived_class_reading : public class_reading{
	public:
  TString pt_list_symbol[10] = {"0.5", "2.0", "3.0", "4","5", "6.0", "7.0", "8", "11", "20"};                                                                                     
	std::pair<float,float> mean_sigma_Zmass; 
  std::pair<float,float>mean_sigma_Zmass_reconstruct;
  std::pair<float,float>mean_sigma_Zmass_smear; 

//  std::pair<float,float>mean_sigma_Zmass_gauss; 
  std::pair<float,float> mean_sigma_Zmass_gauss_smear; 
  

  std::pair<float,float>mean_sigma_Zmass_DSCB; 
  // std::pair<float,float> mean_sigma_Zmass_DSCB_smear; 
  
	std::pair<float,float> mean_sigma_positive_pt[9]; 
  std::pair<float,float> mean_sigma_positive_reco_pt[9]; 
  std::pair<float,float> mean_sigma_negative_pt[9]; 
  std::pair<float,float> mean_sigma_negative_reco_pt[9]; 

  int n_positive_bin_smear[9], n_negative_bin_smear[9], n_positive_bin[9], n_negative_bin[9]; 	
	void plot_histograms();	
	void fitting_histograms();
  void plotting_inclusive_bin();
  void mean_sigma_calculation();	
  void saving_text_file();
};

  // declaring base class functions
	//
	// to plot mean Jpsi after smearing
   void class_reading :: graph_mean_pT_combine_smear(){                                                                                                                                     

   // this is to store the mean and error on mean in a text file names mean_check_smear.txt , so can be used further for comparison
  ofstream check_file_mean_smear;
  check_file_mean_smear.open("mean_check_smear.txt",std::ios_base::app);	

  check_file_mean_smear<<"Mean of positive mu smear"<<std::endl;
  check_file_mean_smear<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean_smear<<mean_positive_reco_pt[i]<<", ";
	}	
  check_file_mean_smear<<"} "<<std::endl;

	check_file_mean_smear<<"Mean of negative mu smear"<<", ";
  check_file_mean_smear<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean_smear<<mean_negative_reco_pt[i]<<", ";
	}
  check_file_mean_smear<<"} "<<std::endl;
 
 	check_file_mean_smear<<"Error on mean of positive mu smear"<<std::endl;
  check_file_mean_smear<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean_smear<<mean_positive_reco_pt_error[i]<<", ";
	}	
  check_file_mean_smear<<"} "<<std::endl;
	
	check_file_mean_smear<<"Error on mean of negative mu smear"<<std::endl;
  check_file_mean_smear<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean_smear<<mean_negative_reco_pt_error[i]<<", ";
	}
  check_file_mean_smear<<"} "<<std::endl;

   TCanvas *graph_canvas_negative = new TCanvas("graph_canvas_negative","mass(#mu^{+}#mu^{-})",900,600);
   graph_canvas_negative->cd(); 
   TMultiGraph *mg1 = new TMultiGraph(); 
   mg1->SetTitle("Mean Zmass after smearing");
    
   TGraphErrors *gr_positive_reco = new TGraphErrors(9, pt_list_mean,mean_positive_reco_pt, 0,mean_positive_reco_pt_error);
   gr_positive_reco->SetMarkerColor(4);
   gr_positive_reco->SetMarkerSize(1.0);
   gr_positive_reco->SetMarkerStyle(21);
   gr_positive_reco->SetName("gr_positive_reco"); 
  
   TGraphErrors *gr_negative_reco = new TGraphErrors(9, pt_list_mean,mean_negative_reco_pt, 0,mean_negative_reco_pt_error);
   gr_negative_reco->SetMarkerColor(kRed+2);
   gr_negative_reco->SetMarkerSize(1.0);
   gr_negative_reco->SetMarkerStyle(21);
   gr_negative_reco->SetName("gr_negative_reco"); 
   
	 mg1->Add(gr_positive_reco);
   mg1->Add(gr_negative_reco);
   mg1->Draw("AP*");
   mg1->GetYaxis()->SetTitle("mass_{(#mu^{+}#mu^{-})} (GeV)");
   mg1->GetXaxis()->SetTitle("pT (GeV)");
   mg1->GetYaxis()->SetRangeUser(3.0955,3.0985);
  
    TLegend* leg1 = new TLegend(0.7, 0.15, 0.88, 0.26);
    leg1->AddEntry("gr_positive_reco","#mu^{+}", "EP");
    leg1->AddEntry("gr_negative_reco","#mu^{-}", "EP");
    leg1->Draw();

    graph_canvas_negative->SaveAs(saving_path+"mean_pT_mu_after_smearing_combine.pdf"); 
  
   }
// this is to plot mean Jpsi mass before smearing

void class_reading :: graph_mean_pT_combine (){                                                                                                                                 
   // this is to store the mean and error on mean in a text file names mean_check.txt , so can be used further for comparison
   ofstream check_file_mean;
   check_file_mean.open("mean_check.txt",std::ios_base::app);	

	 check_file_mean<<"Mean of positive mu smear"<<std::endl;
  check_file_mean<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean<<mean_positive_pt[i]<<", ";
	}	
  check_file_mean<<"} "<<std::endl;

	check_file_mean<<"Mean of negative mu smear"<<std::endl;
  check_file_mean<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean<<mean_negative_pt[i]<<", ";
}
  check_file_mean<<"} "<<std::endl;

	check_file_mean<<"Error on mean of positive mu smear"<<std::endl;
  check_file_mean<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean<<mean_positive_pt_error[i]<<", ";
	}	
  check_file_mean<<"} "<<std::endl;
	
	check_file_mean<<"Error on mean of negative mu smear"<<std::endl;
  check_file_mean<<"{ ";
	for(int i=0; i<9; i++){
  check_file_mean<<mean_negative_pt_error[i]<<", ";
	}
  check_file_mean<<"} "<<std::endl;

 
	 TCanvas *graph_canvas = new TCanvas("graph_canvas","mass(#mu^{+}#mu^{-})",900,600);
 graph_canvas->cd(); 
 TMultiGraph *mg = new TMultiGraph();
 mg->SetTitle("mean Zmass before smearing");
 TGraphErrors *gr_positive = new TGraphErrors(9, pt_list_mean,mean_positive_pt, 0,mean_positive_pt_error);
 gr_positive->SetLineColor(2);
 gr_positive->SetMarkerColor(4); 
 gr_positive->SetMarkerSize(1.0);
 gr_positive->SetMarkerStyle(21); 
 gr_positive->SetName("gr_positive");
 
 TGraphErrors *gr_negative = new TGraphErrors(9, pt_list_mean,mean_negative_pt, 0, mean_negative_pt_error);
 gr_negative->SetMarkerColor(kRed+2);
 gr_negative->SetMarkerSize(1.0);
 gr_negative->SetMarkerStyle(21);
 gr_negative->SetName("gr_negative"); 
                                     
 mg->Add(gr_positive);              
 mg->Add(gr_negative);             
 mg->Draw("AP*");                 
 mg->GetYaxis()->SetTitle("mass_{(#mu^{+}#mu^{-})} (GeV)");  
 mg->GetXaxis()->SetTitle("pT (GeV)"); 
 mg->GetYaxis()->SetRangeUser(3.0967,3.0971); 
                                                                                                                                                                                    
 TLegend* leg = new TLegend(0.7, 0.15, 0.88, 0.26);
 leg->AddEntry("gr_positive","#mu^{+}", "EP");    
 leg->AddEntry("gr_negative", "#mu^{-}", "EP");  
 leg->Draw();                                   
 graph_canvas->SaveAs(saving_path+"mean_pT_mu_before_smearing_combine.pdf"); 
 }

void class_reading :: graph_diff_Zmass_gen_reco(){
   TGraphErrors *gr_diff_gen_reco_positive = new TGraphErrors(9, pt_list_mean,diff_Zmass_positive_gen_reco, 0, diff_Zmass_positive_gen_reco_error);
   TGraphErrors *gr_diff_gen_reco_negative = new TGraphErrors(9, pt_list_mean,diff_Zmass_negative_gen_reco, 0, diff_Zmass_negative_gen_reco_error);
   gr_diff_gen_reco_positive->SetMarkerColor(4);
   gr_diff_gen_reco_positive->SetMarkerSize(1.0);
   gr_diff_gen_reco_positive->SetMarkerStyle(21);

   gr_diff_gen_reco_negative->SetMarkerColor(kRed+2);
   gr_diff_gen_reco_negative->SetMarkerSize(1.0);
   gr_diff_gen_reco_negative->SetMarkerStyle(21);
 
   TMultiGraph *mg = new TMultiGraph(); 
   mg->SetTitle("difference Z mass(mean) before and after smearing");
  
   TCanvas *graph_canvas = new TCanvas("graph_canvas","diff Gen and Reco Z mass",900,600);
   graph_canvas->cd(); 
   gr_diff_gen_reco_positive->SetName("positive");
   gr_diff_gen_reco_negative->SetName("negative");
   mg->Add(gr_diff_gen_reco_positive);
   mg->Add(gr_diff_gen_reco_negative);
   mg->Draw("AP*");
   mg->GetYaxis()->SetTitle("");
   mg->GetXaxis()->SetTitle("pT (GeV)");
   mg->GetYaxis()->SetRangeUser(-0.0004,0.0004);
   TLegend* leg = new TLegend(0.7, 0.15, 0.88, 0.26);
   leg->AddEntry("positive","#mu^{+}", "LP");
   leg->AddEntry("negative","#mu^{-}", "LP");
   leg->Draw();

   graph_canvas->SaveAs(saving_path+"diff_Zmas_gen_combine.pdf"); 

	 ofstream myfile_diff_Zmass;
   // to plot the difference in Z mass with the above smearing percentage in  a text file
   myfile_diff_Zmass.open("diff_Zmass.txt",std::ios_base::app);
   myfile_diff_Zmass<<" difference in Z mass"<<std::endl;
   myfile_diff_Zmass<<" smearing   : 1%"<<std::endl;
   myfile_diff_Zmass<<" for #mu^{+}"<<std::endl;
   
	 myfile_diff_Zmass<<"{  ";
   for(int i=0 ;i<9;i++){
   myfile_diff_Zmass<<diff_Zmass_positive_gen_reco[i]<<", ";
   } 
   myfile_diff_Zmass<<"}  "<<std::endl;

   myfile_diff_Zmass<<" for #mu^{-}"<<std::endl;
	 myfile_diff_Zmass<<"{  ";
   for(int i=0 ;i<9;i++){
   myfile_diff_Zmass<<diff_Zmass_negative_gen_reco[i]<<", ";
  }
   myfile_diff_Zmass<<"}  "<<std::endl;

  // putting the error on difference in the text file
  myfile_diff_Zmass<<" error on the difference in Z mass"<<std::endl;
   myfile_diff_Zmass<<" for #mu^{+}"<<std::endl;
	 myfile_diff_Zmass<<"{  ";
	 for(int i=0 ;i<9;i++){
   myfile_diff_Zmass<<diff_Zmass_positive_gen_reco_error[i]<<", ";
   } 
   myfile_diff_Zmass<<"}  "<<std::endl;

   myfile_diff_Zmass<<" for #mu^{-}"<<std::endl;
	 myfile_diff_Zmass<<"{  ";
   for(int i=0 ;i<9;i++){
   myfile_diff_Zmass<<diff_Zmass_negative_gen_reco_error[i]<<", ";
  }
   myfile_diff_Zmass<<"}  "<<std::endl;



 
 }

   void class_reading :: plotting_hist(TH1F * hist_draw, TString title, TString saving_name, TString title_name_axis){
    TCanvas * c = new TCanvas("c",title,900,600);
    c->cd();
    hist_draw->GetXaxis()->SetTitle(title_name_axis+"");
    hist_draw->Draw();
    c->SaveAs(saving_path+saving_name+".pdf");
   }
 
    void class_reading :: plotting_mass_distribution(TH1F *hist_inclusive, TH1F *hist_bin,float range_min, float range_max,  TString title_name, TString saving_name, TString title_name_axis){

     gStyle->SetOptStat(0);
     TH1F *hist_inclusive_clone = (TH1F*)hist_inclusive->Clone();
     TH1F *hist_bin_clone = (TH1F*)hist_bin->Clone();

 		 hist_inclusive_clone->Scale(1/hist_inclusive_clone->Integral());
     hist_bin_clone->Scale(1/hist_bin_clone->Integral());

    TCanvas * c = new TCanvas("c",title_name,900,600);
    c->cd();
    hist_bin_clone->GetXaxis()->SetTitle(title_name_axis+"");
//    hist_bin_clone->GetYaxis()->SetTitle("N/20(KeV)");
    hist_bin_clone->TAttLine::SetLineColor(kRed+2);
    hist_bin_clone->SetName("hist_bin_clone");
    hist_inclusive_clone->SetName("hist_inclusive_clone");
    hist_bin_clone->SetTitle(title_name);
    hist_bin_clone->GetXaxis()->SetRangeUser(range_min, range_max);
    gStyle->SetOptStat(0);
    hist_bin_clone->Draw("E");
    hist_inclusive_clone->TAttLine::SetLineColor(kBlue);
    hist_inclusive_clone->Draw("SAMES");
 
    gStyle->SetOptStat(0);
     TLegend* leg1 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg1->SetFillColor(kWhite);
     leg1->SetLineColor(kBlack);
     leg1->AddEntry("hist_inclusive_clone","inclusive", "LP");
     leg1->AddEntry("hist_bin_clone","bin","LP");
     leg1->Draw("same");       

     c->SaveAs(saving_path+saving_name+"normalized.pdf");
		 c->Close();
    }

   std::pair<float,float> class_reading :: plotting_fitting_mass_BW(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into Breit Wigner
     RooRealVar mass_var("mass_var","mass_var",3.0968,3.09704);
     RooDataHist histo("histo","mass dataset",mass_var,hist_fit);
     RooRealVar mean_mass("mean_mass","mean of Z mass",3.097,3.0967,3.0971);
     RooRealVar width("width","width of Z mass",0.0001,0.000001,0.09);
   
     RooPlot *xframe=mass_var.frame(Title(title_name));
     histo.plotOn(xframe);
//     histo.statOn(xframe);
   
     RooBreitWigner BW("BW","Breit Wigner fit",mass_var, mean_mass,width);
     BW.fitTo(histo,Range(3.0968,3.09704));
    
     BW.plotOn(xframe,RooFit::LineColor(kRed+2),Name("BW_sig"));
    
     BW.paramOn(xframe,RooFit::Layout(0.6,0.9,0.7));
     
     TCanvas *tmp = new TCanvas("tmp","Gen J/#Psi mass", 900,600);
     tmp->cd();
     gPad->SetLeftMargin(0.15);
     xframe->getAttText()->SetTextSize(0.025);
     xframe->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
   
     xframe->GetYaxis()->SetTitle("N/15 (KeV)");
     xframe->Draw();
   
     TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("BW_sig","BW fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       


  
     gStyle->SetOptStat();
     tmp->SaveAs(saving_path+ saving_name+".pdf");
     tmp->Close(); 
     std::cout<<" Z mass mean : "<<mean_mass.getVal()<<" width : "<<width.getVal()<<std::endl;
     return std::make_pair(mean_mass.getVal(),mean_mass.getError());   
     }
 
  std::pair<float,float> class_reading :: plotting_fitting_mass_BW_smearing(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into Breit Wigner
     RooRealVar mass_var("mass_var","mass_var",3.02,3.18);
     RooDataHist histo("histo","mass dataset",mass_var,hist_fit);
     RooRealVar mean_mass("mean_mass","mean of Z mass",3.097,3.096,3.098);
     RooRealVar width("width","width of Z mass",0.0001,0.000001,0.09);
   
     RooPlot *xframe=mass_var.frame(Title(title_name));
     histo.plotOn(xframe);
//     histo.statOn(xframe);
   
     RooBreitWigner BW("BW","Breit Wigner fit",mass_var, mean_mass,width);
     BW.fitTo(histo,Range(3.02,3.18));
    
     BW.plotOn(xframe,RooFit::LineColor(kRed+2),Name("BW_sig"));
    
     BW.paramOn(xframe,RooFit::Layout(0.6,0.9,0.7));
     
     TCanvas *tmp = new TCanvas("tmp","Gen J/#Psi mass", 900,600);
     tmp->cd();
     gPad->SetLeftMargin(0.15);
     xframe->getAttText()->SetTextSize(0.025);
     xframe->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV)");
   
     xframe->GetYaxis()->SetTitle("N/1 (MeV)");
     xframe->Draw();
   
     TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("BW_sig","BW fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       


  
     gStyle->SetOptStat();
     tmp->SaveAs(saving_path+ saving_name+".pdf");
     tmp->Close(); 
     std::cout<<" Z mass mean : "<<mean_mass.getVal()<<" width : "<<width.getVal()<<std::endl;
     return std::make_pair(mean_mass.getVal(),mean_mass.getError());   
     }

/*   std::pair<float,float> class_reading :: plotting_fitting_mass_gauss(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into gauss
    RooRealVar m4mu("m4mu", "var", 3.0968,3.09704, ""); 
    RooDataHist histo("histo","dataset with var",m4mu,hist_fit);
    RooRealVar Mean("Mean", "Mean",3.097, 3.0968, 3.09704);
    RooRealVar Sigma("#sigma", "#sigma", 0.00001, 0.000001,0.09);//sigma[decay]);
    RooGaussian Gauss("Gauss", "Gausian PDF", m4mu, Mean, Sigma);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
   
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);
//    histo.statOn(xframe);
   
    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    Gauss.fitTo(histo, Range(3.09684,3.097));
    Gauss.plotOn(xframe, RooFit::LineColor(color),Name("Gauss_sig"));
    Gauss.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");

    xframe->GetYaxis()->SetTitle("N/15 (KeV)");

    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
   TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("Gauss_sig","Gauss fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

    gStyle->SetOptStat();
    c_MC->SaveAs((saving_path+ saving_name + "Gauss.pdf"));// + ".pdf");
   
	 c_MC->Close();	
    std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
    return std::make_pair(Mean.getVal(),Mean.getError());   
}
*/
std::pair<float,float> class_reading :: plotting_fitting_mass_gauss_smearing(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into gauss
    RooRealVar m4mu("m4mu", "var", 3.02,3.18, ""); 
    RooDataHist histo("histo","dataset with var",m4mu,hist_fit);
    RooRealVar Mean("Mean", "Mean",3.097, 3.08, 3.14);
    RooRealVar Sigma("#sigma", "#sigma", 0.02, 0.000009,0.09);//sigma[decay]);
    RooGaussian Gauss("Gauss", "Gaussian PDF", m4mu, Mean, Sigma);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
   
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);
//    histo.statOn(xframe);
   
    Int_t color = kRed+2;
    Double_t size_text = 0.020;
//    Gauss.fitTo(histo, Range(3.062,3.13));
    Gauss.fitTo(histo, Range(3.05,3.145));
    Gauss.plotOn(xframe, RooFit::LineColor(color),Name("Gauss_sig"));
    Gauss.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");

    xframe->GetYaxis()->SetTitle("N/1 (MeV)");

    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
   TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("Gauss_sig","Gauss fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

    gStyle->SetOptStat();
    c_MC->SaveAs((saving_path+ saving_name + "Gauss.pdf"));// + ".pdf");
   
	 c_MC->Close();	
    std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
    return std::make_pair(Mean.getVal(),Mean.getError());   
}



   std::pair<float,float> class_reading :: plotting_fitting_mass_DSCB(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into DSCB
    RooRealVar m4mu("m4mu", "var", 3.0968,3.09704, ""); 
    RooDataHist histo("histo","dataset with var",m4mu,hist_fit);
    RooRealVar Mean("Mean", "Mean",3.097, 3.0968, 3.09704);
    RooRealVar Sigma("#sigma", "#sigma", 0.00003, 0.000009,0.09);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 0.8, 0.0001, 100);//alphaL[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 80, 0.1, 500);//expL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1.3, 0.0001, 100);//alphaR[decay]);
    RooRealVar ExpR("n_{R}", "n_{R}", 0.89, 0.1, 500);//expR[decay]);
    RooMyPDF_DSCB DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
   
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);
//    histo.statOn(xframe);
   
    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(3.0968,3.09704));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("DSCB_sig"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");

    xframe->GetYaxis()->SetTitle("N/15 (KeV)");

    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
   TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("DSCB_sig","DSCB fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

    gStyle->SetOptStat();
    c_MC->SaveAs((saving_path+ saving_name + "DSCB.pdf"));// + ".pdf");
   
	 c_MC->Close();	
    std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
    return std::make_pair(Mean.getVal(),Mean.getError());   
}
   std::pair<float,float> class_reading :: plotting_fitting_mass_DSCB_bin(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into DSCB
    RooRealVar m4mu("m4mu", "var", 3.0968,3.09704, ""); 
    RooDataHist histo("histo","dataset with var",m4mu,hist_fit);
    RooRealVar Mean("Mean", "Mean",3.097, 3.0968, 3.09704);
    RooRealVar Sigma("#sigma", "#sigma", 0.00002, 0.000001,0.09);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1, 0.0001, 200);//alphaL[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 13, 0.1, 100);//expL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2, 0.0001, 200);//alphaR[decay]);
    RooRealVar ExpR("n_{R}", "n_{R}", 5, 0.1, 100);//expR[decay]);
    RooMyPDF_DSCB DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
   
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);
//    histo.statOn(xframe);
   
    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(3.0968,3.09704));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("DSCB_sig"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");

    xframe->GetYaxis()->SetTitle("N/15 (KeV)");

    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
   TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("DSCB_sig","DSCB fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

    gStyle->SetOptStat();
    c_MC->SaveAs((saving_path+ saving_name + "DSCB.pdf"));// + ".pdf");
   
	 c_MC->Close();	
    std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
    return std::make_pair(Mean.getVal(),Mean.getError());   
}
/*
std::pair<float,float> class_reading :: plotting_fitting_mass_DSCB_smearing(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into DSCB
    RooRealVar m4mu("m4mu", "var", 3.02,3.18, ""); 
    RooDataHist histo("histo","dataset with var",m4mu,hist_fit);
    RooRealVar Mean("Mean", "Mean",3.097, 3.08, 3.14);
    RooRealVar Sigma("#sigma", "#sigma", 0.02, 0.000009,0.09);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.001, 50);//alphaL[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 80, 0.1, 500);//expL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.0, 0.001, 50);//alphaR[decay]);
    RooRealVar ExpR("n_{R}", "n_{R}", 80, 0.1, 500);//expR[decay]);
    RooMyPDF_DSCB DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
   
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);
//    histo.statOn(xframe);
   
    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(3.04,3.145));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("DSCB_sig"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");

    xframe->GetYaxis()->SetTitle("N/1 (MeV)");

    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
   TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("DSCB_sig","DSCB fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

    gStyle->SetOptStat();
    c_MC->SaveAs((saving_path+ saving_name + "DSCB.pdf"));// + ".pdf");
   
	 c_MC->Close();	
    std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
    return std::make_pair(Mean.getVal(),Mean.getError());   
}

std::pair<float,float> class_reading :: plotting_fitting_mass_DSCB_bin_smearing(TH1F * hist_fit, TString saving_name, TString title_name){
  // Fit the mass into DSCB
    RooRealVar m4mu("m4mu", "var", 3.02,3.18, ""); 
    RooDataHist histo("histo","dataset with var",m4mu,hist_fit);
    RooRealVar Mean("Mean", "Mean",3.097, 3.08, 3.14);
    RooRealVar Sigma("#sigma", "#sigma", 0.02, 0.000009,0.09);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 2.0, 0.1, 50);//alphaL[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 100, 1, 300);//expL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 2.0, 0.1, 50);//alphaR[decay]);
    RooRealVar ExpR("n_{R}", "n_{R}", 20, 1, 300);//expR[decay]);
    RooMyPDF_DSCB DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
   
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo.plotOn(xframe);
//    histo.statOn(xframe);
   
    Int_t color = kRed+2;
    Double_t size_text = 0.020;
    DSCB.fitTo(histo, Range(3.04,3.145));
    DSCB.plotOn(xframe, RooFit::LineColor(color),Name("DSCB_sig"));
    DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    c_MC->cd();
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(color);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");

    xframe->GetYaxis()->SetTitle("N/1 (MeV)");

    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
   TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("histo","Gen J/#Psi", "EP");
     leg2->AddEntry("DSCB_sig","DSCB fit","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     leg2->AddEntry("histo.sumEntries()",Form("Events= %.0f",histo.sumEntries()),"");
     leg2->Draw("same");       

    gStyle->SetOptStat();
    c_MC->SaveAs((saving_path+ saving_name + "DSCB.pdf"));// + ".pdf");
   
	 c_MC->Close();	
    std::cout<<" Z mass mean : "<<Mean.getVal()<<" width : "<<Sigma.getVal()<<std::endl;
    return std::make_pair(Mean.getVal(),Mean.getError());   
}



void class_reading :: plotting_fitting_mass_smearing_DSCB_single(TH1F * hist_fit_gen, TH1F *hist_fit_smear, TString title_name, TString saving_name){
  // Fit the mass into DSCB
 
      TH1F *hist_fit_gen_clone = (TH1F*)hist_fit_gen->Clone();
      TH1F *hist_fit_smear_clone = (TH1F*)hist_fit_smear->Clone();
 
    RooRealVar m4mu("m4mu", "var", 3.02, 3.18, ""); 
    RooDataHist histo_smear("histo_smear","dataset with var",m4mu,hist_fit_smear_clone);
    RooDataHist histo_gen("histo_gen","dataset with var",m4mu,hist_fit_gen_clone);

    RooRealVar Mean("Mean", "Mean",3.11, 3.0, 3.2);
    RooRealVar Sigma("#sigma", "#sigma", 0.0001, 0.00001,0.09);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 0.6, 0.01, 30);//alphaL[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 15, 0.05, 100);//expL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 0.6, 0.01, 30);//alphaR[decay]);
    RooRealVar ExpR("n_{R}", "n_{R}", 25, 0.05, 100);//expR[decay]);
    RooMyPDF_DSCB DSCB1("DSCB1", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    RooMyPDF_DSCB DSCB2("DSCB2", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
   
    TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 1000, 600);
    c_MC->SetFrameFillColor(0);
 
      c_MC->cd();
    RooPlot* xframe = m4mu.frame(RooFit::Title(title_name));
    histo_smear.plotOn(xframe, MarkerColor(kRed+2),Name("hist_smear"));
    // histo.statOn(xframe);
    Int_t color2 = kRed +2;
    Double_t size_text2 = 0.020;
    DSCB2.fitTo(histo_smear, Range(3.02,3.18));
    DSCB2.plotOn(xframe, RooFit::LineColor(color2),Name("DSCB_sig2"));
    DSCB2.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.50));
    xframe->getAttText()->SetTextSize(size_text2);
    xframe->getAttText()->SetTextColor(color2);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");
//    xframe->GetYaxis()->SetTitle("N/0.5 (MeV)");
    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
   
//    c_MC->cd();
  
     Int_t color1 = kBlue;
    Double_t size_text1 = 0.020;
    histo_gen.plotOn(xframe, MarkerColor(kBlue), Name("hist_gen"));
    DSCB1.fitTo(histo_gen, Range(3.0968, 3.09704));
    DSCB1.plotOn(xframe, RooFit::LineColor(color1),Name("DSCB_sig1"));
    DSCB1.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
    xframe->getAttText()->SetTextSize(size_text1);
    xframe->getAttText()->SetTextColor(color1);
    xframe->GetXaxis()->SetTitle("m(#mu^{+} #mu^{-}) (GeV)");
    // xframe->GetYaxis()->SetTitle("N/1 (MeV)");
    xframe->GetXaxis()->SetTitleOffset(1.4);
    xframe->chiSquare();
    xframe->Draw();
 
     TLegend* leg2 = new TLegend(0.7, 0.7, 0.88, 0.88);
     leg2->SetFillColor(kWhite);
     leg2->SetLineColor(kBlack);
     leg2->AddEntry("hist_gen"," Gen J/#Psi before smearing", "EP");
     leg2->AddEntry("DSCB_sig1","DSCB fit before smearing","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     // leg2->AddEntry("histo_gen.sumEntries()",Form("Events= %.0f",histo_gen.sumEntries()),"");
     leg2->Draw("same");       


     leg2->AddEntry("hist_smear"," Gen J/#Psi after smearing", "EP");
     leg2->AddEntry("DSCB_sig2","DSCB fit after smearing","LP");
     leg2->AddEntry("frame->chiSquare()",Form("#chi^{2}/ndf= %.2f",xframe->chiSquare()),"");
     // leg2->AddEntry("histo_smear.sumEntries()",Form("Events= %.0f",histo_smear.sumEntries()),"");
     leg2->Draw("same");      
     gStyle->SetOptStat();
		 
     c_MC->SaveAs((saving_path+ saving_name + "DSCB_before_after_smearing.pdf"));// + ".pdf");
		 c_MC->Close();
}
*/
void class_reading :: saving_histogram_pT(){
   TString string_pT_distribution_negative, string_pT_distribution_positive;  
   TString string_pT_distribution_negative_smear, string_pT_distribution_positive_smear;  
 
   // histograms and corresponding other pT distribution
  for(int i=0; i<9; i++){
  TCanvas *canvas_pT_positive = new TCanvas("canvas_pT_positive", "pT distribution in each bin #mu^{+}", 900,600);
  canvas_pT_positive->Divide(2,1);
  canvas_pT_positive->cd(1);
  histogram_pT_positive_bin[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_positive_bin[i]->GetYaxis()->SetRangeUser(0,300);
  histogram_pT_positive_bin[i]->Draw();
  string_pT_distribution_positive ="pT_positive_" + bin_number_pt[i] + "bin";
  
  canvas_pT_positive->cd(2);
  histogram_pT_negative_corresponding_positive_bin[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_negative_corresponding_positive_bin[i]->GetYaxis()->SetRangeUser(0,300);
  histogram_pT_negative_corresponding_positive_bin[i]->Draw();

  canvas_pT_positive->SaveAs(saving_path+string_pT_distribution_positive+".pdf");
  canvas_pT_positive->Close();

  TCanvas *canvas_pT_negative = new TCanvas("canvas_pT_negative", "pT distribution in each bin #mu^{-}", 900,600);
  canvas_pT_negative->Divide(2,1);
  canvas_pT_negative->cd(1);
  histogram_pT_negative_bin[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_negative_bin[i]->GetYaxis()->SetRangeUser(0,300);
  histogram_pT_negative_bin[i]->Draw();

  canvas_pT_negative->cd(2);
  histogram_pT_positive_corresponding_negative_bin[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_positive_corresponding_negative_bin[i]->GetYaxis()->SetRangeUser(0,300);
  histogram_pT_positive_corresponding_negative_bin[i]->Draw();

  string_pT_distribution_negative ="pT_negative_" + bin_number_pt[i] + "bin";
  canvas_pT_negative->SaveAs(saving_path+string_pT_distribution_negative+".pdf");
  canvas_pT_negative->Close();

  TCanvas *canvas_pT_positive_smear = new TCanvas("canvas_pT_positive_smear", "pT distribution in each bin #mu^{+}", 900,600);
  canvas_pT_positive_smear->Divide(2,1);
  canvas_pT_positive_smear->cd(1);
  histogram_pT_positive_bin_smear[i]->Draw();
  histogram_pT_positive_bin_smear[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_positive_bin_smear[i]->GetYaxis()->SetRangeUser(0,300);
  
  canvas_pT_positive_smear->cd(2);
  histogram_pT_negative_corresponding_positive_bin_smear[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_negative_corresponding_positive_bin_smear[i]->GetYaxis()->SetRangeUser(0,300);
  histogram_pT_negative_corresponding_positive_bin_smear[i]->Draw();

  string_pT_distribution_positive_smear ="pT_positive_" + bin_number_pt[i] +"bin_smear";
  canvas_pT_positive_smear->SaveAs(saving_path+string_pT_distribution_positive_smear+".pdf");
  canvas_pT_positive_smear->Close();

  TCanvas *canvas_pT_negative_smear = new TCanvas("canvas_pT_negative_smear", "pT distribution in each bin #mu^{-}", 900,600);
  canvas_pT_negative_smear->Divide(2,1);
  canvas_pT_negative_smear->cd(1);
  histogram_pT_negative_bin_smear[i]->Draw();
  histogram_pT_negative_bin_smear[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_negative_bin_smear[i]->GetYaxis()->SetRangeUser(0,300);
  
  canvas_pT_negative_smear->cd(2);
  histogram_pT_positive_corresponding_negative_bin_smear[i]->GetXaxis()->SetTitle("pT (GeV)");
  histogram_pT_positive_corresponding_negative_bin_smear[i]->GetYaxis()->SetRangeUser(0,300);
  histogram_pT_positive_corresponding_negative_bin_smear[i]->Draw();

  string_pT_distribution_negative_smear ="pT_negative_" + bin_number_pt[i] + "bin_smear";
  canvas_pT_negative_smear->SaveAs(saving_path+string_pT_distribution_negative_smear+".pdf");
  canvas_pT_negative_smear->Close();
  }
}
/*
void class_reading :: evaluating_mean_pT(){
  myfile_mean_pT<<"mean of the pT distribution before smearing : #mu^{+}"<<std::endl;
  for(int i=0; i<9; i++){
  mean_pT_positive[i] = pT_positive_bin[i]->GetMean(); 
  myfile_mean_pT<<mean_pT_positive[i]<<std::endl;
  }
  myfile_mean_pT<<"mean of the pT distribution after smearing : #mu^{+}"<<std::endl;
  for(int i=0; i<9; i++){
  mean_pT_positive_smear[i] = pT_positive_bin_smear[i]->GetMean(); 
  myfile_mean_pT<<mean_pT_positive_smear[i]<<std::endl;
  }
  myfile_mean_pT<<"mean of the pT distribution before smearing : #mu^{-}"<<std::endl;
  for(int i=0; i<9; i++){
  mean_pT_negative[i] = pT_negative_bin[i]->GetMean(); 
  myfile_mean_pT<<mean_pT_negative[i]<<std::endl;
  }
  myfile_mean_pT<<"mean of the pT distribution after smearing : #mu^{-}"<<std::endl;
  for(int i=0; i<9; i++){
  mean_pT_negative_smear[i] = pT_negative_bin_smear[i]->GetMean(); 
  myfile_mean_pT<<mean_pT_negative_smear[i]<<std::endl;
  }
}


void class_reading :: total_entries(){
  total_entries_positive = 0; 
  total_entries_negative = 0; 
  total_entries_positive_smear = 0; 
  total_entries_negative_smear = 0; 
  for(int i=0; i<9;i++){
  total_entries_positive = n_positive_bin[i] + total_entries_positive;
  total_entries_negative = n_negative_bin[i] + total_entries_negative;
  total_entries_positive_smear = n_positive_bin_smear[i] + total_entries_positive_smear;
  total_entries_negative_smear = n_negative_bin_smear[i] + total_entries_negative_smear;
  }
   if(debug_entries == true){
    std::cout<<"total number of entries #mu^{+} "<<total_entries_positive<<std::endl;
    std::cout<<"total number of entries #mu^{-} "<<total_entries_negative<<std::endl;
    std::cout<<"total number of entries #mu^{+} after smearing "<<total_entries_positive_smear<<std::endl;
    std::cout<<"total number of entries #mu^{-} after smearing "<<total_entries_negative_smear<<std::endl;
   }
    myfile_n_entries.open("total_entries.txt",std::ios_base::app);
    myfile_n_entries<<"total number of entries #mu^{+} "<<total_entries_positive<<std::endl;
    myfile_n_entries<<"total number of entries #mu^{-} "<<total_entries_negative<<std::endl;
    myfile_n_entries<<"total number of entries #mu^{+} after smearing "<<total_entries_positive_smear<<std::endl;
    myfile_n_entries<<"total number of entries #mu^{-} after smearing "<<total_entries_negative_smear<<std::endl;

    myfile_n_entries<<std::endl;

    myfile_n_entries<<"total number of entries above 200 GeV before smearing #mu^{+} "<<count_positive_mu<<std::endl;
    myfile_n_entries<<"total number of entries above 200 GeV after smearing #mu^{+} "<<count_positive_mu_smear<<std::endl;
    myfile_n_entries<<std::endl;
    myfile_n_entries<<"total number of entries #mu^{-} above 200 GeV before smearing "<<count_negative_mu<<std::endl;
    myfile_n_entries<<"total number of entries #mu^{-} above 200 GeV after smearing "<<count_negative_mu_smear<<std::endl;
    
}

void class_reading :: plot_entries(){
  TGraph *gr_positive_entries = new TGraph(9,pt_list_mean,n_positive_bin);
  TGraph *gr_negative_entries = new TGraph(9,pt_list_mean,n_negative_bin);
  TGraph *gr_positive_entries_smear = new TGraph(9,pt_list_mean,n_positive_bin_smear);
  TGraph *gr_negative_entries_smear = new TGraph(9,pt_list_mean,n_negative_bin_smear);

   gr_positive_entries->GetYaxis()->SetTitle("Entries");
   gr_positive_entries->GetXaxis()->SetTitle("pT (GeV)");
   gr_positive_entries->SetTitle("#mu^{+} before smearing");
 
   gr_negative_entries->GetYaxis()->SetTitle("Entries");
   gr_negative_entries->GetXaxis()->SetTitle("pT (GeV)");
   gr_negative_entries->SetTitle("#mu^{-} before smearing");
 
   gr_positive_entries_smear->GetYaxis()->SetTitle("Entries");
   gr_positive_entries_smear->GetXaxis()->SetTitle("pT (GeV)");
   gr_positive_entries_smear->SetTitle("#mu^{+} after smearing");
 
   gr_negative_entries_smear->GetYaxis()->SetTitle("Entries");
   gr_negative_entries_smear->GetXaxis()->SetTitle("pT (GeV)");
   gr_negative_entries_smear->SetTitle("#mu^{-} after smearing");
 
  TCanvas *canvas_entries_positive = new TCanvas("canvas_entries_positive","No. of entries in pT bin : #mu^{+}");
  canvas_entries_positive->Divide(2,1);
  canvas_entries_positive->cd(1);
  gr_positive_entries->Draw("AP*");
  canvas_entries_positive->cd(2);
  gr_positive_entries_smear->Draw("AP*");

  canvas_entries_positive->SaveAs(saving_path+"entries_positive.pdf");

  TCanvas *canvas_entries_negative = new TCanvas("canvas_entries_negative","No. of entries in pT bin : #mu^{-}");
  canvas_entries_negative->Divide(2,1);
  canvas_entries_negative->cd(1);
  gr_negative_entries->Draw("AP*");
  canvas_entries_negative->cd(2);
  gr_negative_entries_smear->Draw("AP*");
  canvas_entries_negative->SaveAs(saving_path+"entries_negative.pdf");

class_reading :: class_reading(TString input_filename){

	TFile * f = TFile::Open(input_filename);
	gen_mass = (TH1F*) f->Get("GENZ_mass"); 
}*/


void class_reading :: initializing(){
//initialize all the 9 histograms positive and negative one before and afer smearing
   histogram_gen_Zmass = (TH1F*) file->Get("hist_gen_Zmass"); 
   histogram_gen_Zmass_reconstruct= (TH1F*) file->Get("hist_gen_Zmass_reconstruct"); 
   histogram_gen_positive_mu= (TH1F*) file->Get("hist_gen_positive_mu") ;
   histogram_id_positive_mu= (TH1F*) file->Get("hist_id_positive_mu"); 
   histogram_gen_negative_mu= (TH1F*) file->Get("hist_gen_negative_mu"); 
   histogram_id_negative_mu= (TH1F*) file->Get("hist_id_negative_mu"); 
	 histogram_eta_gen_positive_mu= (TH1F*) file->Get("hist_eta_gen_positive_mu"); 
	 histogram_eta_gen_negative_mu= (TH1F*) file->Get("hist_eta_gen_negative_mu"); 
	 histogram_phi_gen_positive_mu= (TH1F*) file->Get("hist_phi_gen_positive_mu"); 
	 histogram_phi_gen_negative_mu= (TH1F*) file->Get("hist_phi_gen_negative_mu"); 
   histogram_gen_positive_mu_smear= (TH1F*) file->Get("hist_gen_positive_mu_smear"); 
   histogram_id_positive_mu_smear= (TH1F*) file->Get("hist_id_positive_mu_smear"); 
   histogram_gen_negative_mu_smear= (TH1F*) file->Get("hist_gen_negative_mu_smear"); 
   histogram_id_negative_mu_smear= (TH1F*) file->Get("hist_id_negative_mu_smear"); 
   histogram_reco_Zmass= (TH1F*) file->Get("hist_reco_Zmass"); 
   histogram_pt_gen_Z= (TH1F*) file->Get("hist_pt_gen_Z"); 
   histogram_eta_gen_Z= (TH1F*) file->Get("hist_eta_gen_Z"); 
   histogram_phi_gen_Z= (TH1F*) file->Get("hist_phi_gen_Z"); 

   for(int i=0; i<9; i++){
		TString title_Zmass_positive_pt = TString::Format("histogram_Zmass_positive_pt[%d]",i);
		TString title_Zmass_negative_pt = TString::Format("histogram_Zmass_negative_pt[%d]",i);
		TString title_Zmass_positive_pt_smear = TString::Format("histogram_Zmass_positive_pt_smear[%d]",i);
		TString title_Zmass_negative_pt_smear = TString::Format("histogram_Zmass_negative_pt_smear[%d]",i);
  
		histogram_Zmass_positive_pt[i] = (TH1F*) file->Get(title_Zmass_positive_pt); 
    histogram_Zmass_negative_pt[i] = (TH1F*) file->Get(title_Zmass_negative_pt);
    histogram_Zmass_positive_pt_smear[i] = (TH1F*) file->Get(title_Zmass_positive_pt_smear);
    histogram_Zmass_negative_pt_smear[i] = (TH1F*) file->Get(title_Zmass_negative_pt_smear);
   
    // plot pT distribuiton for each bin before and after smearing and see the average mean of the distribution
   
		TString title_pT_positive_bin = TString::Format("pT_positive_bin[%d]",i);
		TString title_pT_negative_bin = TString::Format("pT_negative_bin[%d]",i);
		TString title_pT_positive_bin_smear = TString::Format("pT_positive_bin_smear[%d]",i);
		TString title_pT_negative_bin_smear = TString::Format("pT_negative_bin_smear[%d]",i);
		TString title_pT_negative_corresponding_positive_bin = TString::Format("pT_negative_corresponding_positive_bin[%d]",i);
		TString title_pT_positive_corresponding_negative_bin = TString::Format("pT_positive_corresponding_negative_bin[%d]",i);
		TString title_pT_negative_corresponding_positive_bin_smear = TString::Format("pT_negative_corresponding_positive_bin_smear[%d]",i);
		TString title_pT_positive_corresponding_negative_bin_smear = TString::Format("pT_positive_corresponding_negative_bin_smear[%d]",i);
	

    histogram_pT_positive_bin[i] = (TH1F*) file->Get(title_pT_positive_bin) ; 
    histogram_pT_negative_bin[i] = (TH1F*) file->Get(title_pT_negative_bin) ; 
    histogram_pT_positive_bin_smear[i] = (TH1F*) file->Get(title_pT_positive_bin_smear) ; 
    histogram_pT_negative_bin_smear[i] = (TH1F*) file->Get(title_pT_negative_bin_smear) ; 
    histogram_pT_negative_corresponding_positive_bin[i] = (TH1F*) file->Get(title_pT_negative_corresponding_positive_bin) ; 
    histogram_pT_positive_corresponding_negative_bin[i] = (TH1F*) file->Get(title_pT_positive_corresponding_negative_bin) ; 
    histogram_pT_negative_corresponding_positive_bin_smear[i] = (TH1F*) file->Get(title_pT_negative_corresponding_positive_bin_smear) ; 
    histogram_pT_positive_corresponding_negative_bin_smear[i] = (TH1F*) file->Get(title_pT_positive_corresponding_negative_bin_smear) ; 

	}
}



  void derived_class_reading :: plot_histograms(){
		class_reading::plotting_hist(histogram_gen_positive_mu," pT : #mu^{+} ","gen_pT_positive_mu","pT (GeV)");
		class_reading::plotting_hist(histogram_eta_gen_positive_mu," #eta : #mu^{+} ","gen_eta_positive_mu","#eta");
		class_reading::plotting_hist(histogram_phi_gen_positive_mu," #phi : #mu^{+} ","gen_phi_positive_mu","#phi");
		class_reading::plotting_hist(histogram_gen_positive_mu," pT : #mu^{+} ","gen_pT_positive_mu","pT (GeV)");
		class_reading::plotting_hist(histogram_id_positive_mu,"id : #mu^{+} ","id_positive_mu", " ");
 
		class_reading::plotting_hist(histogram_gen_negative_mu," pT : #mu^{-} ","gen_pT_negative_mu","pT (GeV)");
		class_reading::plotting_hist(histogram_eta_gen_negative_mu," #eta : #mu^{-} ","gen_eta_negative_mu","#eta");
		class_reading::plotting_hist(histogram_phi_gen_negative_mu," #phi : #mu^{-} ","gen_phi_negative_mu","#phi");

		class_reading::plotting_hist(histogram_id_negative_mu,"id : #mu^{-} ","id_negative_mu", " ");
		class_reading::plotting_hist(histogram_gen_Zmass,"GEN J/#Psi mass ","gen_Z_distribution", " ");
		class_reading::plotting_hist(histogram_gen_Zmass_reconstruct," J/#Psi mass (GEN) reconstructed","gen_Z_reconstructed_distribution", " ");
    class_reading::plotting_hist(histogram_reco_Zmass,"reco J/#Psi mass ","reco_Z_distribution", " ");
 
	 // plotting Jpsi ptogram, eta and phi values
    class_reading::plotting_hist(histogram_pt_gen_Z,"pT GEN J/#Psi  ","gen_Z_pt_distribution", " ");
    class_reading::plotting_hist(histogram_eta_gen_Z,"eta GEN J/#Psi  ","gen_Z_eta_distribution", " ");
    class_reading::plotting_hist(histogram_phi_gen_Z,"phi GEN J/#Psi  ","gen_Z_phi_distribution", " ");
//  class_reading:: if(plotting_on == true

    class_reading::plotting_hist(histogram_gen_positive_mu_smear," pT of #mu^{+} after smearing ","gen_pT_positive_mu_smear","pT (GeV)");
    class_reading::plotting_hist(histogram_id_positive_mu_smear,"id of #mu^{+} after smearing ","id_positive_mu_smear", " ");
 
    class_reading::plotting_hist(histogram_gen_negative_mu_smear," pT of #mu^{-} after smearing ","gen_pT_negative_mu_smear","pT (GeV)");
    class_reading::plotting_hist(histogram_id_negative_mu_smear,"id of #mu^{-} after smearing ","id_negative_mu_smear", " ");


	   // plotting Jpsi mass for positive bins 	
	   for (int i=0; i<9;i++){
     TString title = "#mu^{+} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] +" (GeV)" + " J/#Psi mass before smearing ";
     TString saving_name_hist = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_positive";
     cout<<"saving name of histogram "<<saving_name_hist<<endl;
     class_reading::plotting_hist(histogram_Zmass_positive_pt[i],title,saving_name_hist, "m(#mu^{+}#mu^{-})");
     } 
   
	  // plotting Jpsi mass for negative bins 	
     for (int i=0; i<9;i++){
     TString title = "#mu^{-} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)" +" J/#Psi mass before smearing ";
     TString saving_name_hist = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_negative";
     cout<<"saving name of histogram "<<saving_name_hist<<endl;
     class_reading::plotting_hist(histogram_Zmass_negative_pt[i],title,saving_name_hist, "m(#mu^{+}#mu^{-})");
     }

  }

  // this to do the caluclation for mean and mean error, it will plot all the graphs for mean Jpsi mass before smearing and mean Jpsi mass after smearing 
	// it further will plot the difference in the mass before and aftere smearing
  void derived_class_reading :: mean_sigma_calculation(){	

		class_reading::graph_mean_pT_combine();
		class_reading::graph_mean_pT_combine_smear();
		class_reading::graph_diff_Zmass_gen_reco();
	}	
   void derived_class_reading :: fitting_histograms(){

		 // fitting inclusive distribution
    mean_sigma_Zmass =  class_reading::plotting_fitting_mass_BW(histogram_gen_Zmass,"gen_Zmass_fit", "Gen J/#Psi before smearing");
    mean_sigma_Zmass_reconstruct = class_reading::plotting_fitting_mass_BW(histogram_gen_Zmass_reconstruct,"gen_Zmass_reconstruct_fit", "Gen J/#Psi before smearing (reconstructed)");
    mean_sigma_Zmass_smear =  class_reading::plotting_fitting_mass_BW_smearing(histogram_reco_Zmass,"gen_Zmass_fit_smear", "Gen J/#Psi after smearing");
  
	 // to perform a gaussian fit to distribution	
//	  mean_sigma_Zmass_gauss =  class_reading::plotting_fitting_mass_gauss(histogram_gen_Zmass,"gen_Zmass_fit","Gen J/#Psi before smearing");
    mean_sigma_Zmass_gauss_smear =  class_reading::plotting_fitting_mass_gauss_smearing(histogram_reco_Zmass,"gen_Zmass_fit_smear", "Gen J/#Psi after smearing");

		// To perform a DSCB shift
		mean_sigma_Zmass_DSCB =  class_reading::plotting_fitting_mass_DSCB(histogram_gen_Zmass,"gen_Zmass_fit","Gen J/#Psi before smearing");
//		mean_sigma_Zmass_DSCB_smear =  class_reading::plotting_fitting_mass_DSCB_smearing(histogram_reco_Zmass,"gen_Zmass_fit_smear", "Gen J/#Psi after smearing");



	  // fitting J/psi mass with mu^{+} selection	
     for (int i=0; i<9;i++){
     TString saving_name_hist_fit = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_fit_positive";
     TString title = "#mu^{+} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)"+ " before smearing" ;

		 mean_sigma_positive_pt[i] = class_reading::plotting_fitting_mass_BW(histogram_Zmass_positive_pt[i],saving_name_hist_fit,title);
		 //mean_sigma_positive_pt[i] = class_reading::plotting_fitting_mass_DSCB(histogram_Zmass_positive_pt[i],saving_name_hist_fit,title);

//		 mean_sigma_positive_pt[i] = class_reading::plotting_fitting_mass_DSCB_bin(histogram_Zmass_positive_pt[i],saving_name_hist_fit,title);

     TString saving_name_hist_fit_smear = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_smear_fit_positive";
     TString title_smear = "#mu^{+} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)"+ " after smearing" ;

		 mean_sigma_positive_reco_pt[i] = class_reading::plotting_fitting_mass_gauss_smearing(histogram_Zmass_positive_pt_smear[i],saving_name_hist_fit_smear,title_smear);

//		 mean_sigma_positive_reco_pt[i] = class_reading::plotting_fitting_mass_DSCB_bin_smearing(histogram_Zmass_positive_pt_smear[i],saving_name_hist_fit_smear,title_smear);

    }
 
	  // fitting J/psi mass with mu^{-} selection	
     for (int i=0; i<9;i++){
     TString saving_name_hist_fit = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_fit_negative";
     TString title = "#mu^{-} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)"+ " before smearing" ;

		 mean_sigma_negative_pt[i] = class_reading::plotting_fitting_mass_BW(histogram_Zmass_negative_pt[i],saving_name_hist_fit,title);
//		 mean_sigma_negative_pt[i] = class_reading::plotting_fitting_mass_DSCB(histogram_Zmass_negative_pt[i],saving_name_hist_fit,title);

//		 mean_sigma_negative_pt[i] = class_reading::plotting_fitting_mass_DSCB_bin(histogram_Zmass_negative_pt[i],saving_name_hist_fit,title);

     TString saving_name_hist_fit_smear = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_smear_fit_negative";
     TString title_smear = "#mu^{-} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)"+ " after smearing" ;
     mean_sigma_negative_reco_pt[i] = class_reading::plotting_fitting_mass_gauss_smearing(histogram_Zmass_negative_pt_smear[i],saving_name_hist_fit_smear,title_smear);

//		 mean_sigma_negative_reco_pt[i] = class_reading::plotting_fitting_mass_DSCB_bin_smearing(histogram_Zmass_negative_pt_smear[i],saving_name_hist_fit_smear,title_smear);
  }

}

// this function is specifically to plot the inclusive and binned distribution on the same canvas and seee how well they agree 
void derived_class_reading :: plotting_inclusive_bin(){

	// positive bin
    for (int i=0; i<9;i++){
  	 TString title_plot = "#mu^{+} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)" +" before smearing ";
     TString saving_name_plot_fullrange = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_positive_distribution_fullrange";
		 class_reading::plotting_mass_distribution(histogram_gen_Zmass, histogram_Zmass_positive_pt[i],3.0968, 3.09704 , title_plot, saving_name_plot_fullrange, "m(#mu^{+}#mu^{-})"); 
	 } 

		// negative bin 
    for (int i=0; i<9;i++){
     TString saving_name_plot_fullrange = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_negative_distribution_fullrange";
     TString title_plot = "#mu^{-} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)" +" before smearing ";
		 class_reading::plotting_mass_distribution(histogram_gen_Zmass, histogram_Zmass_negative_pt[i], 3.0968,3.09704, title_plot, saving_name_plot_fullrange, "m(#mu^{+}#mu^{-})"); 
    }

    for (int i=0; i<9;i++){
  	 TString title_plot = "#mu^{+} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)" +" after  smearing ";
     TString saving_name_plot_fullrange = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_positive_distribution_smeared_fullrange";
		 class_reading::plotting_mass_distribution(histogram_reco_Zmass, histogram_Zmass_positive_pt_smear[i],3.02, 3.18 , title_plot, saving_name_plot_fullrange, "m(#mu^{+}#mu^{-})"); 
    }
 
    for (int i=0; i<9;i++){
		 TString title_plot = "#mu^{-} : " + pt_list_symbol[i] + "  #leq pT < " + pt_list_symbol[i+1] + " (GeV)" +" after  smearing ";
     TString saving_name_plot_fullrange = "Zmass_"+class_reading::bin_number_pt[i]+"_pt_negative_distribution_smeared_fullrange";
		 class_reading::plotting_mass_distribution(histogram_reco_Zmass, histogram_Zmass_negative_pt_smear[i],3.02, 3.18 , title_plot, saving_name_plot_fullrange, "m(#mu^{+}#mu^{-})"); 

		 }

}

     void derived_class_reading :: saving_text_file(){

    	ofstream myfile;
//		  ofstream	myfile_n_entries;
      myfile.open("mean_Zmass.txt",std::ios_base::app);
//      myfile_n_entries.open("entries_bin.txt",std::ios_base::app);
      myfile<<"smearing percentage : 1%"<<std::endl;  

       myfile<<"mean and sigma before smearing"<<std::endl;
       myfile<<"gen Z mass fit BW"<<std::endl; 
       myfile<<"mean Z"<<"\t\t  "<<"mean error Z"<<std::endl;
       myfile<<mean_sigma_Zmass.first<<"\t\t  "<<mean_sigma_Zmass.second<<std::endl;
 
       myfile<<"gen Z mass fit DSCB"<<std::endl; 
       myfile<<"mean Z"<<"\t\t  "<<"mean error Z"<<std::endl;
       myfile<<mean_sigma_Zmass_DSCB.first<<"\t\t  "<<mean_sigma_Zmass_DSCB.second<<std::endl;

       myfile<<"mean and sigma after smearing"<<std::endl;
       myfile<<"gen Z mass fit BW"<<std::endl; 
       myfile<<"mean Z"<<"\t\t  "<<"mean error Z"<<std::endl;
       myfile<<mean_sigma_Zmass_smear.first<<"\t\t  "<<mean_sigma_Zmass_smear.second<<std::endl;
 
       myfile<<"gen Z mass fit gauss"<<std::endl; 
       myfile<<"mean Z"<<"\t\t  "<<"mean error Z"<<std::endl;
       myfile<<mean_sigma_Zmass_gauss_smear.first<<"\t\t  "<<mean_sigma_Zmass_gauss_smear.second<<std::endl;

       myfile<<"for positive muons in pT binned"<<std::endl; 
       myfile<<"bin no."<<"\t \t"<<"mean Z (before)"<<"\t \t"<<"mean Z (after)"<<"\t\t  "<<std::endl;

 /*      myfile_n_entries<<"entries in all positive 9 pT bins before smearing"<<endl;
       for(int i = 0 ; i<9; i++){
       myfile_n_entries<<n_positive_bin[i]<<endl;
       }
        myfile_n_entries<<"entries in all negative 9 pT bins before smearing"<<endl;
       for(int i = 0 ; i<9; i++){
       myfile_n_entries<<n_negative_bin[i]<<endl;
       }
  */    
       for(int i= 0 ; i<9;i++){
				 class_reading::mean_positive_pt[i] = mean_sigma_positive_pt[i].first;
				 class_reading::mean_positive_pt_error[i] = mean_sigma_positive_pt[i].second;
				 class_reading::mean_positive_reco_pt[i] = mean_sigma_positive_reco_pt[i].first;
				 class_reading::mean_positive_reco_pt_error[i] = mean_sigma_positive_reco_pt[i].second;
				 class_reading::diff_Zmass_positive_gen_reco[i] = (class_reading::mean_positive_reco_pt[i] - class_reading::mean_positive_pt[i]) / class_reading::mean_positive_pt[i] ; 

        // if y = ab , #sigma(y) = y *( (sigma(a)/a) + (sigma(b)/b) )
				 class_reading::diff_Zmass_positive_gen_reco_error[i] = class_reading::diff_Zmass_positive_gen_reco[i] *sqrt( pow((class_reading::mean_positive_reco_pt_error[i] / class_reading::mean_positive_reco_pt[i]),2) + pow((class_reading::mean_positive_pt_error[i] / class_reading::mean_positive_pt[i]),2 )) ;  	
	
       myfile<<i<<"\t \t "<<class_reading::mean_positive_pt[i]<<"\t\t  "<<class_reading::mean_positive_reco_pt[i]<<"\t \t"<<class_reading::mean_positive_pt_error[i]<<"\t \t"<<class_reading::mean_positive_reco_pt_error[i]<<std::endl; 
       }
     
   /*    myfile_n_entries<<"entries in all positive 9 pt bins after smearing"<<endl;
       for(int i = 0 ; i<9; i++){
       myfile_n_entries<<n_positive_bin_smear[i]<<endl;
       }
        myfile_n_entries<<"entries in all negative 9 pt bins after smearing"<<endl;
       for(int i = 0 ; i<9; i++){
       myfile_n_entries<<n_negative_bin_smear[i]<<endl;
       }
    */   

       myfile<<"for negative muons in pt binned"<<std::endl; 

       myfile<<"bin no."<<"\t \t"<<"mean Z (before)"<<"\t \t"<<"mean Z (after)"<<std::endl;
       for(int i= 0 ; i<9;i++){
				 class_reading::mean_negative_pt[i] = mean_sigma_negative_pt[i].first;
				 class_reading::mean_negative_pt_error[i] = mean_sigma_negative_pt[i].second;
				 class_reading::mean_negative_reco_pt[i] = mean_sigma_negative_reco_pt[i].first;
				 class_reading::mean_negative_reco_pt_error[i] = mean_sigma_negative_reco_pt[i].second;
				 class_reading::diff_Zmass_negative_gen_reco[i] = (class_reading::mean_negative_reco_pt[i] - class_reading::mean_negative_pt[i]) / class_reading::mean_negative_pt[i] ; 
 
        // if y = ab , #sigma(y) = y *( (sigma(a)/a) + (sigma(b)/b) )
				 class_reading::diff_Zmass_negative_gen_reco_error[i] = class_reading::diff_Zmass_negative_gen_reco[i] *sqrt( pow((class_reading::mean_negative_reco_pt_error[i] / class_reading::mean_negative_reco_pt[i]),2) + pow((class_reading::mean_negative_pt_error[i] / class_reading::mean_negative_pt[i]),2 )) ;  
	 
				 myfile<<i<<"\t \t "<<class_reading::mean_negative_pt[i]<<"\t\t  "<<class_reading::mean_negative_reco_pt[i]<<"\t \t"<<class_reading::mean_negative_pt_error[i]<<"\t \t"<<class_reading::mean_negative_reco_pt_error[i]<<std::endl; 
      } 
     }

  void histogram_analysis(){


	 // this code is to analyze all the histograms produced for the Jpsi sample : histograms are Z mass histogram in all the pT categories, for both mu^+ and mu^- 	
   // derived class reading is the class where I am performing operation to fit the histograms , and dervied class reading is using functions from base class reading to perform this operations	
   gROOT->SetBatch(kTRUE); 	 
	
	 derived_class_reading obj;

	 // initializing is to open the root file and to read all the histograms
	 obj.initializing(); 

	 // plotting every histogram on canvas with the plot_histograms function	
   obj.plot_histograms(); 
 
 	 // plotting inclusive and binned distribution on same canvas 
	 obj.plotting_inclusive_bin();

	 // to do BW and DSCB fits and find mean and sigma values  and this mean and sigma values are further stored and used to make final plots in the functio nmean_sigma_calculation
	 obj.fitting_histograms();

	 // saving text file also has the calculation for mean and sigma which will be used furthere
	 obj.saving_text_file();
  
	 // graph to plot mean z mass with pT
   obj.mean_sigma_calculation();
   
	 // to save the pT distributions for individual histograms
	 obj.saving_histogram_pT();

	}
