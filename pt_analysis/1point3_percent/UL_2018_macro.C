#define UL_2018_macro_cxx
#include "UL_2018_macro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void UL_2018_macro::Loop()
{
//	gSystem->Load("RooMyPDF_DSCB.h");
 //gSystem->Load("RooMyPDF_DSCB.cxx");


	 gROOT->SetBatch(kTRUE); 
	 
	 //   In a ROOT session, you can do:
//      root> .L UL_2018_macro.C
//      root> UL_2018_macro t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast(); 
   Long64_t nbytes = 0, nb = 0;
 
//	   for (Long64_t jentry=0; jentry<nentries;jentry++) {
for (Long64_t jentry=0; jentry<50000;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) { 
        std::cout<<"entry is breaking"<<std::endl;
        std::cout<<"entry number "<<jentry<<std::endl; 
        break;
      }
      nb = fChain->GetEntry(jentry);   
      nbytes += nb;
    
		  // information GEN Z 	
      //std::cout<<"Gen Z mass size"<<*GENZ_mass[0]<<std::endl;	
			if(GENZ_mass->size()==0 || GENZ_mass->size()==2) continue;
			count_mass++;
			count_total_events++;
			if(GENZ_DaughtersId->at(0) == 443) count_size_Z++;
/*			if(GENZ_mass->size() ==1) {
       std::cout<<"size"<<GENZ_mass->size()<<" mass : "<<GENZ_mass->at(0)<<"  GEN Z mother and daughters ID "<<GENZ_MomId->at(0)<<"   "<<GENZ_DaughtersId->at(0)<<std::endl;
			 std::cout<<"corr lep "<<lep_pt->size()<<"GEN lep "<<GENlep_pt->size()<<std::endl;
			 std::cout<<"mother id of lep "<<GENlep_MomId->at(0)<<"  "<<GENlep_MomId->at(1)<<std::endl;
			 std::cout<<"mother mother id of lep "<<GENlep_MomMomId->at(0)<<"  "<<GENlep_MomMomId->at(1)<<std::endl;
			}
*/			
      //std::cout<<"Gen Z mass daughters ID"<<GENZ_DaughtersId->at(0)<<" mass "<<GENZ_mass->at(0)<<std::endl;	
			//std::cout<<"when size is "<<GENZ_mass->size()<<" lep id  :"<<Z_pt->size()<<" lep id mass "<<lep_pt->size()<<std::endl;
			//"   "<<lep_id->at(0)<<std::endl; 
			if(GENZ_DaughtersId->at(0) != 443) continue;
			count_daughter++;
      hist_gen_Zmass->Fill(GENZ_mass->at(0));      
 //     hist_MC_Zmass->Fill(Z_mass->at(0));      
	
			// reconstruct gen Z mass 	
		  reconstruct_genZmass();	
			// when positive muon has higher pT than negative one
      if(debug_base) std::cout<<"number of entries in the event"<<std::endl;
  	  if(debug_base) std::cout<<GENlep_id->size()<<std::endl;	
  	  n_muons = GENlep_id->size();

			if(debug_program) std::cout<<"GEN  Z mass in start "<<GENZ_mass->at(0)<<std::endl;

      hist_pt_gen_Z->Fill(GENZ_pt->at(0));      
      hist_eta_gen_Z->Fill(GENZ_eta->at(0));      
      hist_phi_gen_Z->Fill(GENZ_phi->at(0));      
	
//			hist_pt_MC_Z->Fill(Z_pt->at(0));      
//      hist_eta_MC_Z->Fill(Z_eta->at(0));      
//      hist_phi_MC_Z->Fill(Z_phi->at(0));      
	 // selecting positive muons		 
      if(GENlep_id->at(0)==13)    
      {
       if(debug_base==true) { std::cout<<" positive muon lep_id "<<GENlep_id->at(0)<<std::endl;
                            std::cout<<"pt positive muon pt :  "<<GENlep_pt->at(0)<<std::endl; } 

        // dividing in pT bins
       

       // if you just want to analyze pT 
       //
         gen_pT_positive_value = GENlep_pt->at(0);
         gen_pT_negative_value = GENlep_pt->at(1);
         hist_gen_positive_mu->Fill(GENlep_pt->at(0));
         hist_id_positive_mu->Fill(GENlep_id->at(0));        
         hist_gen_negative_mu->Fill(GENlep_pt->at(1));
         hist_id_negative_mu->Fill(GENlep_id->at(1)); 
         
         hist_eta_gen_positive_mu->Fill(GENlep_eta->at(0));
         hist_eta_gen_negative_mu->Fill(GENlep_eta->at(1));
         hist_phi_gen_positive_mu->Fill(GENlep_phi->at(0));
         hist_phi_gen_negative_mu->Fill(GENlep_phi->at(1));


       //gen pT of individual lepton in individual bin

       if(analyze_pT == true){
         if(debug_program) std::cout<<" first muon is positive"<<std::endl; 
			   if(debug_program) std::cout<<"GEN  pT for positive muon  "<<GENlep_pt->at(0)<<std::endl;
         for(int i=0; i<9; i++){
 //      std::cout<<"filling indiviidual bins"<<std::endl;
          if(pt_list[i] <= GENlep_pt->at(0) && GENlep_pt->at(0) < pt_list[i+1]){
			    if(debug_program) std::cout<<"GEN  pT for positive muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" with gen PT "<<GENlep_pt->at(0)<<std::endl;
			    if(debug_program) std::cout<<"GEN  Z mass for positive muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
            histogram_filling_positive_pt(i);
	//	    std::cout<<"gen Z mass   next line    "<<GENZ_mass->at(0)<<std::endl;
           n_positive_bin[i]++;
           pT_positive_bin[i]->Fill(GENlep_pt->at(0));
           pT_negative_corresponding_positive_bin[i] ->Fill(GENlep_pt->at(1));  
          }
         }
    
     
        if(debug_count == true){ 
         if(GENlep_pt->at(0) >20){
          count_positive_mu++;
         }
        }
			   if(debug_program) std::cout<<"GEN  pT for negative muon  "<<GENlep_pt->at(1)<<std::endl;
         // filling Z mass for negative pt bins 
         for(int i=0; i<9; i++){
         if(pt_list[i] <= GENlep_pt->at(1) && GENlep_pt->at(1) < pt_list[i+1]){
 			   if(debug_program) std::cout<<"GEN  pT for negative muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" with gen PT "<<GENlep_pt->at(1)<<std::endl;
			   if(debug_program) std::cout<<"GEN  Z mass for negative muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           histogram_filling_negative_pt(i);
           n_negative_bin[i]++;
           pT_negative_bin[i] ->Fill(GENlep_pt->at(1));
           pT_positive_corresponding_negative_bin[i] ->Fill(GENlep_pt->at(0));  
          }
         }

        if(debug_count == true){ 
          if(GENlep_pt->at(1) >20){
          count_negative_mu++;
          }
        } 
       }
 
       if(smearing == true){
         // smear individual pT
         smearing_pT();

         // reconstruct mass of Z from Lorentzvector for individual lepton 
         reconstruct_Zmass();
         
        smeared_pT_positive_value = GENlep_pt->at(0); 
        smeared_pT_negative_value = GENlep_pt->at(1); 
        // fill the reconstructed Z mass into a histogram
        hist_reco_Zmass ->Fill(massZ_reconstruct);
        hist_gen_positive_mu_smear->Fill(GENlep_pt->at(0));
        hist_id_positive_mu_smear->Fill(GENlep_id->at(0));        
        hist_gen_negative_mu_smear->Fill(GENlep_pt->at(1));
        hist_id_negative_mu_smear->Fill(GENlep_id->at(1)); 

       }
     
      // have to be run after the pT is smeared  
       if(analyze_pT == true && smearing == true){

         for(int i=0; i<9; i++){
          if(pt_list[i] <= GENlep_pt->at(0) && GENlep_pt->at(0) < pt_list[i+1]){
           histogram_filling_positive_pt_smear(i);

					 if(debug_program) std::cout<<std::endl;
					 if(debug_program) std::cout<<" after smearing  Z mass for positive muon is falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           n_positive_bin_smear[i]++;
           pT_positive_bin_smear[i] ->Fill(GENlep_pt->at(0));
           pT_negative_corresponding_positive_bin_smear[i] ->Fill(GENlep_pt->at(1));  
          }
         }
       if(debug_count == true){ 
         if(GENlep_pt->at(0) >20){
          count_positive_mu_smear++;
         }
       }
  
         // filling Z mass for negative pt bins 
         for(int i=0; i<9; i++){
          if(pt_list[i] <= GENlep_pt->at(1) && GENlep_pt->at(1) < pt_list[i+1]){
			   if(debug_program) std::cout<<" after smearing  Z mass for negative muon is falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           histogram_filling_negative_pt_smear(i);
           n_negative_bin_smear[i]++;
           pT_negative_bin_smear[i] ->Fill(GENlep_pt->at(1));
           pT_positive_corresponding_negative_bin_smear[i]->Fill(GENlep_pt->at(0));  
          }
         }
         if(debug_count == true){ 
         if(GENlep_pt->at(1) >20){
          count_negative_mu_smear++;
         }
       }

       }
  }
  // Storing positive muons with pT less than negative muons
       if(GENlep_id->at(1)==13)    
      {
         if(debug_base==true) { std::cout<<" positive muon id "<<GENlep_id->at(1)<<std::endl;
                                std::cout<<"pT positive and pT negative muon :  "<<GENlep_pt->at(1)<<"  "<<GENlep_pt->at(0)<<std::endl; } 
         if(debug_program) std::cout<<" first muon is negative"<<std::endl; 
			   if(debug_program) std::cout<<"GEN  pT for positive muon  "<<GENlep_pt->at(1)<<std::endl;
      // dividing in pT bins
         gen_pT_positive_value = GENlep_pt->at(1); 
         gen_pT_negative_value = GENlep_pt->at(0); 

         hist_gen_positive_mu->Fill(GENlep_pt->at(1));
         hist_id_positive_mu->Fill(GENlep_id->at(1));        
 
         hist_gen_negative_mu->Fill(GENlep_pt->at(0));
         hist_id_negative_mu->Fill(GENlep_id->at(0)); 
 
				 hist_eta_gen_positive_mu->Fill(GENlep_eta->at(1));
         hist_eta_gen_negative_mu->Fill(GENlep_eta->at(0));
         hist_phi_gen_positive_mu->Fill(GENlep_phi->at(1));
         hist_phi_gen_negative_mu->Fill(GENlep_phi->at(0));


 

         if(analyze_pT == true){

         for(int i=0; i<9; i++){
          if(pt_list[i] <= GENlep_pt->at(1) && GENlep_pt->at(1) < pt_list[i+1]){
  		   if(debug_program) std::cout<<"GEN  pT for positive muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" with gen PT "<<GENlep_pt->at(1)<<std::endl;
			   if(debug_program) std::cout<<"GEN  Z mass for positive muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           histogram_filling_positive_pt(i);
           n_positive_bin[i]++;
           pT_positive_bin[i] ->Fill(GENlep_pt->at(1));
           pT_negative_corresponding_positive_bin[i]->Fill(GENlep_pt->at(0));  
          }
         }


        if(debug_count == true){ 
           if(GENlep_pt->at(1) >20){
           count_positive_mu++;
          }
        }
 
			   if(debug_program) std::cout<<"GEN  pT for negative muon  "<<GENlep_pt->at(0)<<std::endl;

         // filling Z mass for negative pt bins 
         for(int i=0; i<9; i++){
          if(pt_list[i] <= GENlep_pt->at(0) && GENlep_pt->at(0) < pt_list[i+1]){
           histogram_filling_negative_pt(i);	
   		   if(debug_program) std::cout<<"GEN  pT for negative muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" with gen PT "<<GENlep_pt->at(0)<<std::endl;
			   if(debug_program) std::cout<<"GEN  Z mass for negative muon falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           n_negative_bin[i]++;
           pT_negative_bin[i] ->Fill(GENlep_pt->at(0));
           pT_positive_corresponding_negative_bin[i]->Fill(GENlep_pt->at(1));  
          }
         }
      
        
        if(debug_count == true){ 
         if(GENlep_pt->at(0) >20){
          count_negative_mu++;
          }
        } 


       }
 
         if(smearing == true){
         // smear individual pT
         smearing_pT();

         // reconstruct mass of Z from Lorentzvector for individual lepton 
         reconstruct_Zmass();

         smeared_pT_positive_value = GENlep_pt->at(1); 
         smeared_pT_negative_value = GENlep_pt->at(0); 


        // fill the reconstructed Z mass into a histogram
        hist_reco_Zmass->Fill(massZ_reconstruct);
//        hist_gen_Zmass_smear->Fill(GENZ_mass->at(0));      
 
        hist_gen_positive_mu_smear->Fill(GENlep_pt->at(1));
        hist_id_positive_mu_smear->Fill(GENlep_id->at(1));        
        hist_gen_negative_mu_smear->Fill(GENlep_pt->at(0));
        hist_id_negative_mu_smear->Fill(GENlep_id->at(0)); 

       }

     
      // have to be run after the pT is smeared  
       if(analyze_pT == true && smearing == true){

         for(int i=0; i<9; i++){
          if(pt_list[i] <= GENlep_pt->at(1) && GENlep_pt->at(1) < pt_list[i+1]){
			   if(debug_program) std::cout<<" after smearing  Z mass for positive muon is falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           histogram_filling_positive_pt_smear(i);
           n_positive_bin_smear[i]++;
           pT_positive_bin_smear[i] ->Fill(GENlep_pt->at(1));
           pT_negative_corresponding_positive_bin_smear[i]->Fill(GENlep_pt->at(0));  
          }
         }

         // filling Z mass for negative pt bins 
         for(int i=0; i<9; i++){
          if(pt_list[i] <= GENlep_pt->at(0) && GENlep_pt->at(0) < pt_list[i+1]){
			   if(debug_program) std::cout<<" after smearing  Z mass for negative muon is falling in  range "<<pt_list[i]<<" , "<<pt_list[i+1]<<" GEN mass "<<GENZ_mass->at(0)<<std::endl;
           histogram_filling_negative_pt_smear(i);
           n_negative_bin_smear[i]++;
           pT_negative_bin_smear[i] ->Fill(GENlep_pt->at(0));
           pT_positive_corresponding_negative_bin_smear[i]->Fill(GENlep_pt->at(1));  
          }
         }
       }

      }
    // end of loop of cmparing gen Id 
    //
    // fill the event pt in a 2D histogram for migration matrix
      hist_events_pTbins_positive->Fill(gen_pT_positive_value, smeared_pT_positive_value);
      hist_events_pTbins_negative->Fill(gen_pT_negative_value, smeared_pT_negative_value);

   }
// writing histogram to the output root file
  percentage_daughter = count_size_Z/count_total_events; 
  std::cout<<" percentage times Z lep has daughter Jpsi "<<percentage_daughter<<" and number of events "<<count_size_Z<<" total events "<<count_total_events<<std::endl;

   // to find the migration bins  
   events_pTbins();
	 write_histograms();

	 
 }
// end of program
