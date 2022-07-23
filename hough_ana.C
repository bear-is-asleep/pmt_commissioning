#include <TStyle.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string> 

#define PI 3.14159265
using std::vector;
//using namespace std;

float distance(int x1, int y1, int x2, int y2){
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) * 1.0);
}

// usual parameter values: 
// threshold = 10, if not interested in anode-cathode (AC) crossing muons, recommended threshold = 20 
// max_gap = 30 
// range = 150 for AC muons, range = 100 otherwise 
// min_length = 500 
// muon_length = 0 if do not want to perform AC cut, otherwise muon_length varies from 2000 to 2500 

void hough(vector<vector<int>> coords, int threshold, int max_gap, int range, float min_length, int muon_length, int nentry, vector<vector<int>>& lines){
   //set global variables 
   TRandom3 rndgen;
   const int h = 3500; const int w = 2000; //range of hit_wire
   constexpr int accu_h = h + w + 1 ; const int accu_w = 180; 
   const int x_c = (w/2); const int y_c = (h/2); 

   //create accumulator/pointer to accumulator 
   int accu[accu_h][accu_w] = {{0}};
   int (*adata)[accu_w];
   adata = &accu[(accu_h-1)/2]; // have pointer point to the middle of the accumulator array (then we can use negative indices)

   //declare necessary vectors 
   vector<vector<int>> data = coords; // points will not be removed 
   vector<vector<int>> deaccu; //deaccumulator
   vector<vector<int>> outlines; 

   // loop over points and perform transform 
   int count = coords.size(); 
   for ( ; count>0; count--){ 
      int idx = rndgen.Uniform(count);
      int max_val = threshold-1;
      if ((coords.at(idx)).empty())
         continue; 
      int x = coords[idx][0], y = coords[idx][1], rho = 0, theta = 0;
      vector<int> v{x,y}; 
      deaccu.push_back(v);
      //loop over all angles and fill the accumulator 
      for (int j=0; j<accu_w; j++){ 
         int r = int(round((x-x_c)*cos(j*PI/accu_w) + (y-y_c)*sin(j*PI/accu_w)));
         int val = ++(adata[r][j]);       
         if (max_val < val){
            max_val = val;
            rho = r;
            theta = j*180/accu_w;
         }
      }
      if (max_val < threshold){
         (coords.at(idx)).clear(); 
         continue;
      }
      //start at point and walk the corridor on both sides 
      vector<vector<int>> endpoint(2, vector<int>(2));
      for (int k=0; k<2;k++){ 
         int i=0, gap=0;
         while (gap < max_gap){ 
            (k==0)? i++ : i--; 
            if ( (idx+i) == data.size() || (idx+i) <0) // if we reach the edges of the data set 
               break;
            if ((data.at(idx+i)).empty()) // if the point has already been removed 
               continue;
            int x1 = data[idx+i][0], y1 = int(data[idx+i][1]); 
            int last_x, diffx, last_y, diffy;
            if (endpoint[k][0]!= 0){ // ensure we don't jump large x-values 
               last_x = endpoint[k][0];
               diffx = abs(last_x - x1);
               if (diffx > 30){
                  break;
               }
            }
            int y_val = int(round((rho - (x1 - x_c)*cos(theta*PI/180.0))/sin(theta*PI/180.0) + y_c));
            if (y1 >= (y_val-range) && y1 <= (y_val+range)){
               gap = 0;
               endpoint[k] = {x1, y1};
               (coords.at(idx+i)).clear();
               (data.at(idx+i)).clear();
            }
            else
               gap++;
         } // end of while loop 
      } // end of k loop 

      // unvote from the accumulator 
      for (int n = (deaccu.size()-1); n>=0; n--){ 
         int x1 = deaccu[n][0], y1 = int(deaccu[n][1]);
         int y_val = int(round((rho - (x1 - x_c)*cos(theta*PI/180.0))/sin(theta*PI/180.0) + y_c));
         if (y1 >= (y_val-range) && y1 <= (y_val+range)){
            for (int m=0; m<accu_w; m++){
               int r = int(round((x1-x_c)*cos(m*PI/accu_w) + (y1-y_c)*sin(m*PI/accu_w)));
               (adata[r][m])--;
            }
            deaccu.erase(deaccu.begin() + n);
         }
      } // end of deaccumulator loop

      int x0_end = endpoint[0][0], y0_end = endpoint[0][1], x1_end = endpoint[1][0], y1_end = endpoint[1][1];
      if ((x0_end==0 && y0_end==0) || (x1_end==0 && y1_end==0)) // don't add the (0,0) points 
         continue;
      vector<int> line = {x0_end, y0_end, x1_end, y1_end, rho, theta};
      outlines.push_back(line);

   } // end of point loop 
   // combine lines that are split 
   // int sameline1[2] = {0,0}; 
   // int sameline2[2] = {0,0};
   for (int i=0; i<outlines.size(); i++){
      bool same = false;
      for (int j=i+1; j<outlines.size() && same == false;j++){ 
         int xi_coords[2] = {outlines[i][0], outlines[i][2]}; int xj_coords[2] = {outlines[j][0], outlines[j][2]};
         int yi_coords[2] = {outlines[i][1], outlines[i][3]}; int yj_coords[2] = {outlines[j][1], outlines[j][3]};
         int rhoi = outlines[i][4], rhoj = outlines[j][4];
         int thetai = outlines[i][5], thetaj = outlines[j][5]; 

         int var = 100;
         int rho_var = 30;
         int theta_var = 20; 
         for (int k=0; k<2 && same == false; k++){
            for (int l=0; l<2 && same == false; l++){
               int counter = 0; 
               if ((xi_coords[k] < (xj_coords[l] + var)) && (xi_coords[k] > (xj_coords[l] - var)))
                  counter++;
               if ((yi_coords[k] < (yj_coords[l] + var)) && (yi_coords[k] > (yj_coords[l] - var)))
                  counter++ ;
               if ((rhoi < (rhoj + rho_var)) && (rhoi > (rhoj - rho_var)))
                  counter++; 
               if ((thetai < (thetaj + theta_var)) && (thetai > (thetaj - theta_var)))
                  counter++;
               if (counter >= 3){ // if at least three of the conditions are fulfilled 

                  if(k==0){
                     if(l==0){
                        outlines[j][2] = outlines[i][0];
                        outlines[j][3] = outlines[i][1];
                     }
                     else{
                        outlines[j][0] = outlines[i][0];
                        outlines[j][1] = outlines[i][1];
                     }
                  }
                  else{
                     if(l==0){
                        outlines[j][2] = outlines[i][2]; 
                        outlines[j][3] = outlines[i][3];                        
                     }
                     else{
                        outlines[j][0] = outlines[i][2];
                        outlines[j][1] = outlines[i][3]; 
                     }  
                  }
                  same = true;
                  (outlines.at(i)).clear();
                  //std::fill ((outlines.at(i)).begin(),(outlines.at(i)).end(),0);
               } 
            }
         }
      } // end of j loop 
   } // end of i loop 
/*    if (same){ // if there is line that's split into two 
      int a = sameline1[0], b = sameline1[1];
         // if b is 1, then we want to replace lines[a][2] and lines[a][3]
         // if b is 0, then we want to replace lines[a][0] and lines[a][1]

      int c = sameline2[0], d = sameline2[1]; 
         // if d is 0, then we want to replace the values with lines[c][2] and lines[c][3]
         // if d is 1, then we want to replace the values with lines[c][0] and lines[c][1]
      if(b==0){
         if(d==0){
            outlines[a][0] = outlines[c][2]; 
            outlines[a][1] = outlines[c][3];
         }
         else{
            outlines[a][0] = outlines[c][0];
            outlines[a][1] = outlines[c][1]; 
         }  
      }
      else{
         if(d==0){
            outlines[a][2] = outlines[c][2]; 
            outlines[a][3] = outlines[c][3];
         }
         else{
            outlines[a][2] = outlines[c][0];
            outlines[a][3] = outlines[c][1]; 
         }  
      }
      outlines.erase(outlines.begin()+c);
   }
 */
   //anode-cathode crossing cut OR min length cut 
   for (int i=0; i<outlines.size(); i++){
      if ((outlines.at(i)).empty())
         continue;
      int x0_end = outlines[i][0], y0_end = outlines[i][1], x1_end = outlines[i][2], y1_end = outlines[i][3];
      if (muon_length!=0){
         int y_diff = abs(y0_end-y1_end);
         if (y_diff > muon_length){
            lines.push_back(outlines.at(i));
         }
      }
      else{
         float length = distance(x0_end,y0_end,x1_end,y1_end);
         if (length > min_length){
            lines.push_back(outlines.at(i));
         }
      }
   }
   //free memory 
   data.clear(); deaccu.clear(); outlines.clear();
} // end of hough 

void ana_one(int threshold, int max_gap, int range, float min_length , int muon_length, int plane, int nentry){

   //--------------- open up hitdumpertree --------------// 
   TFile *infile = new TFile("hitdumper_tree.root");
   TTree *intree = (TTree*)infile->Get("/hitdumper/hitdumpertree");

   //declare tree variables
   int nhits;
   vector<int> *hit_tpc = 0;
   vector<int> *hit_plane = 0;
   vector<int> *hit_wire = 0; 
   vector<double> *hit_peakT = 0;

   int run;
   int subrun;
   int event;
   //get branches from tree
   intree->SetBranchAddress("nhits",&nhits);
   intree->SetBranchAddress("hit_tpc",&hit_tpc);
   intree->SetBranchAddress("hit_plane",&hit_plane);
   intree->SetBranchAddress("hit_wire",&hit_wire);
   intree->SetBranchAddress("hit_peakT",&hit_peakT);
   intree->SetBranchAddress("run",&run);
   intree->SetBranchAddress("subrun",&subrun);
   intree->SetBranchAddress("event",&event);

   intree->GetEntry(nentry);
   
   //-------------- create output tree --------------// 

   TFile *outfile = new TFile("stored_endpoints.root","RECREATE");

   TNtuple *points0 = new TNtuple("points0","hit points (tpc0)","hit_wire:hit_peakT");
   TNtuple *points1 = new TNtuple("points1","hit points (tpc1)","hit_wire:hit_peakT");
   TNtuple *endpoints0 = new TNtuple("end_points0","end points (tpc0)","end_x:end_y");
   TNtuple *endpoints1 = new TNtuple("end_points1","end points (tpc1)","end_x:end_y");

   // initialize vectors 
   vector<vector<int>> coords0;//(nhits,vector<int>(2)); // points will be removed from coords throughout selection; hits in tpc0 
   vector<vector<int>> coords1;//(nhits,vector<int>(2)); // hits in tpc1 

   coords0.clear(); coords1.clear();
   for (int i=0; i<nhits; i++){ // loop to fill points 
      int v0 = hit_wire->at(i), v1 = int(hit_peakT->at(i));
      if (v1 > 10 && v1 < 3490 && hit_plane->at(i)==plane){
         vector<int> v{v0,v1};
         if (hit_tpc->at(i) == 0){
            coords0.push_back(v);
            points0->Fill(v0, v1);
         }
         if (hit_tpc->at(i) == 1){
            coords1.push_back(v);
            points1->Fill(v0,v1);
         }
      } 
   } // end of coordinate loop

   vector<vector<int>> lines0; // not sure if i should worry about dynamically allocating the memory :| 
   vector<vector<int>> lines1; 
   lines0.clear(); lines1.clear();

   hough(coords0, threshold, max_gap, range, min_length, muon_length, nentry, lines0); 
   hough(coords1, threshold, max_gap, range, min_length, muon_length, nentry, lines1);

   // plotting
   for(int i=0; i<lines0.size(); i++){
      int wire0 = lines0[i][0], peakT0 = lines0[i][1], wire1 = lines0[i][2],  peakT1 = lines0[i][3]; 
      endpoints0->Fill(lines0[i][0],lines0[i][1]);
      endpoints0->Fill(lines0[i][2],lines0[i][3]);
      
      double theta_zx, dz = (wire1-wire0)*0.3, dx = (peakT1-peakT0)*0.5*0.16;
      double phi = atan(dx/dz)* 180 / PI;
      if (signbit(dz) == signbit(dx))
         theta_zx = abs(phi);
      else
         theta_zx = 180 - abs(phi);      
      std::cout << "tpc = 0:" << std::endl;
      std::cout << "endpoints: (" << lines0[i][0] << ", " << lines0[i][1] << "), (" << lines0[i][2] << ", " << lines0[i][3] << ") " << endl;
      std::cout << "theta_zx = " << theta_zx << endl;
      //          << "rho: " << lines0[i][4] << ", theta: " << lines0[i][5] << endl; 
   }
   for(int i=0; i<lines1.size(); i++){
      int wire0 = lines1[i][0], peakT0 = lines1[i][1], wire1 = lines1[i][2],  peakT1 = lines1[i][3]; 
      endpoints1->Fill(lines1[i][0],lines1[i][1]);
      endpoints1->Fill(lines1[i][2],lines1[i][3]);

      double theta_zx, dz = (wire1-wire0)*0.3, dx = (peakT1-peakT0)*0.5*0.16;
      double phi = atan(dx/dz)* 180 / PI;
      if (signbit(dz) == signbit(dx))
         theta_zx = abs(phi);
      else
         theta_zx = 180 - abs(phi);      
      std::cout << "tpc = 1:" << std::endl;
      std::cout << "endpoints: (" << lines1[i][0] << ", " << lines1[i][1] << "), (" << lines1[i][2] << ", " << lines1[i][3] << ") " << endl;
      std::cout << "theta_zx = " << theta_zx << endl;
      //          << "rho: " << lines0[i][4] << ", theta: " << lines0[i][5] << endl; 
   }
   //Convert run info to string
   string srun = to_string(run);
   string ssubrun = to_string(subrun);
   string sevent = to_string(event);
   string title = "Run: "+srun + " Subrun: "+ssubrun+" Event: "+sevent;
   //std::cout<<"****"<<run<<"*****"<<srun<<"*****"<<title<<endl;

   TCanvas *c1=new TCanvas("c1","c1",800,600);
   c1->cd();
   points0->SetMarkerStyle(7);
   points0->SetMarkerColor(46);
   points0->Draw("hit_peakT:hit_wire");
      
   TH2F *htemp = (TH2F*)gPad->GetPrimitive("htemp");
   TAxis *xaxis = htemp->GetXaxis();
   TAxis *yaxis = htemp->GetYaxis();
   htemp->SetTitle(title.c_str());
   yaxis->SetLimits(0, 3500);
   xaxis->SetLimits(0, 2000);

   points1->SetMarkerStyle(7);
   points1->SetMarkerColor(38);
   points1->Draw("hit_peakT:hit_wire","","SAME");

   endpoints0->SetMarkerStyle(8);
   endpoints0->SetMarkerSize(1);
   endpoints0->SetMarkerColor(2);
   endpoints0->Draw("end_y:end_x","", "SAME"); 

   endpoints1->SetMarkerStyle(8);
   endpoints1->SetMarkerSize(1);
   endpoints1->SetMarkerColor(4);
   endpoints1->Draw("end_y:end_x","", "SAME"); 

   TFile *tf = new TFile("histos.root","recreate");
   c1->Write();
   tf->Close();

   infile->Close();
   outfile->Close();


} // end of ana_one 

void ana(int threshold, int max_gap, int range, float min_length, int muon_length){
   // open up hitdumpertree
   TChain chain("hitdumper/hitdumpertree");
   chain.Add("hitdumper_sample3.root");
   //chain.Add("hitdumper_planes1.root");

   //declare input tree variables
   int nhits;
   vector<int> *hit_tpc = 0;
   vector<int> *hit_plane = 0;
   vector<int> *hit_wire = 0; 
   vector<double> *hit_peakT = 0;
   //get branches from tree
   chain.SetBranchAddress("nhits",&nhits);
   chain.SetBranchAddress("hit_tpc",&hit_tpc);
   chain.SetBranchAddress("hit_plane",&hit_plane);
   chain.SetBranchAddress("hit_wire",&hit_wire);
   chain.SetBranchAddress("hit_peakT",&hit_peakT);

   // create output tree
   int out_event;
   vector<int> end_wire0; 
   vector<int> end_wire1;
   vector<int> end_peakT0; 
   vector<int> end_peakT1;
   vector<int> end_tpc;
   vector<int> end_theta_zx; 

   TFile *outfile = new TFile("stored_endpoints.root","RECREATE");
   TTree *outtree = new TTree("endpointtree","tree that contains endpoint info");

   outtree->Branch("out_event",&out_event);
   outtree->Branch("end_wire0",&end_wire0);
   outtree->Branch("end_wire1",&end_wire1);
   outtree->Branch("end_peakT0",&end_peakT0);
   outtree->Branch("end_peakT1",&end_peakT1);
   outtree->Branch("end_tpc",&end_tpc);
   outtree->Branch("end_theta_zx",&end_theta_zx);

   //declare variables 
   vector<vector<int>> coords0;//(nhits,vector<int>(2)); // points will be removed from coords throughout selection; hits in tpc0 
   vector<vector<int>> coords1;//(nhits,vector<int>(2)); // hits in tpc1
   vector<vector<int>> lines0; // not sure if i should worry about dynamically allocating the memory :| 
   vector<vector<int>> lines1;
   int nlines0 = 0; 
   int nlines1 = 0; 

   int nevents = chain.GetEntries();

   for (int n=0; n<nevents; n++){ //start event loop //If y completely divides x, the result of the expression is 0.
       if (n%10 == 0 && muon_length==0){
         std::cout << "Processing event # " << n << " out of " << nevents << " events..." << endl;
      } 
      // load event 
      chain.GetEntry(n); 
      // initliaze variables 
      coords0.clear(); coords1.clear(); 
      lines0.clear(); lines1.clear();
      end_wire0.clear(); end_wire1.clear(); end_peakT0.clear(); end_peakT1.clear(); end_tpc.clear(); end_theta_zx.clear(); 


      // loop to fill points 
      for (int i=0; i<nhits; i++){
         int v0 = hit_wire->at(i), v1 = int(hit_peakT->at(i));
         if (v1 > 10 && v1 < 3490 && hit_plane->at(i)==2){
            vector<int> v{v0,v1};
            if (hit_tpc->at(i) == 0){
               coords0.push_back(v);
            }
            if (hit_tpc->at(i) == 1){
               coords1.push_back(v);
            }
         } 
      } // end of coordinate loop

      hough(coords0, threshold, max_gap, range, min_length, muon_length, n, lines0); 
      hough(coords1, threshold, max_gap, range, min_length, muon_length, n, lines1);

      out_event = n;
      if (muon_length!=0 && (lines0.empty() == false || lines1.empty() == false)){
          std::cout << "event #: " << out_event + 1 << endl;
      }
      for(int i=0; i<lines0.size(); i++){
         int wire0 = lines0[i][0], peakT0 = lines0[i][1], wire1 = lines0[i][2],  peakT1 = lines0[i][3]; 
         end_wire0.push_back(wire0);
         end_peakT0.push_back(peakT0);
         end_wire1.push_back(wire1);
         end_peakT1.push_back(peakT1);
         end_tpc.push_back(0);
         if (muon_length!=0){
            std::cout << "tpc = 0 " << endl;
            std::cout << "endpoints: (" << wire0 << ", " << peakT0 << "), (" << wire1 << ", " << peakT1 << ") " << endl; 
            double theta_zx, dz = (wire1-wire0)*0.3, dx = (peakT1-peakT0)*0.5*0.16;
            double phi = atan(dx/dz)* 180 / PI;
            if (signbit(dz) == signbit(dx))
               theta_zx = abs(phi);
            else
               theta_zx = 180 - abs(phi);
            std::cout << "theta_zx = " << theta_zx << endl;
            end_theta_zx.push_back(theta_zx);
            int t_i, t_j; 
            t_i = int((peakT0-500)*0.5);
            t_j = int((peakT1-500)*0.5);
            // int t0 = (t_i<t_j)? t_i:t_j; 
            // std::cout << "interaction times (in us): " << t0 << endl;
         }
      }

      for(int i=0; i<lines1.size(); i++){
         int wire0 = lines1[i][0], peakT0 = lines1[i][1], wire1 = lines1[i][2],  peakT1 = lines1[i][3]; 
         end_wire0.push_back(wire0);
         end_peakT0.push_back(peakT0);
         end_wire1.push_back(wire1);
         end_peakT1.push_back(peakT1);
         end_tpc.push_back(1);
         if (muon_length!=0){
            std::cout << "tpc = 1" << endl; 
            std::cout << "endpoints: (" << wire0 << ", " << peakT0 << "), (" << wire1 << ", " << peakT1 << ") " << endl; 
            double theta_zx, dz = (wire1-wire0)*0.3, dx = (peakT1-peakT0)*0.5*0.16;
            double phi = atan(dx/dz)* 180 / PI;
            if (signbit(dz) == signbit(dx))
               theta_zx = 180 - abs(phi);
            else
               theta_zx = abs(phi);
            std::cout << "theta_zx = " << theta_zx << endl;
            end_theta_zx.push_back(theta_zx);
            // int t_i, t_j; 
            // t_i = int((lines0[i][1]-500)*0.5);
            // t_j = int((lines0[i][3]-500)*0.5);
            // int t0 = (t_i<t_j)? t_i:t_j; 
            // std::cout << "interaction times (in us): " << t0 << endl;
         }
      }
      nlines0 += lines0.size(); 
      nlines1 += lines1.size(); 
      outtree->Fill();
   } //end of event loop
   std::cout << "number of muons detected in TPC==0: " << nlines0 << endl;
   std::cout << "number of muons detected in TPC==1: " << nlines1 << endl;
   std::cout << "total number of lines out of " << nevents << " events: " << nlines0+nlines1 << endl;
   outfile->Write();

   outfile->Close(); 
}
