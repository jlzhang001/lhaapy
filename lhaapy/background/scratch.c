
/* Notice:
 * 1. angular value in function which are included in Astro.h
 *    is in degree !!! ;
 * 2. global and other angular value are in rad;
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>

#include <slalib.h>
#include <slamac.h>
#include <drange.c>
#include <dh2e.c>
//#include <eqeqx.c>
//#include <nutc.c>
#include <dranrm.c>
#include <gmst.c>   

#include <Astro.h>  
#include <Astro.c>  
#include <Moon.c>  

#include <EvRec2.h>  
#include <EvRec2.C>  
#define EvRec2_h   

#define RA 166.114
#define DEC 38.209
#define NHIT_min 100
#define REC_QUAL_min 1 
#define ZENITH_max 40.

using namespace std;

int main( int argc, char **argv)
{

  if (argc <3) {
    printf("Augument error running: output.root filename.lis\n");
    exit(0);
  }

  int Nrec;
  double mjd;
  double azi, zen;
  double ha, dec;

  // exist;
  TFile *fout;
  TTree *tree1;
  fout = new TFile(argv[1],"recreate");
  tree1= new TTree("seawolf","seawolf");
  tree1-> Branch("mjd", &mjd, "mjd/D");
  tree1-> Branch("ha", &ha, "hour angle/D");
  tree1-> Branch("dec", &dec, "declination/D");

  ifstream fp;
  fp.open(argv[2]);
  char inputrootname[1024];
  int ifile= 0;
  while (fp>> inputrootname) {

    ifile++;
    cout<< "input root file"<< ifile<< " ="<< inputrootname<< endl;
    // exist;
    ifstream temp_file;
    temp_file.open(inputrootname);
    if (!temp_file) {
      cout<<"not exist the file"<< endl;
      continue;
    }

    TFile* inputroot= new TFile( inputrootname);
    // size;
    if (inputroot-> GetSize() <1000) {
      inputroot->Close(); cout<<"size err"<< endl; continue;
    }
    // zombie;
    if (inputroot-> IsZombie()) {
      inputroot->Close(); cout<<"zombie err"<< endl;continue;
    }

    TTree* tree2= (TTree*) inputroot-> Get("EvRec2");
    EvRec2 rec(tree2);
    Nrec= rec.fChain-> GetEntries();
    cout<< "Nrec="<< Nrec<< endl;

    //to read reconstruction information;
    for (int ii=0; ii<Nrec; ii++) {

      // cut condition;
      rec.fChain-> GetEntry(ii);
      if ( rec.nHit_ca< NHIT_min) continue;
      if ( rec.Rec_Qual< REC_QUAL_min) continue;

      //mjd= rec.Year + rec.Day* 1.0e-4 
      //   +rec.Second* 1.0e-8 +rec.mSecond* 1.0e-8;
      //zen_event= rec.Rec_Th;
      //azi_event= rec.Rec_Ph;
      mjd= (rec.Year+ 2000.- 1859)*365+ 44.+ int(( rec.Year+ 2000.- 1857)/ 4)+ rec.Day- 1+ (rec.Second+ rec.mSecond/ 1000000.)/ 86400.;
      zen= rec.Rec_Th_cn* DR2D;
      azi= rec.Rec_Ph_cn* DR2D;

      //zen= asin(1.042* sin(zen* DD2R))* DR2D;
      if (isnan(zen)) continue;
      if (zen <0 || zen >ZENITH_max) continue;

      azi= 270- azi- 18.5;
      // -pi, pi;
      azi= slaDrange(azi* DD2R)* DR2D;

      slaDh2e( azi*DD2R, (90-zen)*DD2R, tibet_la*DD2R, &ha, &dec);
      ha*= DR2D; dec*= DR2D;

      // test 1 // 
      //double sidereal_time= slaGmst(mjd)* DR2D;
      //double ra1=sidereal_time+ tibet_lo- ha;
      //ra1= slaDranrm(ra1*DD2R)* DR2D;
      //cout<< "method1: "<< ra1<< endl;
      //cout<< "  dec    "<< dec<< endl;
      //cout<< "         "<< "sidereal is "<< sidereal_time<< " ha is "<< ha<< endl;
      //cout<< "         "<< "mjd is "<< mjd<< endl;
      //double ra;
      //horizon_equator( mjd, zen, azi, &ra, &dec);
      //cout<< "      2: "<< ra<< endl;
      //cout<< "  dec    "<< dec<< endl;


      tree1-> Fill();
    }
    inputroot-> Close();
  }

  // 分列 //
  int num;
  double mjd_min, mjd_max;
  float mjd_divide[100]; 
  float mjd_divide_min, mjd_divide_max;
  num= tree1-> GetEntries();
  tree1-> GetEntry(0);
  mjd_min= mjd;
  mjd_divide_min= ((int)mjd_min) + ((int)((mjd_min-(int)mjd_min)/0.125))* 0.125; 
  tree1-> GetEntry(num-1);
  mjd_max= mjd;
  mjd_divide_max= ((int)mjd_max) + ((int)((mjd_max-(int)mjd_max)/0.125))* 0.125; 
  int i_mjd_divide=0;
  for (float ii= mjd_divide_min; ii<= mjd_divide_max; ii+=0.125) {
    mjd_divide[i_mjd_divide]= ii;
    i_mjd_divide++;
  }

  for (int ii=0; ii< i_mjd_divide; ii++) {
    cout<< "mjd_divide is "<< mjd_divide[ii]<< endl;
    char f_divide_name[1024];
    sprintf(f_divide_name,"%s.mjd%.3f.root",argv[1],mjd_divide[ii]);
    cout<< "f_divide name is "<< f_divide_name<< endl;
    TFile *f_divide= new TFile(f_divide_name,"recreate");
    TTree *t_divide= new TTree("seawolf","seawolf");
    TTree *t_divide_mjd= new TTree("mjd","mjd");
    t_divide-> Branch("mjd", &mjd, "mjd/D");
    t_divide-> Branch("ha", &ha, "hour angle/D");
    t_divide-> Branch("dec", &dec, "declination/D");
    t_divide_mjd-> Branch("mjd", &mjd, "mjd/D");
    for (int jj=0; jj<num; jj++) {
      tree1-> GetEntry(jj);
      if (mjd >=mjd_divide[ii] && mjd< mjd_divide[ii]+0.125) {
	t_divide-> Fill();
	t_divide_mjd-> Fill();
      }
    }
    f_divide-> cd();
    t_divide-> Write();
    t_divide_mjd-> Write();
    f_divide-> Close();
  }
    

  fp.close();
  fout-> cd();
  tree1-> Write();
  fout-> Close();
  cout<< "WeLlDoNe"<< endl;
}


// test 1 // to test out the fourmula from Local equatorial to equatorial;
