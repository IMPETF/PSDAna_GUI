#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TList.h"
#include "TCollection.h"
#include "TString.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TPaletteAxis.h"
#include "TProfile.h"
#include "TLegend.h"

Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");

   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0q");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}

Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

    printf("inf langaupro\n");
   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 100000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
         step = -step/10;

      p += step;
   }
   printf("%d\n",i);
   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;


   return (0);
}

TF1* langaus( TH1F *poHist,float& mpv,float& fwhm,float ped_mean,float ped_sigma)
{
    // Setting fit range and start values
    Double_t fr[2];
    Double_t sv[4],paramlow[4],paramhigh[4],fparam[4],fpe[4];
    Double_t chisqr;
    Int_t ndf;

    fr[0]=ped_mean+25*ped_sigma;
    fr[1]=fr[0]+2500;

    paramlow[0]=0.0;paramlow[1]=400.0;paramlow[2]=10;paramlow[3]=0.0;
    paramhigh[0]=200.0;paramhigh[1]=1500.0;paramhigh[2]=10000000.0;paramhigh[3]=200.0;

    sv[0]=30;
    sv[1]=poHist->GetBinCenter(poHist->GetMaximumBin());
    sv[2]=poHist->Integral(poHist->GetBin(fr[0]),poHist->GetBin(fr[1]));
    sv[3]=50;

    TF1 *fit = langaufit(poHist,fr,sv,paramlow,paramhigh,fparam,fpe,&chisqr,&ndf);

    mpv=fparam[1];
    fwhm=fparam[0];

    return fit;
}

TF1* pedfit(TH1F* poHist,float ped_mean,float ped_sigma,float& pedout_mean,float& pedout_sigma)
{
    //get pedestal from mips histogram
    TF1* fgaus=new TF1(Form("ped_%s",poHist->GetName()),"gaus",ped_mean-30.0,ped_mean+30.0);
    //fgaus->SetParameter(1,ped_mean);
    //fgaus->SetParameter(2,ped_sigma);
    poHist->Fit(Form("ped_%s",poHist->GetName()),"Rq0");
    pedout_mean=fgaus->GetParameter(1);
    pedout_sigma=fgaus->GetParameter(2);

    return fgaus;
}

/*********   This part is for PSD Nanjing testing. There are 4 FEEs  *************/

//---Method 1:convert each event block, synchronization between FEEs----------------------------------------------------------------------------
int convert_normalblock(char* data,Int_t *Xpos,Int_t *Ypos,Int_t *Xneg,Int_t *Yneg)
{
    //printf("in packet\n");
   // const Int_t FEE[4]={0x20,0x24,0x28,0x2c};
    const Int_t FEE[4]={0x20,0x28,0x24,0x2c};
    Int_t* tmp[4];
    tmp[0]=Xpos;
    tmp[1]=Ypos;
    tmp[2]=Xneg;
    tmp[3]=Yneg;
    const Int_t ch_num=90;//each FEE has 90 electronic channels,82 of them are connected to PMTs
    const Int_t header_length=2+6+8;
    const int data_length=186;
    const int block_length=190;
    unsigned char* burst=(unsigned char*)data;

    int read_length=0;
    int trig_couter[4];
    int init=0;
    if(burst[0]==0xE2 && burst[1]==0x25){
        for(int feeid=0;feeid<4;feeid++){
            init=header_length+feeid*block_length;
            if(burst[init]==0xeb && burst[init+1]==0x90){
                if(burst[init+3]==FEE[feeid]){
                    read_length=((burst[init+4]&0xFF)<<8) + (burst[init+5]&0xFF);
                    if(read_length != data_length){
                        printf("error! wrong data length %d != %d \n",read_length,data_length);
                        return -1;
                    }
                    for(int data_id=0;data_id<ch_num;data_id++){
                        tmp[feeid][data_id]=((burst[init+6+data_id*2]&0xFF)<<8) + (burst[init+6+data_id*2+1]&0xFF);
                        if(tmp[feeid][data_id]==0){
                            tmp[feeid][data_id] = -5;
                        }
                        else if(tmp[feeid][data_id] > 0x3FFF){
                            tmp[feeid][data_id] = 16400;
                        }
                    }
                    trig_couter[feeid] = ((burst[init+186]&0xF)<<8 )+ (burst[init+187]&0xFF);
                }
                 else
                {
                    printf("error! wrong sequence in a packet\n");
                    return -1;
                }
            }
            else{
                printf("error! incomplete block\n");
                return -1;
            }
        }
    }
    else{
        printf("error! wrong header(0x%x%x != 0xE225)\n",burst[0],burst[1]);
        return -1;
    }
    if((trig_couter[0] != trig_couter[1]) || (trig_couter[0] != trig_couter[2]) || (trig_couter[0] != trig_couter[3]) ){
        printf("error! inconsistent trigger count in a packet");
        return -1;
    }

    return 0;
}

int convert_normal(const Char_t* parentDir,const Char_t* infile,const Char_t* outDir,Char_t* outfile="raw.root")
{

    //const Int_t FEE[4]={0x20,0x24,0x28,0x2c};//+X,-X,+Y,-Y
    //const Int_t FEE[4]={0x20,0x28,0x24,0x2c};//+X,+Y,-X,-Y
    //const Int_t block_length=194;//55aaeb90****5aa5
    const Int_t block_length=190;//eb90****CRC
    const Int_t header_length=2+6+8;//E225+PrimaryHeader+SecondaryHeader
    const Int_t packet_length=header_length+block_length*4;

    Char_t infname[200],outfname[200];
    sprintf(infname,"%s/%s",parentDir,infile);
    sprintf(outfname,"%s/%s",outDir,outfile);

    TFile *f=new TFile(outfname,"RECREATE");
    if(f->IsZombie()){
        printf("this file will not been recreated!\n");
        delete f;
        return -1;
    }
    else if(f->IsOpen()){
        ifstream in;
        in.open(infname,std::ios_base::in | std::ios_base::binary);
        if(!in.is_open()){
            std::cout << "can't open"<< infname << std::endl;
        }
//-----------------------------------------------------------------------
        TH1F *hxpos[90];
        TH1F *hxneg[90];
        TH1F *hypos[90];
        TH1F *hyneg[90];
        for(int i=0;i<90;i++){
            hxpos[i]=new TH1F(Form("xpos_%d",i+1),Form("xpos_%d",i+1),4000,0,4000);
            hxneg[i]=new TH1F(Form("xneg_%d",i+1),Form("xneg_%d",i+1),4000,0,4000);
            hypos[i]=new TH1F(Form("ypos_%d",i+1),Form("ypos_%d",i+1),4000,0,4000);
            hyneg[i]=new TH1F(Form("yneg_%d",i+1),Form("yneg_%d",i+1),4000,0,4000);
        }

        TTree *tree = new TTree("PSD","PSD Testing event mode");
        Int_t Xpos[90],Ypos[90],Xneg[90],Yneg[90];
        tree->Branch("xpos",Xpos,"xpos[90]/I");
        tree->Branch("ypos",Ypos,"ypos[90]/I");
        tree->Branch("xneg",Xneg,"xneg[90]/I");
        tree->Branch("yneg",Yneg,"yneg[90]/I");
//----------------------------------------------------------------------------------------
    unsigned int event_num=0;
    char burst_block[packet_length];
    while(!in.eof()){
        //printf("test");
        in.read(burst_block,packet_length);
        if(!in.eof()){
            int error_flag=convert_normalblock(burst_block,Xpos,Ypos,Xneg,Yneg);
            if(error_flag == -1){
                break;
            }
            else{
                event_num++;
                for(int i=0;i<90;i++){
                    hxpos[i]->Fill(Xpos[i]);
                    hypos[i]->Fill(Ypos[i]);
                    hxneg[i]->Fill(Xneg[i]);
                    hyneg[i]->Fill(Yneg[i]);
                }
                tree->Fill();
            }
        }
        else{
            break;
        }
    }
    std::cout<< event_num <<" events converted totally!"<<std::endl;

    in.close();
    f->Write(0,TObject::kOverwrite);//overwrite AutoSave keys
    f->Close();
    delete f;
  }
  else{
    printf("error: %s can not be created!\n",outfile);
    delete f;
    return -1;
  }
    return 0;
}

//---compressed mode--------------------------------------
int convert_compressedblock(char* data,Int_t *Xpos,Int_t *Ypos,Int_t *Xneg,Int_t *Yneg)
{
    //printf("in packet\n");
   // const Int_t FEE[4]={0x20,0x24,0x28,0x2c};
    const Int_t FEE[4]={0x60,0x68,0x64,0x6c};
    Int_t* tmp[4];
    tmp[0]=Xpos;
    tmp[1]=Ypos;
    tmp[2]=Xneg;
    tmp[3]=Yneg;
    const Int_t ch_num=90;//each FEE has 90 electronic channels,82 of them are connected to PMTs
    const Int_t header_length=2+6+8;
    const int data_length=186;
    const int block_length=190;
    unsigned char* burst=(unsigned char*)data;

    int read_length=0;
    int trig_couter[4];
    int init=0;
    if(burst[0]==0xE2 && burst[1]==0x25){
        for(int feeid=0;feeid<4;feeid++){
            init=header_length+feeid*block_length;
            if(burst[init]==0xeb && burst[init+1]==0x90){
                if(burst[init+3]==FEE[feeid]){
                    read_length=((burst[init+4]&0xFF)<<8) + (burst[init+5]&0xFF);
                    if(read_length != data_length){
                        printf("error! wrong data length %d != %d \n",read_length,data_length);
                        return -1;
                    }
                    for(int data_id=0;data_id<ch_num;data_id++){
                        tmp[feeid][data_id]=((burst[init+6+data_id*2]&0xFF)<<8) + (burst[init+6+data_id*2+1]&0xFF);
                        if(tmp[feeid][data_id]==0){
                            tmp[feeid][data_id] = -5;
                        }
                        else if(tmp[feeid][data_id] > 0x3FFF){
                            tmp[feeid][data_id] = 16400;
                        }
                    }
                    trig_couter[feeid] = ((burst[init+186]&0xF)<<8 )+ (burst[init+187]&0xFF);
                }
                 else
                {
                    printf("error! wrong sequence in a packet\n");
                    return -1;
                }
            }
            else{
                printf("error! incomplete block\n");
                return -1;
            }
        }
    }
    else{
        printf("error! wrong header(0x%x%x != 0xE225)\n",burst[0],burst[1]);
        return -1;
    }
    if((trig_couter[0] != trig_couter[1]) || (trig_couter[0] != trig_couter[2]) || (trig_couter[0] != trig_couter[3]) ){
        printf("error! inconsistent trigger count in a packet");
        return -1;
    }

    return 0;
}

int convert_compressed(const Char_t* parentDir,const Char_t* infile,const Char_t* outDir,Char_t* outfile="raw.root")
{
    //const Int_t FEE[4]={0x20,0x24,0x28,0x2c};//+X,-X,+Y,-Y
    //const Int_t FEE[4]={0x20,0x28,0x24,0x2c};//+X,+Y,-X,-Y
    //const Int_t block_length=194;//55aaeb90****5aa5
    const Int_t block_length=190;//eb90****CRC
    const Int_t header_length=2+6+8;//E225+PrimaryHeader+SecondaryHeader
    const Int_t packet_length=header_length+block_length*4;

    Char_t infname[200],outfname[200];
    sprintf(infname,"%s/%s",parentDir,infile);
    sprintf(outfname,"%s/%s",outDir,outfile);

    TFile *f=new TFile(outfname,"RECREATE");
    if(f->IsZombie()){
        printf("this file will not been recreated!\n");
        delete f;
        return -1;
    }
    else if(f->IsOpen()){
        ifstream in;
        in.open(infname,std::ios_base::in | std::ios_base::binary);
        if(!in.is_open()){
            std::cout << "can't open"<< infname << std::endl;
        }
//-----------------------------------------------------------------------
        TH1F *hxpos[90];
        TH1F *hxneg[90];
        TH1F *hypos[90];
        TH1F *hyneg[90];
        for(int i=0;i<90;i++){
            hxpos[i]=new TH1F(Form("xpos_%d",i+1),Form("xpos_%d",i+1),4000,0,4000);
            hxneg[i]=new TH1F(Form("xneg_%d",i+1),Form("xneg_%d",i+1),4000,0,4000);
            hypos[i]=new TH1F(Form("ypos_%d",i+1),Form("ypos_%d",i+1),4000,0,4000);
            hyneg[i]=new TH1F(Form("yneg_%d",i+1),Form("yneg_%d",i+1),4000,0,4000);
        }

        TTree *tree = new TTree("PSD","PSD Testing event mode");
        Int_t Xpos[90],Ypos[90],Xneg[90],Yneg[90];
        tree->Branch("xpos",Xpos,"xpos[90]/I");
        tree->Branch("ypos",Ypos,"ypos[90]/I");
        tree->Branch("xneg",Xneg,"xneg[90]/I");
        tree->Branch("yneg",Yneg,"yneg[90]/I");
//----------------------------------------------------------------------------------------
    unsigned int event_num=0;
    char burst_block[packet_length];
    while(!in.eof()){
        //printf("test");
        in.read(burst_block,packet_length);
        if(!in.eof()){
            int error_flag=convert_compressedblock(burst_block,Xpos,Ypos,Xneg,Yneg);
            if(error_flag == -1){
                break;
            }
            else{
                event_num++;
                for(int i=0;i<90;i++){
                    hxpos[i]->Fill(Xpos[i]);
                    hypos[i]->Fill(Ypos[i]);
                    hxneg[i]->Fill(Xneg[i]);
                    hyneg[i]->Fill(Yneg[i]);
                }
                tree->Fill();
            }
        }
        else{
            break;
        }
    }
    std::cout<< event_num <<" events converted totally!"<<std::endl;

    in.close();
    f->Write(0,TObject::kOverwrite);//overwrite AutoSave keys
    f->Close();
    delete f;
  }
  else{
    printf("error: %s can not be created!\n",outfile);
    delete f;
    return -1;
  }
    return 0;
}

//---calibration mode-------------------------------------
int convert_calibblock(char* data,Int_t *Xpos,Int_t *Ypos,Int_t *Xneg,Int_t *Yneg)
{
    //printf("in packet\n");
   // const Int_t FEE[4]={0x20,0x24,0x28,0x2c};
    const Int_t FEE[4]={0xa0,0xa8,0xa4,0xac};
    Int_t* tmp[4];
    tmp[0]=Xpos;
    tmp[1]=Ypos;
    tmp[2]=Xneg;
    tmp[3]=Yneg;
    const Int_t ch_num=90;//each FEE has 90 electronic channels,82 of them are connected to PMTs
    const Int_t header_length=2+6+8;
    const int data_length=186;
    const int block_length=190;
    unsigned char* burst=(unsigned char*)data;

    int read_length=0;
    int trig_couter[4];
    int init=0;
    if(burst[0]==0xE2 && burst[1]==0x25){
        for(int feeid=0;feeid<4;feeid++){
            init=header_length+feeid*block_length;
            if(burst[init]==0xeb && burst[init+1]==0x90){
                if(burst[init+3]==FEE[feeid]){
                    read_length=((burst[init+4]&0xFF)<<8) + (burst[init+5]&0xFF);
                    if(read_length != data_length){
                        printf("error! wrong data length %d != %d \n",read_length,data_length);
                        return -1;
                    }
                    for(int data_id=0;data_id<ch_num;data_id++){
                        tmp[feeid][data_id]=((burst[init+6+data_id*2]&0xFF)<<8) + (burst[init+6+data_id*2+1]&0xFF);
                        if(tmp[feeid][data_id]==0){
                            tmp[feeid][data_id] = -5;
                        }
                        else if(tmp[feeid][data_id] > 0x3FFF){
                            tmp[feeid][data_id] = 16400;
                        }
                    }
                    trig_couter[feeid] = ((burst[init+186]&0xF)<<8 )+ (burst[init+187]&0xFF);
                }
                 else
                {
                    //printf("error! wrong sequence in a packet\n");
                    //break;
                    return -2;
                }
            }
            else{
                printf("error! incomplete block\n");
                return -1;
            }
        }
    }
    else{
        printf("error! wrong header(0x%x%x != 0xE225)\n",burst[0],burst[1]);
        return -1;
    }
    if((trig_couter[0] != trig_couter[1]) || (trig_couter[0] != trig_couter[2]) || (trig_couter[0] != trig_couter[3]) ){
        printf("error! inconsistent trigger count in a packet");
        return -1;
    }

    return 0;
}

int convert_calib(const Char_t* parentDir,const Char_t* infile,const Char_t* outDir,Char_t* outfile="calib_raw.root")
{

    //const Int_t FEE[4]={0x20,0x24,0x28,0x2c};//+X,-X,+Y,-Y
    //const Int_t FEE[4]={0x20,0x28,0x24,0x2c};//+X,+Y,-X,-Y
    //const Int_t block_length=194;//55aaeb90****5aa5
    const Int_t block_length=190;//eb90****CRC
    const Int_t header_length=2+6+8;//E225+PrimaryHeader+SecondaryHeader
    const Int_t packet_length=header_length+block_length*4;

    Char_t infname[200],outfname[200];
    sprintf(infname,"%s/%s",parentDir,infile);
    sprintf(outfname,"%s/%s",outDir,outfile);

    TFile *f=new TFile(outfname,"RECREATE");
    if(f->IsZombie()){
        printf("this file will not been recreated!\n");
        delete f;
        return -1;
    }
    else if(f->IsOpen()){
        ifstream in;
        in.open(infname,std::ios_base::in | std::ios_base::binary);
        if(!in.is_open()){
            std::cout << "can't open"<< infname << std::endl;
        }
//-----------------------------------------------------------------------
        TH1F *hxpos[90];
        TH1F *hxneg[90];
        TH1F *hypos[90];
        TH1F *hyneg[90];
        for(int i=0;i<90;i++){
            hxpos[i]=new TH1F(Form("xpos_%d",i+1),Form("xpos_%d",i+1),17000,0.5,17000.5);
            hxneg[i]=new TH1F(Form("xneg_%d",i+1),Form("xneg_%d",i+1),17000,0.5,17000.5);
            hypos[i]=new TH1F(Form("ypos_%d",i+1),Form("ypos_%d",i+1),17000,0.5,17000.5);
            hyneg[i]=new TH1F(Form("yneg_%d",i+1),Form("yneg_%d",i+1),17000,0.5,17000.5);
        }

        TTree *tree = new TTree("PSD_calib","PSD Electronic Calibration");
        Int_t Xpos[90],Ypos[90],Xneg[90],Yneg[90];
        tree->Branch("xpos",Xpos,"xpos[90]/I");
        tree->Branch("ypos",Ypos,"ypos[90]/I");
        tree->Branch("xneg",Xneg,"xneg[90]/I");
        tree->Branch("yneg",Yneg,"yneg[90]/I");
//----------------------------------------------------------------------------------------
    unsigned int event_num=0;
    unsigned int packet_num=0;
    char burst_block[packet_length];
    while(!in.eof()){
        //printf("test");
        in.read(burst_block,packet_length);
        if(!in.eof()){
            packet_num++;
            int error_flag=convert_calibblock(burst_block,Xpos,Ypos,Xneg,Yneg);
            if(error_flag == -1){
                printf("error in the middle of the file(%d packets have been successfully decoded)\n",packet_num-1);
                break;
            }
            else{
                if(error_flag == 0){
                    event_num++;
                    for(int i=0;i<90;i++){
                        hxpos[i]->Fill(Xpos[i]);
                        hypos[i]->Fill(Ypos[i]);
                        hxneg[i]->Fill(Xneg[i]);
                        hyneg[i]->Fill(Yneg[i]);
                    }
                    tree->Fill();
                }
            }
        }
        else{
            printf("success!:end of file reached \n");
            break;
        }
    }
    std::cout<< event_num <<" events converted totally!"<<std::endl;

    in.close();
    f->Write(0,TObject::kOverwrite);//overwrite AutoSave keys
    f->Close();
    delete f;
  }
  else{
    printf("error: %s can not be created!\n",outfile);
    delete f;
    return -1;
  }
    return 0;
}

//---all mode----------------------------------------------
/*
int convert_event(const Char_t* parentDir,const Char_t* infile,const Char_t* outDir,Char_t* outfile="raw.root")
{
    const Int_t block_length=190;//eb90****CRC
    const Int_t header_length=2+6+8;//E225+PrimaryHeader+SecondaryHeader
    const Int_t packet_length=header_length+block_length*4;

    Char_t infname[200],outfname[200];
    sprintf(infname,"%s/%s",parentDir,infile);
    sprintf(outfname,"%s/%s",outDir,outfile);

    TFile *f=new TFile(outfname,"RECREATE");
    if(f->IsZombie()){
        printf("this file will not been recreated!\n");
        delete f;
        return -1;
    }
    else if(f->IsOpen()){
        ifstream in;
        in.open(infname,std::ios_base::in | std::ios_base::binary);
        if(!in.is_open()){
            std::cout << "can't open"<< infname << std::endl;
        }
//-----------------------------------------------------------------------
        TH1F *hxpos[90];
        TH1F *hxneg[90];
        TH1F *hypos[90];
        TH1F *hyneg[90];
        for(int i=0;i<90;i++){
            hxpos[i]=new TH1F(Form("xpos_%d",i+1),Form("xpos_%d",i+1),4000,0,4000);
            hxneg[i]=new TH1F(Form("xneg_%d",i+1),Form("xneg_%d",i+1),4000,0,4000);
            hypos[i]=new TH1F(Form("ypos_%d",i+1),Form("ypos_%d",i+1),4000,0,4000);
            hyneg[i]=new TH1F(Form("yneg_%d",i+1),Form("yneg_%d",i+1),4000,0,4000);
        }

        TTree *tree = new TTree("PSD","PSD Testing event mode");
        Int_t Xpos[90],Ypos[90],Xneg[90],Yneg[90];
        tree->Branch("xpos",Xpos,"xpos[90]/I");
        tree->Branch("ypos",Ypos,"ypos[90]/I");
        tree->Branch("xneg",Xneg,"xneg[90]/I");
        tree->Branch("yneg",Yneg,"yneg[90]/I");
//----------------------------------------------------------------------------------------
    unsigned int event_num=0;
    unsigned char buffer;
    while(!in.fail()){
        if((((unsigned char)in.get())==0xe2) && !in.fail()){
            if((((unsigned char)in.get())==0x25) && !in.fail()){
                if((((unsigned char)in.get())==0x08) && !in.fail()){
                    if((((unsigned char)in.get())==0x13) && !in.fail()){
                        if((((unsigned char)in.get())==0xeb) && !in.fail()){
                            if((((unsigned char)in.get())==0x90) && !in.fail())
                        in.ignore(header_length);
                        if(!in.fail()){
                            buffer=in.get();
                        }
                        }

                    }
                }
            }
        }
    }
    std::cout<< event_num <<" events converted totally!"<<std::endl;

    in.close();
    f->Write(0,TObject::kOverwrite);//overwrite AutoSave keys
    f->Close();
    delete f;
  }
  else{
    printf("error: %s can not be created!\n",outfile);
    delete f;
    return -1;
  }

    return 0;
}
*/
//---Method 2: a much more robust decoding method---------------------------------------------------------
//---a. identify each packet using packet_header(0xE2 25 08 13)
//---b. insensitive to the order of FEEs in a single block
//---c. insensitive to the type_id of FEEs in a single block
//---d. decode all the info from the packet,including packet_id,time_code,trigger_id,trigger_status,science data
//---e. require the same trigger_id and type_id of all the FEEs in the same single packet
//---f. the following three functions are the basic utility functions used in the decoding process.
//---   i.e. convert_feeblock,convert_psd_packet_scidata,find_nextheader

int convert_feeblock(char* data,int remaining_length,int& expected_length,int& type_id,Int_t& fee_id,Int_t& trigger_counter,Int_t& trigger_status,UShort_t& status_code,Int_t *buffer)
{
    unsigned char* block=(unsigned char*)data;
    int channel_num,ch_id;
    if(block[0]==0xeb && block[1]==0x90){
        status_code = block[2];
        fee_id = block[3]&0xF;
        type_id = block[3]&0xF0;
        expected_length = (block[4]<<8)+block[5]+4;//data_length+header
        if(expected_length > remaining_length){
            printf("error!less bytes found than expected\n");
            return -2;
        }
        switch(type_id)
        {
        case 0x20://normal
            if(expected_length!=190){
                printf("error! wrong data_length in a normal block\n");
                return -1;
            }
            for(int i=0;i<90;i++){
                buffer[i]=(block[6+2*i]<<8)+block[6+2*i+1];
                if(buffer[i]>0x3FFF){
                    buffer[i]=16400;
                }
            }
            trigger_counter = ((block[186]&0xF)<<8) + (block[187]);
            trigger_status =block[186]>>4;
            break;
        case 0x60://compressed
            channel_num=(expected_length-10)/3;
            for(int i=0;i<channel_num;i++){
                ch_id=block[6+3*i];
                buffer[ch_id]=(block[6+3*i+1]<<8)+block[6+3*i+2];
                if(buffer[ch_id]>0x3FFF){
                    buffer[ch_id]=16400;
                }
            }
            trigger_counter = ((block[expected_length-4]&0xF)<<8) + block[expected_length-3];
            trigger_status = block[expected_length-4]>>4;
            break;
        case 0xa0://calibration
            if(expected_length!=190){
                printf("error! wrong data_length in a calibration block\n");
                return -1;
            }
            for(int i=0;i<90;i++){
                buffer[i]=(block[6+2*i]<<8)+block[6+2*i+1];
                if(buffer[i]>0x3FFF){
                    buffer[i]=16400;
                }
            }
            trigger_counter = ((block[186]&0xF)<<8) + (block[187]);
            trigger_status =block[186]>>4;
            break;
        default:
            printf("error! not a valid type_id(FEE_typeid = 0x%x)\n",type_id);
            return -3;
        }
    }
    else{
        printf("error! This is not a valid block header(0xEB90 != 0x%x%x)\n",block[0],block[1]);
        return -1;
    }

    return 0;
}

int convert_psd_packet_scidata(bool minimum_flag,char* data,int total_length,Int_t& type_id,Int_t* actual_trigger_id,Int_t* actual_trigger_status,UShort_t* actual_status_code,Int_t *Xpos,Int_t *Ypos,Int_t *Xneg,Int_t *Yneg)
{
    const int FEEID[4]={0x0,0x8,0x4,0xc};//+X,+Y,-X,-Y
    const int DATATYPE[4]={0x20,0x60,0xa0};//type_id need to be one of them:normal,compressed,calibration

    //default value,if the FEE do not exist in this packet
    std::fill_n(Xpos,90,-5);
    std::fill_n(Ypos,90,-5);
    std::fill_n(Xneg,90,-5);
    std::fill_n(Yneg,90,-5);
    std::fill_n(actual_trigger_id,4,-5);
    std::fill_n(actual_trigger_status,4,-5);
    std::fill_n(actual_status_code,4,-5);
    //--------------
    int tmp_trigger_id;
    int tmp_trigger_status;
    int actual_typeid[4];
    int tmp_feeid,tmp_typeid;
    int buffer[90];
    int remaining_length,expected_length;
    UShort_t tmp_status_code;
    int fee_num=0;

    int flag;
    char* next_data;
    remaining_length=total_length;
    while(remaining_length > 0){
        std::fill_n(buffer,90,-5);//if no ADC value decoded,the default value of each channel will be -5
        next_data=data+total_length-remaining_length;
        flag=convert_feeblock(next_data,remaining_length,expected_length,tmp_typeid,tmp_feeid,tmp_trigger_id,tmp_trigger_status,tmp_status_code,buffer);
        if(flag<0){
            return -1;
        }
        else if(flag==0){
            fee_num++;
            remaining_length-=expected_length;
            actual_typeid[fee_num-1]=tmp_typeid;
            switch(tmp_feeid)
            {
            case 0x0:
                for(int i=0;i<90;i++){
                    Xpos[i]=buffer[i];
                }
                actual_trigger_id[0]=tmp_trigger_id;
                actual_trigger_status[0]=tmp_trigger_status;
                actual_status_code[0]=tmp_status_code;
                break;
            case 0x8:
                for(int i=0;i<90;i++){
                    Ypos[i]=buffer[i];
                }
                actual_trigger_id[1]=tmp_trigger_id;
                actual_trigger_status[1]=tmp_trigger_status;
                actual_status_code[1]=tmp_status_code;
                break;
            case 0x4:
                for(int i=0;i<90;i++){
                    Xneg[i]=buffer[i];
                }
                actual_trigger_id[2]=tmp_trigger_id;
                actual_trigger_status[2]=tmp_trigger_status;
                actual_status_code[2]=tmp_status_code;
                break;
            case 0xc:
                for(int i=0;i<90;i++){
                    Yneg[i]=buffer[i];
                }
                actual_trigger_id[3]=tmp_trigger_id;
                actual_trigger_status[3]=tmp_trigger_status;
                actual_status_code[3]=tmp_status_code;
                break;
            default:
                return -1;
            }
        }
    }

    if(fee_num == 0){
        printf("error! this patch don't contain any PSD_FEE data!\n");
        return -1;
    }
    else{
        if(!minimum_flag){
            for(int i=1;i<fee_num;i++){
                if(actual_typeid[0] != actual_typeid[i]){
                    printf("error! unmatched type_id in a single packet\n");
                    return -1;
                }
            }
            type_id=actual_typeid[0];
        }
        else{
            for(int i=0;i<fee_num;i++){
                if((actual_typeid[i])!=0x20 && (actual_typeid[i]!=0x60)){
                    printf("error! unidetified type_id(0x%x) in this packet\n",actual_typeid[i]);
                    return -1;
                }
            }
            type_id=0x30;//minimum_mode
        }
    }

    return 0;
}

int check_triggerid(Int_t* trigger_id)
{
    int tmp_triggerid=-5;
    for(int i=0;i<4;i++){
        if(trigger_id[i] != -5){
            if(tmp_triggerid == -5){
                tmp_triggerid=trigger_id[i];
            }
            else{
                if(tmp_triggerid != trigger_id[i]){
                    return -1;
                }
            }
        }
    }
    return 0;
}

int check_triggerstatus(Int_t* trigger_status,unsigned int packet_num=-1,bool verbose=false)
{
    char fee_id[][10]={"Xpos","Ypos","Xneg","Yneg"};
    int flag=0;

    for(int i=0;i<4;i++){
        if(trigger_status[i] != -5){
            if(trigger_status[i] == 0x2){
                if(verbose){
                    if(!flag){
                        printf("packet %u:\n",packet_num);
                        flag=1;
                    }
                    printf("\twarning: FEE_%s trigger status error in this packet\n",fee_id[i]);
                }
            }
        }
    }
    return 0;
}

int check_statuscode(UShort_t* status_code,unsigned int packet_num=-1,bool verbose=false)
{
    char fee_id[][10]={"Xpos","Ypos","Xneg","Yneg"};
    int flag=0;

    for(int i=0;i<4;i++){
        if(status_code[i] != -5){
            if(status_code[i] != 0xFF){
                if(verbose){
                    if(!flag){
                        printf("packet %u:\n",packet_num);
                        flag=1;
                    }

                    if(((status_code[i]>>7)&0x1) == 0){
                        if(verbose){
                            printf("\twarning: FEE_%s was cut power in this packet\n",fee_id[i]);
                        }
                    }
                    if(((status_code[i]>>6)&0x1) == 0){
                        if(verbose){
                            printf("\twarning: FEE_%s high_thresh error in this packet\n",fee_id[i]);
                        }
                    }
                    if(((status_code[i]>>5)&0x1) == 0){
                        if(verbose){
                            printf("\twarning: FEE_%s low_thresh error in this packet\n",fee_id[i]);
                        }
                    }
                    if(((status_code[i]>>3)&0x1) == 0){
                        if(verbose){
                            printf("\twarning: FEE_%s ADC_976 busyn error in this packet\n",fee_id[i]);
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int get_decodedfee(Int_t* trigger_id,Int_t* flag)
{
    int fee_num=0;
    std::fill_n(flag,4,0);
    for(int i=0;i<4;i++){
        if(trigger_id[i] > -1){
            flag[i]=1;
            fee_num++;
        }
    }

    return fee_num;
}
    //return value: 0:a complete packet has been found,1:the last packet will be either complete or not   2:the last packet is definitely incomplete,-1:error occured,3:no byte read,eof reached
int find_nextheader(std::ifstream& in,int& packet_length)
{
    const char HEADER[4]={0xe2,0x25,0x08,0x13};
    std::streampos init_pos;
    init_pos=in.tellg();
    packet_length=0;

    int count;
    char tmp[4];
    std::fill_n(tmp,4,0);
    in.read(tmp,4);
    count=in.gcount();
    if(in.eof()){
        if(count){
            in.clear();
            if(!memcmp(HEADER,tmp,count)){
                packet_length=count;
                printf("the last packet has incomplete header,no content\n");
                in.seekg(init_pos);
                return 2;
            }
            else{
                packet_length=count;
                printf("the last packet has wrong header,no content\n");
                in.seekg(init_pos);
                return 2;
            }
        }
        else{
            printf("end of file reached\n");
            return 3;
        }
    }
    else {
        if(!memcmp(HEADER,tmp,4)){
            packet_length+=4;
        }
        else{
            in.seekg(-3,in.cur);
            packet_length+=1;
        }
    }

    while(in.read(tmp,4)){
        if(!memcmp(HEADER,tmp,4)){
            in.seekg(init_pos);
            return 0;
        }
        else{
            in.seekg(-3,in.cur);
            packet_length++;
        }
    }

    if(in.eof()){
        count=in.gcount();
        in.clear();
        if(count<3){
            packet_length+=count;
            in.seekg(init_pos);
            printf("the last packet has right header,but the content is definitely incomplete\n");
            return 2;
        }
        else{
            if(!memcmp(HEADER,tmp,3)){
                in.seekg(init_pos);
                return 0;
            }
            else if(!memcmp(HEADER,tmp+1,2)){
                packet_length+=1;
                in.seekg(init_pos);
                return 0;
            }
            else if(!memcmp(HEADER,tmp+2,1)){
                packet_length+=2;
                in.seekg(init_pos);
                return 0;
            }
            else{
                packet_length+=3;
                in.seekg(init_pos);
                printf("the last packet has right header,but the content maybe incomplete(can't find next header)\n");
                return 1;
            }
        }
    }
    else{
        printf("fatal error!an error occured in reading the file:");
        //in.seekg(init_pos);
        //exit(EXIT_FAILURE);
         return -1;
    }

}

int find_nextheader_u(std::ifstream& in,int& packet_length,const char* HEADER,int header_length)
{
    std::streampos init_pos;
    init_pos=in.tellg();
    packet_length=0;

    int count;
    char tmp[4];
    std::fill_n(tmp,4,0);
    in.read(tmp,header_length);
    count=in.gcount();
    if(in.eof()){
        if(count){
            in.clear();
            if(!memcmp(HEADER,tmp,count)){
                packet_length=count;
                printf("the last packet has incomplete header,no content\n");
                in.seekg(init_pos);
                return 2;
            }
            else{
                packet_length=count;
                printf("the last packet has wrong header,no content\n");
                in.seekg(init_pos);
                return 2;
            }
        }
        else{
            printf("end of file reached\n");
            return 3;
        }
    }
    else {
        if(!memcmp(HEADER,tmp,header_length)){
            packet_length+=header_length;
        }
        else{
            in.seekg(-header_length+1,in.cur);
            packet_length+=1;
        }
    }

    while(in.read(tmp,header_length)){
        if(!memcmp(HEADER,tmp,header_length)){
            in.seekg(init_pos);
            return 0;
        }
        else{
            in.seekg(-header_length+1,in.cur);
            packet_length++;
        }
    }

    if(in.eof()){
        count=in.gcount();
        in.clear();
        if(count<(header_length-1)){
            packet_length+=count;
            in.seekg(init_pos);
            printf("the last packet has right header,but the content is definitely incomplete\n");
            return 2;
        }
        else{
            for(int i=0;i<(header_length-1);i++){
                if(!memcmp(HEADER,tmp+i,header_length-1-i)){
                    packet_length+=i;
                    in.seekg(init_pos);
                    return 0;
                }
            }

            packet_length+=(header_length-1);
            in.seekg(init_pos);
            printf("the last packet has right header,but the content maybe incomplete(can't find next header)\n");
            return 1;
        }
    }
    else{
        printf("fatal error!an error occured in reading the file:");
        return -1;
    }

}

int find_nextheader_uu(std::ifstream &in, int &packet_length, const char *HEADER, int header_length)
{
    std::streampos init_pos;
    init_pos=in.tellg();
    packet_length=0;

    int count;
    char tmp[4];
    std::fill_n(tmp,4,0);
    in.read(tmp,header_length);
    count=in.gcount();
    if(in.eof()){
        if(count){
            in.clear();
            if(!memcmp(HEADER,tmp,count)){
                packet_length=count;
                printf("the last packet has incomplete header,no content\n");
                in.seekg(init_pos);
                return 2;
            }
            else{
                packet_length=count;
                printf("the last packet has wrong header,no content\n");
                in.seekg(init_pos);
                return 2;
            }
        }
        else{
            printf("end of file reached\n");
            return 3;
        }
    }
    else {
        if(!memcmp(HEADER,tmp,header_length)){
            packet_length+=header_length;
        }
        else{
            in.seekg(-header_length+1,in.cur);
            packet_length+=1;
        }
    }

    while(in.read(tmp,header_length)){
        if(!memcmp(HEADER,tmp,header_length)){
            in.seekg(init_pos);
            return 0;
        }
        else{
            in.seekg(-header_length+1,in.cur);
            packet_length++;
        }
    }

    if(in.eof()){
        count=in.gcount();
        in.clear();
        if(count<(header_length-1)){
            packet_length+=count;
            in.seekg(init_pos);
            printf("the last packet has right header,but the content is definitely incomplete\n");
            return 2;
        }
        else{
            for(int i=0;i<(header_length-1);i++){
                if(!memcmp(HEADER,tmp+i,header_length-1-i)){
                    packet_length+=i;
                    in.seekg(init_pos);
                    return 0;
                }
            }

            packet_length+=(header_length-1);
            in.seekg(init_pos);
            printf("the last packet has right header,but the content maybe incomplete(can't find next header)\n");
            return 1;
        }
    }
    else{
        printf("fatal error!an error occured in reading the file:");
        return -1;
    }

}

UInt_t convert_dateTotimecode(const char* date)//example:20140410162029
{
    if(strlen(date) != 14){
        printf("error: invalid date format!\n");
        printf("The expected format is: 20140410162029\n");
        return -1;
    }

    UInt_t timecode;
    int leapyear_num;
    int totalyear_num;
    char buffer[10];
    int year,init_year;
    int month,init_month;
    int day,init_day;
    int hour,init_hour;
    int miniute,init_miniute;
    int second,init_second;

    init_year=2013;
    init_month=1;
    init_day=1;
    init_hour=0;
    init_miniute=0;
    init_second=0;
    int monthday[12]={31,28,31,30,31,30,31,31,30,31,30,31};
    int leap_monthday[12]={31,29,31,30,31,30,31,31,30,31,30,31};
    int leap_flag=0;
    //-------
    strncpy(buffer,date,4);
    buffer[4]='\0';
    year=atoi(buffer);

    strncpy(buffer,date+4,2);
    buffer[2]='\0';
    month=atoi(buffer);

    strncpy(buffer,date+6,2);
    buffer[2]='\0';
    day=atoi(buffer);

    strncpy(buffer,date+8,2);
    buffer[2]='\0';
    hour=atoi(buffer);

    strncpy(buffer,date+10,2);
    buffer[2]='\0';
    miniute=atoi(buffer);

    strncpy(buffer,date+12,2);
    buffer[2]='\0';
    second=atoi(buffer);
    //printf("input:%s\n",date);
    //printf("output:%4.d%2.d%2.d%2.d%2.d%2.d\n",year,month,day,hour,miniute,second);
    //--------
    totalyear_num=year-init_year;
    leapyear_num=0;
    int current_year;
    for(int i=0;i<totalyear_num;i++){
        current_year=init_year+i;
        if(((current_year%4)== 0)&&((current_year%100)>0)){
            leapyear_num++;
        }
        else if((current_year%400)==0){
            leapyear_num++;
        }
    }
    if(((year%4)== 0)&&((year%100)!=0)){
        leap_flag=1;
    }
    else if((year%400)==0){
        leap_flag=1;
    }
    //----------
    int month_num=month-init_month;
    int day_num=day-1;
    timecode=leapyear_num*366*24*60*60 + (totalyear_num-leapyear_num)*365*24*60*60;
    if(leap_flag){
        for(int i=0;i<month_num;i++){
            day_num+=leap_monthday[i];
        }
    }
    else{
        for(int i=0;i<month_num;i++){
            day_num+=monthday[i];
        }
    }
    timecode+=day_num*24*60*60 + hour*60*60 + miniute*60 + second;

    return timecode;
}

bool check_leapyear(int year)
{
    if((year%4==0)&&(year%100)!=0){
        return true;
    }
    else if((year%400)==0){
        return true;
    }
    else{
        return false;
    }
}

char* convert_timecodeTodate(UInt_t timecode)
{
    static char date[50];

    int totalyear_num;
    int current_year;
    int current_month;
    int current_day;
    int current_hour;
    int current_miniute;
    int current_second;
    bool currentyear_leapflag=false;
    UInt_t totalsecond;
    int init_year=2013;
    int normal_monthday[12]={31,28,31,30,31,30,31,31,30,31,30,31};
    int leap_monthday[12]={31,29,31,30,31,30,31,31,30,31,30,31};
    int noramlyear_second=365*24*3600;
    int leapyear_second=366*24*3600;

    int temp_year;
    totalyear_num=-1;
    totalsecond=0;
    int tmp_totalsecond;
    while(timecode>=totalsecond){
        totalyear_num++;
        temp_year=init_year+totalyear_num;
        tmp_totalsecond=totalsecond;
        if(check_leapyear(temp_year)){
            totalsecond+=leapyear_second;
        }
        else{
            totalsecond+=noramlyear_second;
        }
    }
    current_year=temp_year;

    if(check_leapyear(current_year)){
        currentyear_leapflag=true;
    }
    int current_daynum=(timecode-tmp_totalsecond)/(24*3600);
    current_month=0;
    int temp_daynum=0;
    int temp_daynum_before;
    for(int i=0;i<12;i++){
        current_month++;
        if(currentyear_leapflag){
            temp_daynum+=leap_monthday[i];
        }
        else{
            temp_daynum+=normal_monthday[i];
        }
        if(temp_daynum > current_daynum){
            current_day=current_daynum-temp_daynum_before+1;
            break;
        }
        temp_daynum_before=temp_daynum;
    }

    int remaining_second=(timecode-tmp_totalsecond)%(24*3600);
    current_hour=remaining_second/3600;
    remaining_second-=current_hour*3600;
    current_miniute=remaining_second/60;
    current_second=remaining_second%60;

    sprintf(date,"%d/%d/%d/%d/%d/%d",current_year,current_month,current_day,current_hour,current_miniute,current_second);

    return date;
}

//---you need to decode the header in the main decoding function---------------------------------------
int convert_psd_scidata(FILE* fp,const int datatype,const Char_t* parentDir,const Char_t* infile,const Char_t* outDir,const Char_t* outfile="raw.root",const char *start_date=0,const char* stop_date=0)
{
    bool minimum_flag;
    int trigger_id_before=0;
    int trigger_id_current=0;
    const char HEADER[4]={0xe2,0x25,0x08,0x13};//science data application id
    const int TYPEID[4]={0x20,0x60,0xa0,0x30};//normal,compressed,calibration,smaller
    const Int_t header_length=2+6+8;//E225+PrimaryHeader+SecondaryHeader

    Char_t infname[400],outfname[400];
    sprintf(infname,"%s/%s",parentDir,infile);
    sprintf(outfname,"%s/%s",outDir,outfile);

    fprintf(fp,"Decoding PSD science data:\n");
    fprintf(fp,"Binary Source File:\n\t%s\n",infname);
    fprintf(fp,"Root Target File:\n\t%s\n",outfname);
    fprintf(fp,"Log Information:\n");

    TFile *f=new TFile(outfname,"RECREATE");
    if(f->IsZombie()){
        fprintf(fp,"\t%s can not been recreated!\n",outfname);
        delete f;
        return -1;
    }
    else if(f->IsOpen()){
        std::ifstream in;
        in.open(infname,std::ios_base::in | std::ios_base::binary);
        if(!in.is_open()){
            std::cout << "can't open"<< infname << std::endl;
            fprintf(fp,"\tcan't open %s\n",infname);
        }
//-----------------------------------------------------------------------
        TH1F *hxpos[90];
        TH1F *hxneg[90];
        TH1F *hypos[90];
        TH1F *hyneg[90];
        for(int i=0;i<90;i++){
            if(datatype == 0xa0){
                hxpos[i]=new TH1F(Form("xpos_%d",i+1),Form("xpos_%d",i+1),17000,-0.5,16999.5);
                hxneg[i]=new TH1F(Form("xneg_%d",i+1),Form("xneg_%d",i+1),17000,-0.5,16999.5);
                hypos[i]=new TH1F(Form("ypos_%d",i+1),Form("ypos_%d",i+1),17000,-0.5,16999.5);
                hyneg[i]=new TH1F(Form("yneg_%d",i+1),Form("yneg_%d",i+1),17000,-0.5,16999.5);
            }
            else{
                hxpos[i]=new TH1F(Form("xpos_%d",i+1),Form("xpos_%d",i+1),4000,-0.5,3999.5);
                hxneg[i]=new TH1F(Form("xneg_%d",i+1),Form("xneg_%d",i+1),4000,-0.5,3999.5);
                hypos[i]=new TH1F(Form("ypos_%d",i+1),Form("ypos_%d",i+1),4000,-0.5,3999.5);
                hyneg[i]=new TH1F(Form("yneg_%d",i+1),Form("yneg_%d",i+1),4000,-0.5,3999.5);
            }
        }

        TTree *tree = new TTree("PSD","PSD Testing event mode");
        Int_t Xpos[90],Ypos[90],Xneg[90],Yneg[90];
        int trigger_id[4],trigger_status[4],packet_id;
        UShort_t status_code[4];
        UShort_t error_code;
        ULong64_t packet_num;
        UInt_t time_second;
        UShort_t time_millisecond;
        tree->Branch("xpos",Xpos,"xpos[90]/I");
        tree->Branch("ypos",Ypos,"ypos[90]/I");
        tree->Branch("xneg",Xneg,"xneg[90]/I");
        tree->Branch("yneg",Yneg,"yneg[90]/I");
        tree->Branch("trigger_id",trigger_id,"trigger_id[4]/I");
        tree->Branch("trigger_status",trigger_status,"trigger_status[4]/I");
        tree->Branch("packet_id",&packet_id,"packet_id/I");
        tree->Branch("time_second",&time_second,"time_second/i");
        tree->Branch("time_millisecond",&time_millisecond,"time_millisecond/s");
        tree->Branch("error_code",&error_code,"error_code/s");
        tree->Branch("status_code",status_code,"status_code[4]/s");
        tree->Branch("packet_num",&packet_num,"packet_num/l");
//----------------------------------------------------------------------------------------
    unsigned int total_event_num=0;
    unsigned int event_num=0;
    packet_num=0;
    unsigned int init_packet_num=0;
    ULong64_t total_eventlength=0;
    bool init_packet_flag=false;
    UInt_t start_timecode,stop_timecode;
    if(start_date && stop_date){
        start_timecode=convert_dateTotimecode(start_date);
        stop_timecode=convert_dateTotimecode(stop_date);
    }
    else{
        start_timecode=-1;
        stop_timecode=-1;
    }
    unsigned int unprocessed_num=0;
    int type_id;
    char tmp[4];
    char packet_buffer[1000];
    int error_flag,packet_length,event_flag;
    //--if the file header has a incomplete packet,this packet will not be processed.
    fprintf(fp,"\tStart Converting:\n");
    if(in.read(tmp,4)){
        if(!memcmp(HEADER,tmp,4)){
            in.seekg(0,in.beg);
        }
        else{
            unprocessed_num++;
            packet_num++;
            fprintf(fp,"\tpacket_%d not processed: first packet is incomplete\n",packet_num);
            in.seekg(0,in.beg);
            error_flag=find_nextheader(in,packet_length);
            if(error_flag == -1){
                fprintf(fp,"\terror! reading first packet\n");
                exit(EXIT_FAILURE);
            }
            in.read(packet_buffer,packet_length);
        }
    }
    //--the last packet will not be processed,whether complete or not
    while(!in.eof()){
        error_flag=find_nextheader(in,packet_length);
        if(error_flag == 3){
            //no bytes available in the input stream,eof reached
            //i.e. the previous read has already get all the bytes in the file
            continue;
        }
        else if(error_flag == -1){
            fprintf(fp,"\tpacket_num=%d,event_num=%d,has been processed successfully\n",packet_num,event_num);
            exit(EXIT_FAILURE);
        }
        else{
            packet_num++;
            if(packet_num%3000==0)
            {
                printf("%d packets converted!\n",packet_num);
            }
            in.read(packet_buffer,packet_length);
            if(error_flag==0){
                packet_id=((packet_buffer[4]&0x3F)<<8) + (packet_buffer[5]&0xFF);
                time_second=((packet_buffer[10]&0xFF)<<24)+((packet_buffer[11]&0xFF)<<16)+((packet_buffer[12]&0xFF)<<8)+(packet_buffer[13]&0xFF);
                time_millisecond=((packet_buffer[14]&0xFF)<<8)+(packet_buffer[15]&0xFF);

                if((start_timecode==-1 && stop_timecode==-1) || ((start_timecode <= time_second) && (stop_timecode >= time_second))){
                    if(!init_packet_flag){
                        init_packet_num=packet_num;
                        init_packet_flag=true;
                    }
                    if(datatype==0x30){
                        minimum_flag=true;
                        event_flag=convert_psd_packet_scidata(minimum_flag,packet_buffer+header_length,packet_length-header_length,type_id,trigger_id,trigger_status,status_code,Xpos,Ypos,Xneg,Yneg);
                    }
                    else{
                        minimum_flag=false;
                        event_flag=convert_psd_packet_scidata(minimum_flag,packet_buffer+header_length,packet_length-header_length,type_id,trigger_id,trigger_status,status_code,Xpos,Ypos,Xneg,Yneg);
                    }
                    error_code=event_flag;
                    if(!event_flag){
                        total_event_num++;
                        if(type_id == datatype){
                            event_num++;
			    total_eventlength+=(packet_length-header_length);
                            for(int i=0;i<90;i++){
                                hxpos[i]->Fill(Xpos[i]);
                                hypos[i]->Fill(Ypos[i]);
                                hxneg[i]->Fill(Xneg[i]);
                                hyneg[i]->Fill(Yneg[i]);
                            }
                            tree->Fill();

                            //check all the status code
                            //and print out the warnings
                            if(check_triggerid(trigger_id)){
                                fprintf(fp,"\tpacket %d:\n",packet_num);
                                fprintf(fp,"\t\twarning: unmatched trigger_id in this packet\n");
                            }
                            check_triggerstatus(trigger_status,packet_num,true);
                            check_statuscode(status_code,packet_num);

                            trigger_id_current=trigger_id[0];
                            if(trigger_id_before==0xFFF){
                                trigger_id_before=0;
                            }
                            else{
                                trigger_id_before++;
                            }
                            if(packet_num>init_packet_num){
                                if(trigger_id_before != trigger_id_current){
                                    fprintf(fp,"\tpacket %d:\n",packet_num);
                                    fprintf(fp,"\tincontinuous trigger_id in this packet(%d,%d)\n",trigger_id_before,trigger_id_current);
                                }
                            }
                            trigger_id_before=trigger_id_current;

                            //--
                            /*
                            if(event_num%3000==0)
                            {
                                printf("%d events converted!\n",event_num);
                            }
                            */
                        }
                    }
                    else{
                        unprocessed_num++;
                        fprintf(fp,"\tpacket_%d not processed: bad event data,default value will be filled\n",packet_num);

                        for(int i=0;i<90;i++){
                            hxpos[i]->Fill(Xpos[i]);
                            hypos[i]->Fill(Ypos[i]);
                            hxneg[i]->Fill(Xneg[i]);
                            hyneg[i]->Fill(Yneg[i]);
                        }
                        tree->Fill();

                        continue;
                    }
                }
            }
            else{
                //by default,the last packet will not be processed
                unprocessed_num++;
                fprintf(fp,"\tpacket_%d not processed: the last packet\n",packet_num);
            }
        }
    }

    fprintf(fp,"\tConverting End:\n");
    fprintf(fp,"\tTotally,%d packets has been found.\n",packet_num);
    fprintf(fp,"\t%d of them are not processed\n",unprocessed_num);
    fprintf(fp,"\tTotally,%d events has been converted successfully which begins at packet_%d.\n",total_event_num,init_packet_num);
    fprintf(fp,"\t%d of them are 0x%x events\n",event_num,datatype);
    fprintf(fp,"\ttotal event length: %d\n",total_eventlength);
    fprintf(fp,"\tCompression Ratio: %.4f\n",1.0-(total_eventlength+event_num*4*4)/(194.0*4.0*event_num));
    
    in.close();
    f->Write(0,TObject::kOverwrite);//overwrite AutoSave keys
    f->Close();
    //printf("f=%x\n",f);
    //delete f;
  }
  else{
    fprintf(fp,"\terror: %s can not be created!\n",outfile);
    delete f;
    return -1;
  }

    return 0;
}

//---Method 2: end--------------------------------------------------------------------------------------------------------------------------------------------

/*********** END ************************************/

/************ Analyze Part*********************************/
//----calibration------------------------
int sorting(float *a,const int nfound)
{
    float tmp[15];
    int index[15];
    for(int i=0;i<nfound;i++){
        tmp[i]=a[i];
    }
    TMath::Sort(nfound,a,index,false);
    for(int i=0;i<nfound;i++){
        a[i]=tmp[index[i]];
    }

    return 0;
}

TF1 *calibration_linearfit(TGraph* g,float &INL,float *par,float* calib_voltage,Int_t datapoints)
{
    char FunName[100];
    sprintf(FunName,"ffit_%s",g->GetName());

    TF1* ffit=(TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
    if(ffit) delete ffit;

    g->Fit("pol1","q0","goff",50,1250);
    ffit=(TF1*)g->GetFunction("pol1");
    ffit->SetName(FunName);

    //float calib_voltage[15]={100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500};
    double* adc_value;
    adc_value=g->GetY();
    par[0]=ffit->GetParameter(0);
    par[1]=ffit->GetParameter(1);
    float delta[15];
    for(int i=0;i<datapoints;i++){
        //delta[i]=TMath::Abs(ffit->Eval(calib_voltage[i]) - adc_value[i]);
        delta[i]=TMath::Abs(par[0]+par[1]*calib_voltage[i]-adc_value[i]);
    }

    INL=TMath::MaxElement(datapoints,delta)/15000;

    return ffit;
}

int fit_calibration(const Char_t* parentDir,const Char_t* infile,const Char_t* outDir,const Char_t* outfile="calib_result",const char* config_file="calib.config")
{
    gStyle->SetOptStat(0);
    const int id_PMT[82]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,
                      23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};

    TH1F *hxpos[90];
    TH1F *hxneg[90];
    TH1F *hypos[90];
    TH1F *hyneg[90];
     TFile* file_in=new TFile(Form("%s/%s",parentDir,infile));
     for(int i=0;i<90;i++){
         hxpos[i]=(TH1F*)file_in->Get(Form("xpos_%d",i+1));
         hxpos[i]=(TH1F*)file_in->Get(Form("xpos_%d",i+1));

         hxneg[i]=(TH1F*)file_in->Get(Form("xneg_%d",i+1));
         hxneg[i]=(TH1F*)file_in->Get(Form("xneg_%d",i+1));

         hypos[i]=(TH1F*)file_in->Get(Form("ypos_%d",i+1));
         hypos[i]=(TH1F*)file_in->Get(Form("ypos_%d",i+1));

         hyneg[i]=(TH1F*)file_in->Get(Form("yneg_%d",i+1));
         hyneg[i]=(TH1F*)file_in->Get(Form("yneg_%d",i+1));
     }

     TFile *file_out=new TFile(Form("%s/%s.root",outDir,outfile),"recreate");
     TTree* tree_out=new TTree("calib_fit","Fitting result of each PMT channel");
     Float_t xpos_mean[15],xpos_sigma[15],*xpos_tmp;
     Float_t xneg_mean[15],xneg_sigma[15],*xneg_tmp;
     Float_t ypos_mean[15],ypos_sigma[15],*ypos_tmp;
     Float_t yneg_mean[15],yneg_sigma[15],*yneg_tmp;
     tree_out->Branch("xpos_mean",xpos_mean,"xpos_mean[15]/F");
     tree_out->Branch("xpos_sigma",xpos_sigma,"xpos_sigma[15]/F");
     tree_out->Branch("xneg_mean",xneg_mean,"xneg_mean[15]/F");
     tree_out->Branch("xneg_sigma",xneg_sigma,"xneg_sigma[15]/F");
     tree_out->Branch("ypos_mean",ypos_mean,"ypos_mean[15]/F");
     tree_out->Branch("ypos_sigma",ypos_sigma,"ypos_sigma[15]/F");
     tree_out->Branch("yneg_mean",yneg_mean,"yneg_mean[15]/F");
     tree_out->Branch("yneg_sigma",yneg_sigma,"yneg_sigma[15]/F");


     //Float_t range_tmp[15]={1600,2700,3900,5000,6100,7200,8300,9400,10500,11700,12900,14000,14960,15600,15830};
     TSpectrum *s=new TSpectrum(18);
     TF1* ffit=0;
     TH1F *hxpos_new,*hxneg_new,*hypos_new,*hyneg_new;
     int xpos_nfound[90],xneg_nfound[90],ypos_nfound[90],yneg_nfound[90];
     Float_t sigma=15.0;

     TCanvas* can_fit=new TCanvas("can_fit","can_fit",1000,1000);
     can_fit->Divide(2,2);
     can_fit->Print(Form("%s/%s_fit.pdf[",outDir,outfile));
     for(int i=0;i<90;i++){
         can_fit->cd(1);
         hxpos_new=(TH1F*)hxpos[i]->Rebin(10,"hxpos_new");
         xpos_nfound[i]=s->Search(hxpos_new);
         /*if(xpos_nfound[i] !=15 ){
             //printf("error! not 15 peaks(hxpos[%d],nfound=%d)\n",i,xpos_nfound[i]);
             //exit(EXIT_FAILURE);
         }*/
         xpos_tmp=s->GetPositionX();
         for(int j=0;j<xpos_nfound[i];j++){
             xpos_mean[j]=xpos_tmp[j];
         }
         sorting(xpos_mean,xpos_nfound[i]);
         for(int j=0;j<xpos_nfound[i];j++){
             //xpos_mean[j]=range_tmp[j];
             hxpos[i]->Fit("gaus","q0","",xpos_mean[j]-3*sigma,xpos_mean[j]+3*sigma);
             ffit=(TF1*)hxpos[i]->GetFunction("gaus");
             xpos_mean[j]=ffit->GetParameter(1);
             xpos_sigma[j]=ffit->GetParameter(2);
             //delete ffit;
         }
         //delete hxpos_new;

         //s=new TSpectrum(15);
         can_fit->cd(2);
         hxneg_new=(TH1F*)hxneg[i]->Rebin(10,"hxneg_new");
         xneg_nfound[i]=s->Search(hxneg_new);
         /*if(xneg_nfound[i] !=15 ){
             //printf("error! not 15 peaks(hxneg[%d],nfound=%d)\n",i,xneg_nfound[i]);
             //exit(EXIT_FAILURE);
         }*/
         xneg_tmp=s->GetPositionX();
         for(int j=0;j<xneg_nfound[i];j++){
             xneg_mean[j]=xneg_tmp[j];
         }
         sorting(xneg_mean,xneg_nfound[i]);
         for(int j=0;j<xneg_nfound[i];j++){
             //xneg_mean[j]=range_tmp[j];
             hxneg[i]->Fit("gaus","q0","",xneg_mean[j]-3*sigma,xneg_mean[j]+3*sigma);
             ffit=(TF1*)hxneg[i]->GetFunction("gaus");
             xneg_mean[j]=ffit->GetParameter(1);
             xneg_sigma[j]=ffit->GetParameter(2);
             //delete ffit;
         }
         //delete hxneg_new;

        // s=new TSpectrum(15);
         can_fit->cd(3);
         hypos_new=(TH1F*)hypos[i]->Rebin(10,"hypos_new");
         ypos_nfound[i]=s->Search(hypos_new);
         /*if(ypos_nfound[i] !=15 ){
             //printf("error! not 15 peaks(hypos[%d],nfound=%d)\n",i,ypos_nfound[i]);
             //exit(EXIT_FAILURE);
         }*/
         ypos_tmp=s->GetPositionX();
         for(int j=0;j<ypos_nfound[i];j++){
             ypos_mean[j]=ypos_tmp[j];
         }
         sorting(ypos_mean,ypos_nfound[i]);
         for(int j=0;j<ypos_nfound[i];j++){
             //ypos_mean[j]=range_tmp[j];
             hypos[i]->Fit("gaus","q0","",ypos_mean[j]-3*sigma,ypos_mean[j]+3*sigma);
             ffit=(TF1*)hypos[i]->GetFunction("gaus");
             ypos_mean[j]=ffit->GetParameter(1);
             ypos_sigma[j]=ffit->GetParameter(2);
             //delete ffit;
         }
         //delete hypos_new;

         can_fit->cd(4);
         hyneg_new=(TH1F*)hyneg[i]->Rebin(10,"hyneg_new");
         yneg_nfound[i]=s->Search(hyneg_new);
         /*if(yneg_nfound[i] !=15 ){
             //printf("error! not 15 peaks(hyneg[%d],nfound=%d)\n",i,yneg_nfound[i]);
             //exit(EXIT_FAILURE);
         }*/
         yneg_tmp=s->GetPositionX();
         for(int j=0;j<yneg_nfound[i];j++){
             yneg_mean[j]=yneg_tmp[j];
         }
         sorting(yneg_mean,yneg_nfound[i]);
         for(int j=0;j<yneg_nfound[i];j++){
             //yneg_mean[j]=range_tmp[j];
             hyneg[i]->Fit("gaus","q0","",yneg_mean[j]-3*sigma,yneg_mean[j]+3*sigma);
             ffit=(TF1*)hyneg[i]->GetFunction("gaus");
             yneg_mean[j]=ffit->GetParameter(1);
             yneg_sigma[j]=ffit->GetParameter(2);
             //delete ffit;
         }
         //delete hyneg_new;

         can_fit->Print(Form("%s/%s_fit.pdf",outDir,outfile));
         tree_out->Fill();
         delete hxpos_new;
         delete hxneg_new;
         delete hypos_new;
         delete hyneg_new;
     }
     can_fit->Print(Form("%s/%s_fit.pdf]",outDir,outfile));
     //--------------------------------------------------------------------------------------------------------------
     char buffer[100];
     sprintf(buffer,"%s/%s.csv",outDir,outfile);
     FILE* fp=fopen(buffer,"w");
     fprintf(fp,"xpos\tINL\tp0\tp1\t\txneg\t\t\t\typos\t\t\t\tyneg\t\t\t\t\n");
     tree_out->SetBranchAddress("xpos_mean",xpos_mean);
     tree_out->SetBranchAddress("xneg_mean",xneg_mean);
     tree_out->SetBranchAddress("ypos_mean",ypos_mean);
     tree_out->SetBranchAddress("yneg_mean",yneg_mean);

     TGraph* gxpos[90];
     TGraph* gxneg[90];
     TGraph* gypos[90];
     TGraph* gyneg[90];
     //Float_t calib_voltage[10]={100,200,400,500,700,800,1000,1100,1300,1500};
     Float_t calib_voltage[20];
     Int_t calib_datapoints;
     Float_t xpos_INL[90],xneg_INL[90],ypos_INL[90],yneg_INL[90];
     Float_t par[3];

     //--read calib.config
     FILE* fp_calib=fopen(config_file,"r");
     if(!fp_calib){
         printf("can't open calibration config file\n");
         return -1;
     }
     fscanf(fp_calib,"datapoints: %d\n",&calib_datapoints);
     printf("working dir: %s\n",gSystem->WorkingDirectory());
     printf("datapoints: %d\n",calib_datapoints);
     for(int i=0;i<calib_datapoints;i++){
         fscanf(fp_calib,"%f\n",calib_voltage+i);
         printf("%f\n",calib_voltage[i]);
     }
     fclose(fp_calib);
     //--
     Int_t max_datapoints;
     for(int i=0;i<calib_datapoints;i++){
         if(calib_voltage[calib_datapoints-1-i]<=1100){
             max_datapoints=calib_datapoints-i;
             break;
         }
     }
     TCanvas* can=new TCanvas("can","can",1000,1000);
     can->Divide(2,2);
     can->Print(Form("%s/%s.pdf[",outDir,outfile));
     for(int i=0;i<90;i++){
         tree_out->GetEntry(i);
         gxpos[i]=new TGraph(xpos_nfound[i],calib_voltage,xpos_mean);
         gxpos[i]->SetNameTitle(Form("gxpos_%d",i+1),Form("gxpos_%d",i+1));
         ffit=calibration_linearfit(gxpos[i],xpos_INL[i],par,calib_voltage,max_datapoints);
         fprintf(fp,"%d,%.5f,%.2f,%.2f,,",i+1,xpos_INL[i],par[0],par[1]);
         can->cd(1);
         gxpos[i]->Draw("A*");
         ffit->Draw("lsame");

         gxneg[i]=new TGraph(xneg_nfound[i],calib_voltage,xneg_mean);
         gxneg[i]->SetNameTitle(Form("gxneg_%d",i+1),Form("gxneg_%d",i+1));
         ffit=calibration_linearfit(gxneg[i],xneg_INL[i],par,calib_voltage,max_datapoints);
         fprintf(fp,"%.5f,%.2f,%.2f,,",xneg_INL[i],par[0],par[1]);
         can->cd(2);
         gxneg[i]->Draw("A*");
         ffit->Draw("lsame");

         gypos[i]=new TGraph(ypos_nfound[i],calib_voltage,ypos_mean);
         gypos[i]->SetNameTitle(Form("gypos_%d",i+1),Form("gypos_%d",i+1));
         ffit=calibration_linearfit(gypos[i],ypos_INL[i],par,calib_voltage,max_datapoints);
         fprintf(fp,"%.5f,%.2f,%.2f,,",ypos_INL[i],par[0],par[1]);
         can->cd(3);
         gypos[i]->Draw("A*");
         ffit->Draw("lsame");

         gyneg[i]=new TGraph(yneg_nfound[i],calib_voltage,yneg_mean);
         gyneg[i]->SetNameTitle(Form("gyneg_%d",i+1),Form("gyneg_%d",i+1));
         ffit=calibration_linearfit(gyneg[i],yneg_INL[i],par,calib_voltage,max_datapoints);
         fprintf(fp,"%.5f,%.2f,%.2f\n",yneg_INL[i],par[0],par[1]);
         can->cd(4);
         gyneg[i]->Draw("A*");
         ffit->Draw("lsame");

         can->Print(Form("%s/%s.pdf",outDir,outfile));
     }
     can->Print(Form("%s/%s.pdf]",outDir,outfile));
     fclose(fp);

     file_out->cd();
     tree_out->Write();
     for(int i=0;i<90;i++){
         gxpos[i]->Write();
         gxneg[i]->Write();
         gypos[i]->Write();
         gyneg[i]->Write();
     }
     delete file_in;
     delete file_out;

    delete can_fit;
     delete can;

    return 0;
}
//----end---------------------------------

//----pedestal----------------------------
int getpedseed_event(const char* parentDir,const char* infile,const char* outDir,const char* outfile)
{
    gStyle->SetOptFit(11);

    float xpos_mean[90]={351.177,311.881,281.001,336.517,258.712,313.785,367.483,306.686,365.061,400.932,389.930,414.645,381.898,411.290,457.176,472.376,463.821,539.654,495.049,496.161,441.246,514.411,469.810,480.274,360.767,386.358,352.114,449.317,322.460,344.191,420.976,343.341,240.283,314.475,384.396,319.270,298.015,267.687,222.540,274.143,163.888,253.895,150.799,200.748,175.790,145.786,524.328,515.608,526.098,480.832,485.753,494.254,464.564,435.246,343.468,312.436,341.823,331.957,229.788,180.386,167.866,143.693,157.214,148.327,113.189,226.599,165.258,195.581,410.142,399.905,425.209,444.228,467.265,398.202,483.074,502.183,482.189,552.946,492.655,465.126,451.151,434.841,447.516,462.986,427.054,431.159,441.765,404.810,370.627,377.385};
    float xneg_mean[90]={446.400,520.242,514.340,474.583,496.813,511.273,451.954,431.075,434.599,440.726,489.202,435.936,499.171,509.421,467.255,535.061,496.474,563.510,574.983,541.864,537.884,496.657,596.976,366.517,391.959,375.442,351.044,405.371,442.040,453.049,475.211,478.180,435.166,421.951,420.496,458.119,457.642,411.283,510.404,547.936,484.632,467.074,450.268,434.619,498.571,393.499,420.707,422.288,433.151,393.430,429.871,415.825,436.522,474.949,460.671,484.791,485.606,453.007,471.585,500.918,497.301,502.769,396.582,485.579,465.391,514.925,499.670,426.066,424.271,463.562,454.416,464.023,409.735,455.579,354.684,375.137,407.772,453.645,461.221,462.112,465.888,447.161,534.270,467.812,467.011,472.202,423.727,475.979,441.521,414.418};
    float ypos_mean[90]={452.297,396.719,465.625,479.740,528.990,451.267,480.754,451.002,424.205,458.843,438.643,442.601,436.487,414.472,413.737,369.849,492.507,466.578,337.399,441.559,374.001,292.447,310.534,376.861,384.398,464.878,458.101,409.130,432.157,468.180,525.836,426.869,496.900,449.711,540.653,521.539,493.202,489.634,462.673,486.457,399.784,497.665,358.482,479.102,343.498,421.276,399.754,497.368,455.586,454.538,413.200,431.225,420.818,345.452,323.384,384.369,422.946,351.598,326.737,311.787,384.763,351.652,284.502,269.942,281.036,241.743,251.200,222.739,462.892,475.289,442.673,449.252,420.579,475.795,479.726,412.265,460.324,410.219,446.704,432.104,474.956,529.460,500.751,453.041,417.644,420.573,377.888,386.825,331.463,392.985};
    float yneg_mean[90]={417.906,472.811,455.049,452.441,453.630,426.665,413.650,458.842,491.942,415.372,438.802,455.469,401.188,413.670,430.686,444.817,489.856,400.641,473.606,387.638,458.840,453.257,439.884,315.324,256.090,265.857,255.013,298.694,286.077,283.444,285.370,303.172,315.953,275.825,234.171,293.651,249.059,276.830,270.094,341.402,282.482,244.047,257.179,300.565,285.312,233.734,427.574,432.091,328.272,419.002,395.342,355.594,297.112,425.499,384.624,395.578,408.127,322.187,396.572,409.241,430.655,393.992,374.721,431.261,373.104,342.076,386.307,329.733,499.101,476.227,497.179,514.796,493.029,511.467,495.702,494.944,478.117,528.368,518.425,441.704,552.574,457.136,484.718,451.821,499.494,413.811,472.878,330.334,307.009,269.479};
    float* tmp_mean;

    TString label[4]={"xpos","xneg","ypos","yneg"};
    Int_t idneg[8]={2,21,25,44,48,67,70,89};

    Char_t infname[200],outfname[200],txtfname[200];
    sprintf(infname,"%s/%s",parentDir,infile);
    TFile* file_in=new TFile(infname);
    sprintf(outfname,"%s/%s.root",outDir,outfile);
    TFile* file_out=new TFile(outfname,"RECREATE");

    Float_t mean,sigma,sum,mean_sigma;
    Int_t channel;
    Float_t xpeaks;
    Float_t xmin,xmax;
    Char_t channel_label[10];

    TH1F* hist;
    TF1 *fgaus;
    FILE* fp=0;
    TTree* tree_in;
    TTree* tree_out;

    tree_in=(TTree*)file_in->Get("PSD");
    for(Int_t fee_index=0;fee_index<4;fee_index++){

        sprintf(txtfname,"%s/%s_%s.csv",outDir,outfile,label[fee_index].Data());
        //printf("open file:%s\n",txtfname);
        fp=fopen(txtfname,"w");
        fprintf(fp,"Label/C:Channel/I:Mean/F:Sigma/F\n");

        tree_out=new TTree(Form("%s_ped",label[fee_index].Data()),Form("%s channel pedestals",label[fee_index].Data()));
        tree_out->Branch("label",channel_label,"label/C");
        tree_out->Branch("channel",&channel,"channel/I");
        tree_out->Branch("mean",&mean,"mean/F");
        tree_out->Branch("sigma",&sigma,"sigma/F");

        TCanvas* canvas=new TCanvas(Form("%s_canvas",label[fee_index].Data()),Form("%s: Pedestal fitting of All channels",label[fee_index].Data()),900,1000);
        canvas->Print(Form("%s/%s_%s.pdf[",outDir,outfile,label[fee_index].Data()));

        sum=0;
        switch(fee_index)
        {
        case 0:
            tmp_mean=xpos_mean;
            break;
        case 1:
            tmp_mean=xneg_mean;
            break;
        case 2:
            tmp_mean=ypos_mean;
            break;
        case 3:
            tmp_mean=yneg_mean;
            break;
        default:
            return -1;
        }
        for(Int_t ch_id=0;ch_id<90;ch_id++){
            channel=ch_id+1;
            sprintf(channel_label,"%s_%d",label[fee_index].Data(),channel);
            hist=new TH1F(Form("%s_%d",label[fee_index].Data(),channel),Form("%s_%d",label[fee_index].Data(),channel),300,0,1000);
            tree_in->Project(Form("%s_%d",label[fee_index].Data(),channel),Form("%s[%d]",label[fee_index].Data(),ch_id));

            xpeaks=tmp_mean[ch_id];
            xmin=xpeaks-30;
            xmax=xpeaks+30;
            //xpeaks=hist->GetBinContent(hist->GetMaximumBin());
            //xmin=xpeaks-100;
            //xmax=xpeaks+100;
            fgaus=new TF1("fgaus","gaus",xmin,xmax);
            fgaus->SetNpx(1000);
            hist->Fit("fgaus","rq");

            canvas->Print(Form("%s/%s_%s.pdf",outDir,outfile,label[fee_index].Data()));
            mean=fgaus->GetParameter(1);
            sigma=fgaus->GetParameter(2);
            tree_out->Fill();
            fprintf(fp,"%s,\t%d,\t%f,\t%f\n",channel_label,channel,mean,sigma);
            delete hist;

            sum+=sigma;
            for(int i=0;i<8;i++){
                if(ch_id == idneg[i]){
                    sum-=sigma;
                    break;
                }
            }
        }

        mean_sigma=sum/82;
        fprintf(fp,"mean_sig,\t,\t,%f\n",mean_sigma);
        //printf("%s\n",outDir);
        canvas->Print(Form("%s/%s_%s.pdf]",outDir,outfile,label[fee_index].Data()));

        fclose(fp);
        file_out->cd();
        tree_out->Write();
        delete canvas;

    }
    delete file_in;
    delete file_out;

    return 0;
}

int getped_u(TTree* tree_in,const UInt_t startEntry,const UInt_t eNum,Float_t *xpos_mean_fit,Float_t *xpos_sigma_fit,Float_t* xneg_mean_fit,Float_t* xneg_sigma_fit,
                         Float_t *ypos_mean_fit,Float_t *ypos_sigma_fit,Float_t *yneg_mean_fit,Float_t *yneg_sigma_fit)
{
    Float_t xpos_mean[90]={351.177,311.881,281.001,336.517,258.712,313.785,367.483,306.686,365.061,400.932,389.930,414.645,381.898,411.290,457.176,472.376,463.821,539.654,495.049,496.161,441.246,514.411,469.810,480.274,360.767,386.358,352.114,449.317,322.460,344.191,420.976,343.341,240.283,314.475,384.396,319.270,298.015,267.687,222.540,274.143,163.888,253.895,150.799,200.748,175.790,145.786,524.328,515.608,526.098,480.832,485.753,494.254,464.564,435.246,343.468,312.436,341.823,331.957,229.788,180.386,167.866,143.693,157.214,148.327,113.189,226.599,165.258,195.581,410.142,399.905,425.209,444.228,467.265,398.202,483.074,502.183,482.189,552.946,492.655,465.126,451.151,434.841,447.516,462.986,427.054,431.159,441.765,404.810,370.627,377.385};
    Float_t xneg_mean[90]={446.400,520.242,514.340,474.583,496.813,511.273,451.954,431.075,434.599,440.726,489.202,435.936,499.171,509.421,467.255,535.061,496.474,563.510,574.983,541.864,537.884,496.657,596.976,366.517,391.959,375.442,351.044,405.371,442.040,453.049,475.211,478.180,435.166,421.951,420.496,458.119,457.642,411.283,510.404,547.936,484.632,467.074,450.268,434.619,498.571,393.499,420.707,422.288,433.151,393.430,429.871,415.825,436.522,474.949,460.671,484.791,485.606,453.007,471.585,500.918,497.301,502.769,396.582,485.579,465.391,514.925,499.670,426.066,424.271,463.562,454.416,464.023,409.735,455.579,354.684,375.137,407.772,453.645,461.221,462.112,465.888,447.161,534.270,467.812,467.011,472.202,423.727,475.979,441.521,414.418};
    Float_t ypos_mean[90]={452.297,396.719,465.625,479.740,528.990,451.267,480.754,451.002,424.205,458.843,438.643,442.601,436.487,414.472,413.737,369.849,492.507,466.578,337.399,441.559,374.001,292.447,310.534,376.861,384.398,464.878,458.101,409.130,432.157,468.180,525.836,426.869,496.900,449.711,540.653,521.539,493.202,489.634,462.673,486.457,399.784,497.665,358.482,479.102,343.498,421.276,399.754,497.368,455.586,454.538,413.200,431.225,420.818,345.452,323.384,384.369,422.946,351.598,326.737,311.787,384.763,351.652,284.502,269.942,281.036,241.743,251.200,222.739,462.892,475.289,442.673,449.252,420.579,475.795,479.726,412.265,460.324,410.219,446.704,432.104,474.956,529.460,500.751,453.041,417.644,420.573,377.888,386.825,331.463,392.985};
    Float_t yneg_mean[90]={417.906,472.811,455.049,452.441,453.630,426.665,413.650,458.842,491.942,415.372,438.802,455.469,401.188,413.670,430.686,444.817,489.856,400.641,473.606,387.638,458.840,453.257,439.884,315.324,256.090,265.857,255.013,298.694,286.077,283.444,285.370,303.172,315.953,275.825,234.171,293.651,249.059,276.830,270.094,341.402,282.482,244.047,257.179,300.565,285.312,233.734,427.574,432.091,328.272,419.002,395.342,355.594,297.112,425.499,384.624,395.578,408.127,322.187,396.572,409.241,430.655,393.992,374.721,431.261,373.104,342.076,386.307,329.733,499.101,476.227,497.179,514.796,493.029,511.467,495.702,494.944,478.117,528.368,518.425,441.704,552.574,457.136,484.718,451.821,499.494,413.811,472.878,330.334,307.009,269.479};

    Float_t *tmp_mean;
    Float_t *tmp_mean_fit,*tmp_sigma_fit;
    TString label[4]={"xpos","xneg","ypos","yneg"};

    Int_t channel;
    Float_t xmin,xmax;
    Char_t channel_label[10];

    TH1F* hist;
    TF1 *fgaus;
    for(Int_t fee_index=0;fee_index<4;fee_index++){
        switch(fee_index)
        {
        case 0:
            tmp_mean=xpos_mean;
            tmp_mean_fit=xpos_mean_fit;
            tmp_sigma_fit=xpos_sigma_fit;
            break;
        case 1:
            tmp_mean=xneg_mean;
            tmp_mean_fit=xneg_mean_fit;
            tmp_sigma_fit=xneg_sigma_fit;
            break;
        case 2:
            tmp_mean=ypos_mean;
            tmp_mean_fit=ypos_mean_fit;
            tmp_sigma_fit=ypos_sigma_fit;
            break;
        case 3:
            tmp_mean=yneg_mean;
            tmp_mean_fit=yneg_mean_fit;
            tmp_sigma_fit=yneg_sigma_fit;
            break;
        default:
            return -1;
        }

        for(Int_t ch_id=0;ch_id<90;ch_id++){
            channel=ch_id+1;
            sprintf(channel_label,"%s_%d",label[fee_index].Data(),channel);
            hist=new TH1F(channel_label,channel_label,2000,-0.5,1999.5);
            tree_in->Project(channel_label,Form("%s[%d]",label[fee_index].Data(),ch_id),"","",eNum,startEntry);

            xmin=tmp_mean[ch_id]-200;
            xmax=tmp_mean[ch_id]+200;
            fgaus=new TF1("fgaus","gaus",xmin,xmax);
            //fgaus->SetNpx(1000);
            hist->Fit("fgaus","rq0","goff");
            tmp_mean_fit[ch_id]=fgaus->GetParameter(1);
            tmp_sigma_fit[ch_id]=fgaus->GetParameter(2);

            delete hist;
            delete fgaus;
        }
    }

    return 0;
}

Float_t getmean_u(const Float_t *data,const Int_t *valid_channel,const Int_t channel_num,Bool_t exclude_flag=false)
{
    Float_t sum=0.0;
    Float_t mean=0.0;
    if(!exclude_flag){
        for(int i=0;i<channel_num;i++){
            sum+=data[valid_channel[i]];
        }
        mean=sum/channel_num;
    }
    else{
        for(int i=0;i<90;i++){
            sum+=data[i];
        }
        for(int i=0;i<channel_num;i++){
            sum-=data[valid_channel[i]];
        }
        mean=sum/(90-channel_num);
    }

    return mean;
}

int draw_pedVStime(const char* parentDir,const char * infile,const char* outDir,const char* outfile,const Int_t eNum=2000)
{
    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};
    int idempty[8]={2,21,25,44,48,67,70,89};
    TString label[4]={"xpos","xneg","ypos","yneg"};

    Float_t mean_fit[4][90];
    Float_t sigma_fit[4][90];
    Float_t tmp_mean_mean,tmp_sigma_mean;

    char infname[200],outfname[200];
    sprintf(infname,"%s/%s",parentDir,infile);
    TFile* file_in=new TFile(infname);
    TString prefix;
    prefix=outfile;
    prefix.ReplaceAll(".root","");
    sprintf(outfname,"%s/%s.root",outDir,prefix.Data());
    TFile* file_out=new TFile(outfname,"update");
    TGraph* gr_dy8_mean[4],*gr_dy8_sigma[4];
    TGraph* gr_dy5_mean[4],*gr_dy5_sigma[4];
    TGraph* gr_validch_mean[4],*gr_validch_sigma[4];
    TGraph* gr_emptych_mean[4],*gr_emptych_sigma[4];
    for(int i=0;i<4;i++){
        gr_dy5_mean[i]=new TGraph();
        gr_dy5_sigma[i]=new TGraph();
        gr_dy5_mean[i]->SetName(Form("%s_dy5_pedMeanVStime",label[i].Data()));
        gr_dy5_sigma[i]->SetName(Form("%s_dy5_pedSigmaVStime",label[i].Data()));
        gr_dy5_mean[i]->SetTitle(Form("%s_dy5: pedestal_Mean VS time",label[i].Data()));
        gr_dy5_sigma[i]->SetTitle(Form("%s_dy5: pedestal_Sigma VS time",label[i].Data()));
        gr_dy5_mean[i]->SetMaximum(650.);gr_dy5_mean[i]->SetMinimum(250.);
        gr_dy5_sigma[i]->SetMaximum(60.);gr_dy5_sigma[i]->SetMinimum(0);
        gr_dy8_mean[i]=new TGraph();
        gr_dy8_sigma[i]=new TGraph();
        gr_dy8_mean[i]->SetName(Form("%s_dy8_pedMeanVStime",label[i].Data()));
        gr_dy8_sigma[i]->SetName(Form("%s_dy8_pedSigmaVStime",label[i].Data()));
        gr_dy8_mean[i]->SetTitle(Form("%s_dy8: pedestal_Mean VS time",label[i].Data()));
        gr_dy8_sigma[i]->SetTitle(Form("%s_dy8: pedestal_Sigma VS time",label[i].Data()));
        gr_dy8_mean[i]->SetMaximum(650.);gr_dy8_mean[i]->SetMinimum(250.);
        gr_dy8_sigma[i]->SetMaximum(60.);gr_dy8_sigma[i]->SetMinimum(0);
        gr_validch_mean[i]=new TGraph();
        gr_validch_sigma[i]=new TGraph();
        gr_validch_mean[i]->SetName(Form("%s_validch_pedMeanVStime",label[i].Data()));
        gr_validch_sigma[i]->SetName(Form("%s_validch_pedSigmaVStime",label[i].Data()));
        gr_validch_mean[i]->SetTitle(Form("%s_validch: pedestal_Mean VS time",label[i].Data()));
        gr_validch_sigma[i]->SetTitle(Form("%s_validch: pedestal_Sigma VS time",label[i].Data()));
        gr_validch_mean[i]->SetMaximum(650.);gr_validch_mean[i]->SetMinimum(250.);
        gr_validch_sigma[i]->SetMaximum(60.);gr_validch_sigma[i]->SetMinimum(0);
        gr_emptych_mean[i]=new TGraph();
        gr_emptych_sigma[i]=new TGraph();
        gr_emptych_mean[i]->SetName(Form("%s_emptych_pedMeanVStime",label[i].Data()));
        gr_emptych_sigma[i]->SetName(Form("%s_emptych_pedSigmaVStime",label[i].Data()));
        gr_emptych_mean[i]->SetTitle(Form("%s_emptych: pedestal_Mean VS time",label[i].Data()));
        gr_emptych_sigma[i]->SetTitle(Form("%s_emptych: pedestal_Sigma VS time",label[i].Data()));
        gr_emptych_mean[i]->SetMaximum(650.);gr_emptych_mean[i]->SetMinimum(250.);
        gr_emptych_sigma[i]->SetMaximum(20.);gr_emptych_sigma[i]->SetMinimum(0);
    }

    TTree* tree_in=(TTree*)file_in->Get("PSD");
    UInt_t entries=tree_in->GetEntries();
    Int_t nPoints;
    Int_t entries_last;
    if((entries%eNum) > 900){
        nPoints=entries/eNum + 1;
        entries_last=entries%eNum;
    }
    else{
        nPoints=entries/eNum;
        entries_last=entries%eNum + eNum;
    }
    for(int i=0;i<nPoints;i++){
        if(i == (nPoints-1)){
            getped_u(tree_in,i*eNum,entries_last,mean_fit[0],sigma_fit[0],mean_fit[1],sigma_fit[1],mean_fit[2],sigma_fit[2],mean_fit[3],sigma_fit[3]);
        }
        else{
            getped_u(tree_in,i*eNum,eNum,mean_fit[0],sigma_fit[0],mean_fit[1],sigma_fit[1],mean_fit[2],sigma_fit[2],mean_fit[3],sigma_fit[3]);
        }

        for(int fee_id=0;fee_id<4;fee_id++){
            tmp_mean_mean=getmean_u(mean_fit[fee_id],id5,41);
            tmp_sigma_mean=getmean_u(sigma_fit[fee_id],id5,41);
            gr_dy5_mean[fee_id]->SetPoint(i,i+1,tmp_mean_mean);
            gr_dy5_sigma[fee_id]->SetPoint(i,i+1,tmp_sigma_mean);

            tmp_mean_mean=getmean_u(mean_fit[fee_id],id8,41);
            tmp_sigma_mean=getmean_u(sigma_fit[fee_id],id8,41);
            gr_dy8_mean[fee_id]->SetPoint(i,i+1,tmp_mean_mean);
            gr_dy8_sigma[fee_id]->SetPoint(i,i+1,tmp_sigma_mean);

            tmp_mean_mean=getmean_u(mean_fit[fee_id],idempty,8,true);
            tmp_sigma_mean=getmean_u(sigma_fit[fee_id],idempty,8,true);
            gr_validch_mean[fee_id]->SetPoint(i,i+1,tmp_mean_mean);
            gr_validch_sigma[fee_id]->SetPoint(i,i+1,tmp_sigma_mean);

            tmp_mean_mean=getmean_u(mean_fit[fee_id],idempty,8);
            tmp_sigma_mean=getmean_u(sigma_fit[fee_id],idempty,8);
            gr_emptych_mean[fee_id]->SetPoint(i,i+1,tmp_mean_mean);
            gr_emptych_sigma[fee_id]->SetPoint(i,i+1,tmp_sigma_mean);
        }
    }

    TCanvas* can=new TCanvas("can","can",1000,500);
    can->Divide(2,1);
    can->Print(Form("%s/%s_pedVStime.pdf[",outDir,prefix.Data()));
    file_out->cd();
    for(int i=0;i<4;i++){
        gr_dy5_mean[i]->Write(0,TObject::kOverwrite);
        gr_dy5_sigma[i]->Write(0,TObject::kOverwrite);
        gr_dy8_mean[i]->Write(0,TObject::kOverwrite);
        gr_dy8_sigma[i]->Write(0,TObject::kOverwrite);
        gr_validch_mean[i]->Write(0,TObject::kOverwrite);
        gr_validch_sigma[i]->Write(0,TObject::kOverwrite);
        gr_emptych_mean[i]->Write(0,TObject::kOverwrite);
        gr_emptych_sigma[i]->Write(0,TObject::kOverwrite);

        can->cd(1);
        gr_validch_mean[i]->Draw("A*");
        can->cd(2);
        gr_validch_sigma[i]->Draw("A*");
        can->Print(Form("%s/%s_pedVStime.pdf",outDir,prefix.Data()));
    }
    can->Print(Form("%s/%s_pedVStime.pdf]",outDir,prefix.Data()));

    delete can;
    delete file_in;
    delete file_out;

    return 0;
}

int draw_rel(const Char_t* testfile,const Char_t* reffile,const Char_t* outdir,const char* outName)
{
    TString label[4]={"xpos","xneg","ypos","yneg"};
    Float_t mean,sigma,mean_ref,sigma_ref,ratio_mean,ratio_sigma,delta_mean,delta_sigma;
    Float_t mean_max,sigma_max,mean_min,sigma_min,tmp_mean_max,tmp_mean_min,tmp_sigma_max,tmp_sigma_min;
    Int_t channel;
    //Char_t channel_label[10];

    TFile *file_in=new TFile(Form("%s",testfile));
    TFile *file_ref=new TFile(Form("%s",reffile));

    TTree* tree_test;
    TTree* tree_ref;
    TGraph* graph_mean[4];
    TGraph* graph_sigma[4];
    TGraph* graph_ref_mean[4];
    TGraph* graph_ref_sigma[4];
    TGraph* graph_delta_mean[4];
    TGraph* graph_delta_sigma[4];

    TCanvas* can=new TCanvas("can","Pedestals",900,1200);
    can->Divide(2,3);
    can->Print(Form("%s/pedrel.pdf[",outdir));

    unsigned int nentries;
    for(int fee_id=0;fee_id<4;fee_id++){
        tmp_mean_max=1.02;tmp_sigma_max=1.05;
        tmp_mean_min=0.98;tmp_sigma_min=0.95;

        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id].Data()));
        //tree_test->SetBranchAddress("label",channel_label);
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("mean",&mean);
        tree_test->SetBranchAddress("sigma",&sigma);

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id].Data()));
        tree_ref->SetBranchAddress("mean",&mean_ref);
        tree_ref->SetBranchAddress("sigma",&sigma_ref);
        tree_ref->BuildIndex("channel");

        nentries=tree_test->GetEntries();
        mean_max=tree_test->GetMaximum("mean");
        mean_min=tree_test->GetMinimum("mean");
        sigma_max=tree_test->GetMaximum("sigma");
        sigma_min=tree_test->GetMinimum("sigma");
//--------------------------------------------------------------------------------------------
        graph_mean[fee_id]=new TGraph(nentries);
        graph_mean[fee_id]->SetName(Form("g_%s_mean",label[fee_id].Data()));
        graph_mean[fee_id]->SetTitle(Form("%s: mean value of pedestal",label[fee_id].Data()));
        graph_mean[fee_id]->SetMaximum(mean_max+50);
        graph_mean[fee_id]->SetMinimum(0);
        graph_ref_mean[fee_id]=new TGraph(nentries);
        graph_ref_mean[fee_id]->SetName(Form("g_%s_ref_mean",label[fee_id].Data()));
        graph_ref_mean[fee_id]->SetTitle(Form("%s: ratio of mean value of ped",label[fee_id].Data()));
        graph_ref_mean[fee_id]->SetMarkerStyle(26);
        graph_delta_mean[fee_id]=new TGraph(nentries);
        graph_delta_mean[fee_id]->SetName(Form("g_%s_delta_mean",label[fee_id].Data()));
        graph_delta_mean[fee_id]->SetTitle(Form("%s: delta of mean value of pedestal",label[fee_id].Data()));
        graph_delta_mean[fee_id]->SetMarkerStyle(25);
        //graph_delta_mean->SetMaximum(mean_max+50);
        //graph_delta_mean->SetMinimum(0);

        graph_sigma[fee_id]=new TGraph(nentries);
        graph_sigma[fee_id]->SetName(Form("g_%s_sigma",label[fee_id].Data()));
        graph_sigma[fee_id]->SetTitle(Form("%s: sigma of pedestal",label[fee_id].Data()));
        graph_sigma[fee_id]->SetMaximum(sigma_max+1);
        graph_sigma[fee_id]->SetMinimum(0);
        graph_ref_sigma[fee_id]=new TGraph(nentries);
        graph_ref_sigma[fee_id]->SetName(Form("g_%s_ref_sigma",label[fee_id].Data()));
        graph_ref_sigma[fee_id]->SetTitle(Form("%s: ratio of sigma value of ped",label[fee_id].Data()));
        graph_ref_sigma[fee_id]->SetMarkerStyle(26);
        graph_delta_sigma[fee_id]=new TGraph(nentries);
        graph_delta_sigma[fee_id]->SetName(Form("g_%s_delta_sigma",label[fee_id].Data()));
        graph_delta_sigma[fee_id]->SetTitle(Form("%s: delta of sigma of pedestal",label[fee_id].Data()));
        graph_delta_sigma[fee_id]->SetMarkerStyle(25);

        for(int ch_id=0;ch_id<nentries;ch_id++){
            tree_test->GetEntry(ch_id);
            tree_ref->GetEntryWithIndex(channel);

            graph_mean[fee_id]->SetPoint(ch_id,ch_id+1,mean);
            graph_sigma[fee_id]->SetPoint(ch_id,ch_id+1,sigma);

            ratio_mean=mean/mean_ref;
            ratio_sigma=sigma/sigma_ref;
            graph_ref_mean[fee_id]->SetPoint(ch_id,ch_id+1,ratio_mean);
            graph_ref_sigma[fee_id]->SetPoint(ch_id,ch_id+1,ratio_sigma);
            if(ratio_mean > tmp_mean_max)   tmp_mean_max=ratio_mean;
            if(ratio_mean < tmp_mean_min)   tmp_mean_min=ratio_mean;
            if(ratio_sigma > tmp_sigma_max) tmp_sigma_max=ratio_sigma;
            if(ratio_sigma < tmp_sigma_min) tmp_sigma_min=ratio_sigma;

            delta_mean=mean-mean_ref;
            delta_sigma=sigma-sigma_ref;
            graph_delta_mean[fee_id]->SetPoint(ch_id,ch_id+1,delta_mean);
            graph_delta_sigma[fee_id]->SetPoint(ch_id,ch_id+1,delta_sigma);
        }

        can->cd(1);
        graph_mean[fee_id]->Draw("A*");
        can->cd(2);
        graph_sigma[fee_id]->Draw("A*");
        can->cd(3);
        graph_ref_mean[fee_id]->SetMaximum(tmp_mean_max+0.01);
        graph_ref_mean[fee_id]->SetMinimum(tmp_mean_min-0.01);
        graph_ref_mean[fee_id]->Draw("AP");
        can->cd(4);
        graph_ref_sigma[fee_id]->SetMaximum(tmp_sigma_max+0.05);
        graph_ref_sigma[fee_id]->SetMinimum(tmp_sigma_min-0.05);
        graph_ref_sigma[fee_id]->Draw("AP");
        can->cd(5);
        graph_delta_mean[fee_id]->Draw("AP");
        can->cd(6);
        graph_delta_sigma[fee_id]->Draw("AP");

        can->Print(Form("%s/pedrel.pdf",outdir));
    }

    can->Print(Form("%s/pedrel.pdf]",outdir));
    TFile *file_out=new TFile(Form("%s/%s",outdir,outName),"update");
    for(int i=0;i<4;i++){
        graph_mean[i]->Write(0,TObject::kOverwrite);
        graph_sigma[i]->Write(0,TObject::kOverwrite);
        graph_ref_mean[i]->Write(0,TObject::kOverwrite);
        graph_ref_sigma[i]->Write(0,TObject::kOverwrite);
        graph_delta_mean[i]->Write(0,TObject::kOverwrite);
        graph_delta_sigma[i]->Write(0,TObject::kOverwrite);
    }
    delete file_out;
    delete file_in;
    delete file_ref;
    delete can;

    return 0;
}

int draw_relp(const Char_t* testfile,const Char_t* reffile,const Char_t* outdir,const char* outName)
{
    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};

    TString label[4]={"xpos","xneg","ypos","yneg"};
    Float_t sigma,sigma_ref,ratio_sigma;
    Float_t mean,mean_ref,ratio_mean;
    Int_t channel;

    TFile *file_in=new TFile(Form("%s",testfile));
    TFile *file_ref=new TFile(Form("%s",reffile));
    TFile *file_out=new TFile(Form("%s/%s",outdir,outName),"update");

    TTree* tree_test;
    TTree* tree_ref;

    TMultiGraph *mg_sigma_x;
    TMultiGraph *mg_ref_sigma_x;
    TMultiGraph *mg_sigma_y;
    TMultiGraph *mg_ref_sigma_y;
    TGraph* graphp_sigma_xpos;
    TGraph* graphp_ref_sigma_xpos;
    TGraph* graphp_sigma_ypos;
    TGraph* graphp_ref_sigma_ypos;

    TGraph* graphp_sigma_xneg;
    TGraph* graphp_ref_sigma_xneg;
    TGraph* graphp_sigma_yneg;
    TGraph* graphp_ref_sigma_yneg;


    //-----------------------------------------------------------------------------------------------
    graphp_sigma_xpos=new TGraph(41);
    graphp_sigma_xpos->SetTitle("X_layer_dy8: sigma of pedestal");
    graphp_ref_sigma_xpos=new TGraph(41);
    graphp_ref_sigma_xpos->SetTitle("X_layer_dy8: ratio of sigma value of ped");
    graphp_ref_sigma_xpos->SetMarkerStyle(26);
    graphp_sigma_ypos=new TGraph(41);
    graphp_sigma_ypos->SetTitle("Y_layer_dy8: sigma of pedestal");
    graphp_ref_sigma_ypos=new TGraph(41);
    graphp_ref_sigma_ypos->SetTitle("Y_layer_dy8: ratio of sigma value of ped");
    graphp_ref_sigma_ypos->SetMarkerStyle(26);

    graphp_sigma_xneg=new TGraph(41);
    graphp_sigma_xneg->SetTitle("X_layer_dy8: sigma of pedestal");
    graphp_ref_sigma_xneg=new TGraph(41);
    graphp_ref_sigma_xneg->SetTitle("X_layer_dy8: ratio of sigma value of ped");
    graphp_ref_sigma_xneg->SetMarkerStyle(26);
    graphp_sigma_yneg=new TGraph(41);
    graphp_sigma_yneg->SetTitle("Y_layer_dy8: sigma of pedestal");
    graphp_ref_sigma_yneg=new TGraph(41);
    graphp_ref_sigma_yneg->SetTitle("Y_layer_dy8: ratio of sigma value of ped");
    graphp_ref_sigma_yneg->SetMarkerStyle(26);

    mg_sigma_x = new TMultiGraph();
    mg_sigma_x->SetName("mg_x_sigma_dy8");
    mg_sigma_x->SetTitle("X_layer_dy8: sigma of pedestal");
    mg_sigma_x->Add(graphp_sigma_xpos);
    mg_sigma_x->Add(graphp_sigma_xneg);
    mg_sigma_y = new TMultiGraph();
    mg_sigma_y->SetName("mg_y_sigma_dy8");
    mg_sigma_y->SetTitle("Y_layer_dy8: sigma of pedestal");
    mg_sigma_y->Add(graphp_sigma_ypos);
    mg_sigma_y->Add(graphp_sigma_yneg);
    mg_ref_sigma_x = new TMultiGraph();
    mg_ref_sigma_x->SetName("mg_x_ref_sigma_dy8");
    mg_ref_sigma_x->SetTitle("X_layer_dy8: ratio of sigma value of ped");
    mg_ref_sigma_x->Add(graphp_ref_sigma_xpos);
    mg_ref_sigma_x->Add(graphp_ref_sigma_xneg);
    mg_ref_sigma_y = new TMultiGraph();
    mg_ref_sigma_y->SetName("mg_y_ref_sigma_dy8");
    mg_ref_sigma_y->SetTitle("Y_layer_dy8: ratio of sigma value of ped");
    mg_ref_sigma_y->Add(graphp_ref_sigma_ypos);
    mg_ref_sigma_y->Add(graphp_ref_sigma_yneg);

    for(int fee_id=0;fee_id<2;fee_id++){
        //-----------------X--------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("sigma",&sigma);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id].Data()));
        tree_ref->SetBranchAddress("sigma",&sigma_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id8[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_xpos->SetPoint(ch_id,ch_id+1,sigma);
            }
            else{
                graphp_sigma_xneg->SetPoint(ch_id,41-ch_id,-sigma);
            }
            ratio_sigma=sigma/sigma_ref;
            if(fee_id==0){
                graphp_ref_sigma_xpos->SetPoint(ch_id,ch_id+1,ratio_sigma);
            }
            else{
                graphp_ref_sigma_xneg->SetPoint(ch_id,41-ch_id,-ratio_sigma);
            }
        }
        //-----------------Y---------------------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("sigma",&sigma);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_ref->SetBranchAddress("sigma",&sigma_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id8[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_ypos->SetPoint(ch_id,41-ch_id,sigma);
            }
            else{
                graphp_sigma_yneg->SetPoint(ch_id,ch_id+1,-sigma);
            }
            ratio_sigma=sigma/sigma_ref;
            if(fee_id==0){
                graphp_ref_sigma_ypos->SetPoint(ch_id,41-ch_id,ratio_sigma);
            }
            else{
                graphp_ref_sigma_yneg->SetPoint(ch_id,ch_id+1,-ratio_sigma);
            }
        }
    }

    file_out->cd();
    TCanvas* canp=new TCanvas("canp","Rel_ped VS positon",900,900);
    canp->Divide(2,2);
    canp->cd(1);
    mg_sigma_x->Draw("A*");
    mg_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(2);
    mg_ref_sigma_x->Draw("AP");
    mg_ref_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(3);
    mg_sigma_y->Draw("A*");
    mg_sigma_y->Write(0,TObject::kOverwrite);
    canp->cd(4);
    mg_ref_sigma_y->Draw("AP");
    mg_ref_sigma_y->Write(0,TObject::kOverwrite);
    canp->Print(Form("%s/pedrel_sigma_dy8.pdf",outdir));

    delete mg_sigma_x;
    delete mg_sigma_y;
    delete mg_ref_sigma_x;
    delete mg_ref_sigma_y;
    delete canp;
    //------------------------------------------------
    graphp_sigma_xpos=new TGraph(41);
    graphp_sigma_xpos->SetTitle("X_layer_dy5: sigma of pedestal");
    graphp_ref_sigma_xpos=new TGraph(41);
    graphp_ref_sigma_xpos->SetTitle("X_layer_dy5: ratio of sigma value of ped");
    graphp_ref_sigma_xpos->SetMarkerStyle(26);
    graphp_sigma_ypos=new TGraph(41);
    graphp_sigma_ypos->SetTitle("Y_layer_dy5: sigma of pedestal");
    graphp_ref_sigma_ypos=new TGraph(41);
    graphp_ref_sigma_ypos->SetTitle("Y_layer_dy5: ratio of sigma value of ped");
    graphp_ref_sigma_ypos->SetMarkerStyle(26);

    graphp_sigma_xneg=new TGraph(41);
    graphp_sigma_xneg->SetTitle("X_layer_dy5: sigma of pedestal");
    graphp_ref_sigma_xneg=new TGraph(41);
    graphp_ref_sigma_xneg->SetTitle("X_layer_dy5: ratio of sigma value of ped");
    graphp_ref_sigma_xneg->SetMarkerStyle(26);
    graphp_sigma_yneg=new TGraph(41);
    graphp_sigma_yneg->SetTitle("Y_layer_dy5: sigma of pedestal");
    graphp_ref_sigma_yneg=new TGraph(41);
    graphp_ref_sigma_yneg->SetTitle("Y_layer_dy5: ratio of sigma value of ped");
    graphp_ref_sigma_yneg->SetMarkerStyle(26);

    mg_sigma_x = new TMultiGraph();
    mg_sigma_x->SetName("mg_x_sigma_dy5");
    mg_sigma_x->SetTitle("X_layer_dy5: sigma of pedestal");
    mg_sigma_x->Add(graphp_sigma_xpos);
    mg_sigma_x->Add(graphp_sigma_xneg);
    mg_sigma_y = new TMultiGraph();
    mg_sigma_y->SetName("mg_y_sigma_dy5");
    mg_sigma_y->SetTitle("Y_layer_dy5: sigma of pedestal");
    mg_sigma_y->Add(graphp_sigma_ypos);
    mg_sigma_y->Add(graphp_sigma_yneg);
    mg_ref_sigma_x = new TMultiGraph();
    mg_ref_sigma_x->SetName("mg_x_ref_sigma_dy5");
    mg_ref_sigma_x->SetTitle("X_layer_dy5: ratio of sigma value of ped");
    mg_ref_sigma_x->Add(graphp_ref_sigma_xpos);
    mg_ref_sigma_x->Add(graphp_ref_sigma_xneg);
    mg_ref_sigma_y = new TMultiGraph();
    mg_ref_sigma_y->SetName("mg_y_ref_sigma_dy5");
    mg_ref_sigma_y->SetTitle("Y_layer_dy5: ratio of sigma value of ped");
    mg_ref_sigma_y->Add(graphp_ref_sigma_ypos);
    mg_ref_sigma_y->Add(graphp_ref_sigma_yneg);

    for(int fee_id=0;fee_id<2;fee_id++){
        //-----------------X--------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("sigma",&sigma);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id].Data()));
        tree_ref->SetBranchAddress("sigma",&sigma_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id5[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_xpos->SetPoint(ch_id,ch_id+1,sigma);
            }
            else{
                graphp_sigma_xneg->SetPoint(ch_id,41-ch_id,-sigma);
            }
            ratio_sigma=sigma/sigma_ref;
            if(fee_id==0){
                graphp_ref_sigma_xpos->SetPoint(ch_id,ch_id+1,ratio_sigma);
            }
            else{
                graphp_ref_sigma_xneg->SetPoint(ch_id,41-ch_id,-ratio_sigma);
            }
        }
        //-----------------Y---------------------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("sigma",&sigma);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_ref->SetBranchAddress("sigma",&sigma_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id5[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_ypos->SetPoint(ch_id,41-ch_id,sigma);
            }
            else{
                graphp_sigma_yneg->SetPoint(ch_id,ch_id+1,-sigma);
            }
            ratio_sigma=sigma/sigma_ref;
            if(fee_id==0){
                graphp_ref_sigma_ypos->SetPoint(ch_id,41-ch_id,ratio_sigma);
            }
            else{
                graphp_ref_sigma_yneg->SetPoint(ch_id,ch_id+1,-ratio_sigma);
            }
        }
    }

    file_out->cd();
    canp=new TCanvas("canp","Rel_ped VS positon",900,900);
    canp->Divide(2,2);
    canp->cd(1);
    mg_sigma_x->Draw("A*");
    mg_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(2);
    mg_ref_sigma_x->Draw("AP");
    mg_ref_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(3);
    mg_sigma_y->Draw("A*");
    mg_sigma_y->Write(0,TObject::kOverwrite);
    canp->cd(4);
    mg_ref_sigma_y->Draw("AP");
    mg_ref_sigma_y->Write(0,TObject::kOverwrite);
    canp->Print(Form("%s/pedrel_sigma_dy5.pdf",outdir));

    delete mg_sigma_x;
    delete mg_sigma_y;
    delete mg_ref_sigma_x;
    delete mg_ref_sigma_y;
    delete canp;
    //----------------------------------------------------------------------------------------------
    graphp_sigma_xpos=new TGraph(41);
    graphp_sigma_xpos->SetTitle("X_layer_dy8: mean of pedestal");
    graphp_ref_sigma_xpos=new TGraph(41);
    graphp_ref_sigma_xpos->SetTitle("X_layer_dy8: ratio of mean value of ped");
    graphp_ref_sigma_xpos->SetMarkerStyle(26);
    graphp_sigma_ypos=new TGraph(41);
    graphp_sigma_ypos->SetTitle("Y_layer_dy8: mean of pedestal");
    graphp_ref_sigma_ypos=new TGraph(41);
    graphp_ref_sigma_ypos->SetTitle("Y_layer_dy8: ratio of mean value of ped");
    graphp_ref_sigma_ypos->SetMarkerStyle(26);

    graphp_sigma_xneg=new TGraph(41);
    graphp_sigma_xneg->SetTitle("X_layer_dy8: mean of pedestal");
    graphp_ref_sigma_xneg=new TGraph(41);
    graphp_ref_sigma_xneg->SetTitle("X_layer_dy8: ratio of mean value of ped");
    graphp_ref_sigma_xneg->SetMarkerStyle(26);
    graphp_sigma_yneg=new TGraph(41);
    graphp_sigma_yneg->SetTitle("Y_layer_dy8: mean of pedestal");
    graphp_ref_sigma_yneg=new TGraph(41);
    graphp_ref_sigma_yneg->SetTitle("Y_layer_dy8: ratio of mean value of ped");
    graphp_ref_sigma_yneg->SetMarkerStyle(26);

    mg_sigma_x = new TMultiGraph();
    mg_sigma_x->SetName("mg_x_mean_dy8");
    mg_sigma_x->SetTitle("X_layer_dy8: mean of pedestal");
    mg_sigma_x->Add(graphp_sigma_xpos);
    mg_sigma_x->Add(graphp_sigma_xneg);
    mg_sigma_y = new TMultiGraph();
    mg_sigma_y->SetName("mg_y_mean_dy8");
    mg_sigma_y->SetTitle("Y_layer_dy8: mean of pedestal");
    mg_sigma_y->Add(graphp_sigma_ypos);
    mg_sigma_y->Add(graphp_sigma_yneg);
    mg_ref_sigma_x = new TMultiGraph();
    mg_ref_sigma_x->SetName("mg_x_ref_mean_dy8");
    mg_ref_sigma_x->SetTitle("X_layer_dy8: ratio of mean value of ped");
    mg_ref_sigma_x->Add(graphp_ref_sigma_xpos);
    mg_ref_sigma_x->Add(graphp_ref_sigma_xneg);
    mg_ref_sigma_y = new TMultiGraph();
    mg_ref_sigma_y->SetName("mg_y_ref_mean_dy8");
    mg_ref_sigma_y->SetTitle("Y_layer_dy8: ratio of mean value of ped");
    mg_ref_sigma_y->Add(graphp_ref_sigma_ypos);
    mg_ref_sigma_y->Add(graphp_ref_sigma_yneg);

    for(int fee_id=0;fee_id<2;fee_id++){
        //-----------------X--------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("mean",&mean);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id].Data()));
        tree_ref->SetBranchAddress("mean",&mean_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id8[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_xpos->SetPoint(ch_id,ch_id+1,mean);
            }
            else{
                graphp_sigma_xneg->SetPoint(ch_id,41-ch_id,-mean);
            }
            ratio_mean=mean/mean_ref;
            if(fee_id==0){
                graphp_ref_sigma_xpos->SetPoint(ch_id,ch_id+1,ratio_mean);
            }
            else{
                graphp_ref_sigma_xneg->SetPoint(ch_id,41-ch_id,-ratio_mean);
            }
        }
        //-----------------Y---------------------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("mean",&mean);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_ref->SetBranchAddress("mean",&mean_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id8[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_ypos->SetPoint(ch_id,41-ch_id,mean);
            }
            else{
                graphp_sigma_yneg->SetPoint(ch_id,ch_id+1,-mean);
            }
            ratio_mean=mean/mean_ref;
            if(fee_id==0){
                graphp_ref_sigma_ypos->SetPoint(ch_id,41-ch_id,ratio_mean);
            }
            else{
                graphp_ref_sigma_yneg->SetPoint(ch_id,ch_id+1,-ratio_mean);
            }
        }
    }

    file_out->cd();
    canp=new TCanvas("canp","Rel_ped VS positon",900,900);
    canp->Divide(2,2);
    canp->cd(1);
    mg_sigma_x->Draw("A*");
    mg_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(2);
    mg_ref_sigma_x->Draw("AP");
    mg_ref_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(3);
    mg_sigma_y->Draw("A*");
    mg_sigma_y->Write(0,TObject::kOverwrite);
    canp->cd(4);
    mg_ref_sigma_y->Draw("AP");
    mg_ref_sigma_y->Write(0,TObject::kOverwrite);
    canp->Print(Form("%s/pedrel_mean_dy8.pdf",outdir));

    delete mg_sigma_x;
    delete mg_sigma_y;
    delete mg_ref_sigma_x;
    delete mg_ref_sigma_y;
    delete canp;
    //----------------------------------------------------------------
    graphp_sigma_xpos=new TGraph(41);
    graphp_sigma_xpos->SetTitle("X_layer_dy5: mean of pedestal");
    graphp_ref_sigma_xpos=new TGraph(41);
    graphp_ref_sigma_xpos->SetTitle("X_layer_dy5: ratio of mean value of ped");
    graphp_ref_sigma_xpos->SetMarkerStyle(26);
    graphp_sigma_ypos=new TGraph(41);
    graphp_sigma_ypos->SetTitle("Y_layer_dy5: mean of pedestal");
    graphp_ref_sigma_ypos=new TGraph(41);
    graphp_ref_sigma_ypos->SetTitle("Y_layer_dy5: ratio of mean value of ped");
    graphp_ref_sigma_ypos->SetMarkerStyle(26);

    graphp_sigma_xneg=new TGraph(41);
    graphp_sigma_xneg->SetTitle("X_layer_dy5: mean of pedestal");
    graphp_ref_sigma_xneg=new TGraph(41);
    graphp_ref_sigma_xneg->SetTitle("X_layer_dy5: ratio of mean value of ped");
    graphp_ref_sigma_xneg->SetMarkerStyle(26);
    graphp_sigma_yneg=new TGraph(41);
    graphp_sigma_yneg->SetTitle("Y_layer_dy5: mean of pedestal");
    graphp_ref_sigma_yneg=new TGraph(41);
    graphp_ref_sigma_yneg->SetTitle("Y_layer_dy5: ratio of mean value of ped");
    graphp_ref_sigma_yneg->SetMarkerStyle(26);

    mg_sigma_x = new TMultiGraph();
    mg_sigma_x->SetName("mg_x_mean_dy5");
    mg_sigma_x->SetTitle("X_layer_dy5: mean of pedestal");
    mg_sigma_x->Add(graphp_sigma_xpos);
    mg_sigma_x->Add(graphp_sigma_xneg);
    mg_sigma_y = new TMultiGraph();
    mg_sigma_y->SetName("mg_y_mean_dy5");
    mg_sigma_y->SetTitle("Y_layer_dy5: mean of pedestal");
    mg_sigma_y->Add(graphp_sigma_ypos);
    mg_sigma_y->Add(graphp_sigma_yneg);
    mg_ref_sigma_x = new TMultiGraph();
    mg_ref_sigma_x->SetName("mg_x_ref_mean_dy5");
    mg_ref_sigma_x->SetTitle("X_layer_dy5: ratio of mean value of ped");
    mg_ref_sigma_x->Add(graphp_ref_sigma_xpos);
    mg_ref_sigma_x->Add(graphp_ref_sigma_xneg);
    mg_ref_sigma_y = new TMultiGraph();
    mg_ref_sigma_y->SetName("mg_y_ref_mean_dy5");
    mg_ref_sigma_y->SetTitle("Y_layer_dy5: ratio of mean value of ped");
    mg_ref_sigma_y->Add(graphp_ref_sigma_ypos);
    mg_ref_sigma_y->Add(graphp_ref_sigma_yneg);

    for(int fee_id=0;fee_id<2;fee_id++){
        //-----------------X--------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("mean",&mean);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id].Data()));
        tree_ref->SetBranchAddress("mean",&mean_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id5[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_xpos->SetPoint(ch_id,ch_id+1,mean);
            }
            else{
                graphp_sigma_xneg->SetPoint(ch_id,41-ch_id,-mean);
            }
            ratio_mean=mean/mean_ref;
            if(fee_id==0){
                graphp_ref_sigma_xpos->SetPoint(ch_id,ch_id+1,ratio_mean);
            }
            else{
                graphp_ref_sigma_xneg->SetPoint(ch_id,41-ch_id,-ratio_mean);
            }
        }
        //-----------------Y---------------------------------------------------------------
        tree_test=(TTree*)file_in->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_test->SetBranchAddress("channel",&channel);
        tree_test->SetBranchAddress("mean",&mean);
        tree_test->BuildIndex("channel");

        tree_ref=(TTree*)file_ref->Get(Form("%s_ped",label[fee_id+2].Data()));
        tree_ref->SetBranchAddress("mean",&mean_ref);
        tree_ref->BuildIndex("channel");

        for(int ch_id=0;ch_id<41;ch_id++){
            tree_test->GetEntryWithIndex(id5[ch_id]+1);
            tree_ref->GetEntryWithIndex(channel);
            if(fee_id==0){
                graphp_sigma_ypos->SetPoint(ch_id,41-ch_id,mean);
            }
            else{
                graphp_sigma_yneg->SetPoint(ch_id,ch_id+1,-mean);
            }
            ratio_mean=mean/mean_ref;
            if(fee_id==0){
                graphp_ref_sigma_ypos->SetPoint(ch_id,41-ch_id,ratio_mean);
            }
            else{
                graphp_ref_sigma_yneg->SetPoint(ch_id,ch_id+1,-ratio_mean);
            }
        }
    }

    file_out->cd();
    canp=new TCanvas("canp","Rel_ped VS positon",900,900);
    canp->Divide(2,2);
    canp->cd(1);
    mg_sigma_x->Draw("A*");
    mg_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(2);
    mg_ref_sigma_x->Draw("AP");
    mg_ref_sigma_x->Write(0,TObject::kOverwrite);
    canp->cd(3);
    mg_sigma_y->Draw("A*");
    mg_sigma_y->Write(0,TObject::kOverwrite);
    canp->cd(4);
    mg_ref_sigma_y->Draw("AP");
    mg_ref_sigma_y->Write(0,TObject::kOverwrite);
    canp->Print(Form("%s/pedrel_mean_dy5.pdf",outdir));

    delete mg_sigma_x;
    delete mg_sigma_y;
    delete mg_ref_sigma_x;
    delete mg_ref_sigma_y;
    delete canp;
    //////////////////////////////////////////

    delete file_in;
    delete file_ref;
    delete file_out;

    return 0;
}

int draw_ped(const char* peddir,const char* pedfile,const char* outDir,const char* outName)
{
    gStyle->SetOptStat(0);

    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};
    TString label[4]={"xpos","xneg","ypos","yneg"};

    Int_t xped[90],yped[90];
    Int_t x,y;
    Int_t strip_id;

    TFile* fped=new TFile(Form("%s/%s",peddir,pedfile));
    TTree* tped=(TTree*)fped->Get("PSD");
    TFile* fout=new TFile(Form("%s/%s",outDir,outName),"update");

    TH2F *hx,*hy;
    hx=new TH2F("hx_ped_dy8","Layer_X_Dy8",41,0.5,41.5,700,0,1400);
    hy=new TH2F("hy_ped_dy8","Layer_Y_Dy8",41,0.5,41.5,700,0,1400);

    Int_t nentries=tped->GetEntries();
    for(Int_t fee_id=0;fee_id<2;fee_id++){
        tped->ResetBranchAddresses();
        tped->SetBranchAddress(label[fee_id].Data(),xped);
        tped->SetBranchAddress(label[fee_id+2].Data(),yped);
        for(int event_id=0;event_id<nentries;event_id++){
            tped->GetEntry(event_id);

            for(int ch_id=0;ch_id<41;ch_id++){
                    if(fee_id==0){
                        x=xped[id8[ch_id]];
                        strip_id=ch_id+1;
                        hx->Fill(strip_id,x);
                    }
                    else{
                        x= 1400-xped[id8[ch_id]];
                        strip_id=41-ch_id;
                        hx->Fill(strip_id,x);
                    }

                    if(fee_id==0){
                        y=yped[id8[ch_id]];
                        strip_id=41-ch_id;
                        hy->Fill(strip_id,y);
                    }
                    else{
                        y=1400-yped[id8[ch_id]];
                        strip_id=ch_id+1;
                        hy->Fill(strip_id,y);
                    }

            }
        }


    }

    TCanvas *can=new TCanvas("can","can",800,500);
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogz();
    hx->Draw("colz");
    can->cd(2);
    gPad->SetLogz();
    hy->Draw("colz");
    can->Print(Form("%s/pedVSstrip_dy8.pdf",outDir));
    fout->cd();
    hx->Write(0,TObject::kOverwrite);
    hy->Write(0,TObject::kOverwrite);
    //----------------------------------------------------
    hx->Reset();
    hy->Reset();
    hx->SetName("hx_ped_dy5");
    hx->SetTitle("Layer_X_Dy5");
    hy->SetName("hy_ped_dy5");
    hy->SetTitle("Layer_Y_Dy5");
    for(Int_t fee_id=0;fee_id<2;fee_id++){
        tped->ResetBranchAddresses();
        tped->SetBranchAddress(label[fee_id].Data(),xped);
        tped->SetBranchAddress(label[fee_id+2].Data(),yped);
        for(int event_id=0;event_id<nentries;event_id++){
            tped->GetEntry(event_id);

            for(int ch_id=0;ch_id<41;ch_id++){
                    if(fee_id==0){
                        x=xped[id5[ch_id]];
                        strip_id=ch_id+1;
                        hx->Fill(strip_id,x);
                    }
                    else{
                        x= 1400-xped[id5[ch_id]];
                        strip_id=41-ch_id;
                        hx->Fill(strip_id,x);
                    }

                    if(fee_id==0){
                        y=yped[id5[ch_id]];
                        strip_id=41-ch_id;
                        hy->Fill(strip_id,y);
                    }
                    else{
                        y=1400-yped[id5[ch_id]];
                        strip_id=ch_id+1;
                        hy->Fill(strip_id,y);
                    }

            }
        }

    }

    can->Clear();
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogz();
    hx->Draw("colz");
    can->cd(2);
    gPad->SetLogz();
    hy->Draw("colz");
    can->Print(Form("%s/pedVSstrip_dy5.pdf",outDir));
    hx->Write(0,TObject::kOverwrite);
    hy->Write(0,TObject::kOverwrite);

    delete fout;
    delete can;
    delete fped;

    return 0;
}
//-----end--------------------------------

//----MIPs--------------------------------
int draw_mip(const char* mipfile,const char* pedfile,const char* outDir,const char* outName)
{
    gStyle->SetOptStat(0);

    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};
    TString label[4]={"xpos","xneg","ypos","yneg"};

    Float_t xped_mean[41],xped_sigma[41],yped_mean[41],yped_sigma[41];
    Float_t xmean_buffer,xsigma_buffer,ymean_buffer,ysigma_buffer;
    Float_t xlimit[41],ylimit[41];
    Int_t xmip[90],ymip[90];
    Int_t xmip_tmp,ymip_tmp;
    Float_t x,y;
    Int_t strip_id;

    TFile* fped=new TFile(pedfile);
    TTree* txped=0;
    TTree* typed=0;

    TFile* fmip=new TFile(mipfile);
    TTree* tmip=(TTree*)fmip->Get("PSD");


    TH2F *hx,*hy;
    hx=new TH2F("hx_mip","",41,0.5,41.5,1600,0,8000);
    hy=new TH2F("hy_mip","",41,0.5,41.5,1600,0,8000);

    Int_t nentries=tmip->GetEntries();
    for(Int_t fee_id=0;fee_id<2;fee_id++){
        txped=(TTree*)fped->Get(Form("%s_ped",label[fee_id].Data()));
        txped->BuildIndex("channel");
        txped->SetBranchAddress("mean",&xmean_buffer);
        txped->SetBranchAddress("sigma",&xsigma_buffer);

        typed=(TTree*)fped->Get(Form("%s_ped",label[fee_id+2].Data()));
        typed->BuildIndex("channel");
        typed->SetBranchAddress("mean",&ymean_buffer);
        typed->SetBranchAddress("sigma",&ysigma_buffer);
        for(Int_t i=0;i<41;i++){
            if(fee_id==0){
                txped->GetEntryWithIndex(id8[i]+1);
                typed->GetEntryWithIndex(id8[40-i]+1);
            }
            else{
                txped->GetEntryWithIndex(id8[40-i]+1);
                typed->GetEntryWithIndex(id8[i]+1);
            }

            xped_mean[i]=xmean_buffer;
            xped_sigma[i]=xsigma_buffer;
            xlimit[i]=2*xped_sigma[i];
            yped_mean[i]=ymean_buffer;
            yped_sigma[i]=ysigma_buffer;
            ylimit[i]=2*yped_sigma[i];
        }

        tmip->ResetBranchAddresses();
        tmip->SetBranchAddress(label[fee_id].Data(),xmip);
        tmip->SetBranchAddress(label[fee_id+2].Data(),ymip);
        for(int event_id=0;event_id<nentries;event_id++){
            tmip->GetEntry(event_id);

            for(int ch_id=0;ch_id<41;ch_id++){
                if(fee_id==0){
                    xmip_tmp=xmip[id8[ch_id]];
                    ymip_tmp=ymip[id8[40-ch_id]];
                }
                else{
                    xmip_tmp=xmip[id8[40-ch_id]];
                    ymip_tmp=ymip[id8[ch_id]];
                }
                if(xmip_tmp> (xped_mean[ch_id]+xlimit[ch_id])){
                    if(fee_id==0){
                        x=xmip_tmp-xped_mean[ch_id];
                        strip_id=ch_id+1;
                        hx->Fill(strip_id,x);
                    }
                    else{
                        x= 8000-(xmip_tmp-xped_mean[ch_id]);
                        strip_id=ch_id+1;
                        hx->Fill(strip_id,x);
                    }
                }

                if(ymip_tmp > yped_mean[ch_id]+ylimit[ch_id]){
                    if(fee_id==0){
                        y=ymip_tmp-yped_mean[ch_id];
                        strip_id=ch_id+1;
                        hy->Fill(strip_id,y);
                    }
                    else{
                        y=8000-(ymip_tmp-yped_mean[ch_id]);
                        strip_id=ch_id+1;
                        hy->Fill(strip_id,y);
                    }
                }
            }
        }


    }
    TCanvas *can=new TCanvas("can_mip","can_mip",800,500);
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogz();
    hx->Draw("colz");
    gPad->Update();
    Double_t xmin=gPad->GetUxmin();
    Double_t xmax=gPad->GetUxmax();
    Double_t ymin=gPad->GetUymin();
    Double_t ymax=gPad->GetUymax();
    TGaxis *ax1=new TGaxis(xmin,ymax,xmax,ymax,0.5,41.5,510,"-L");
    ax1->SetLabelSize(hx->GetXaxis()->GetLabelSize());
    ax1->SetLabelFont(hx->GetXaxis()->GetLabelFont());
    ax1->Draw();
    hx->GetXaxis()->SetTitle("+X");
    hx->GetXaxis()->CenterTitle();
    //hx->GetXaxis()->SetTitleColor(kRed);
    //hx->GetXaxis()->SetLineColor(kRed);

    ax1->SetTitle("-X");
    ax1->CenterTitle();
    ax1->SetTitleFont(hx->GetXaxis()->GetTitleFont());
    TGaxis *ax2=new TGaxis(xmax,ymax,xmax-0.01,ymin,0,8000,510,"-");
    ax2->SetLabelSize(hx->GetYaxis()->GetLabelSize());
    ax2->SetLabelFont(hx->GetYaxis()->GetLabelFont());
    ax2->Draw();
    TPaletteAxis *paletx=(TPaletteAxis*)hx->GetListOfFunctions()->FindObject("palette");
    TGaxis* paletx_axis=paletx->GetAxis();
    paletx_axis->SetLabelSize(0);
    paletx_axis->SetTickSize(0);

    can->cd(2);
    gPad->SetLogz();
    hy->GetXaxis()->SetTitle("+Y");
    hy->GetXaxis()->CenterTitle();
    hy->Draw("colz");
    gPad->Update();
    TGaxis *ax3=new TGaxis(xmin,ymax,xmax,ymax,0.5,41.5,510,"-L");
    ax3->SetLabelSize(hy->GetXaxis()->GetLabelSize());
    ax3->SetLabelFont(hy->GetXaxis()->GetLabelFont());
    ax3->SetTitle("-Y");
    ax3->CenterTitle();
    ax3->SetTitleFont(hy->GetXaxis()->GetTitleFont());
    ax3->Draw();
    ax2->Draw();

    TPaletteAxis *palety=(TPaletteAxis*)hy->GetListOfFunctions()->FindObject("palette");
    TGaxis* palety_axis=palety->GetAxis();

    palety_axis->SetLabelSize(0);
    palety_axis->SetTickSize(0);

    can->Print(Form("%s/mipVSstrip.pdf",outDir));
    TFile *file_out=new TFile(Form("%s/%s",outDir,outName),"update");
    file_out->cd();
    can->Write(0,TObject::kOverwrite);
    hx->Write(0,TObject::kOverwrite);
    hy->Write(0,TObject::kOverwrite);
    delete file_out;

    delete can;
    delete fped;
    delete fmip;
    return 0;
}

int draw_mapping(const char* pardir,const char* filename,const char* outDir,unsigned int max=600000)
{
    int id8_pos[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id8_neg[41]={66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,47,46,22,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,1,0};
    int id5_pos[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};
    int id5_neg[41]={88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,69,68,45,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,24,23};
	int id_empty[8]={2,21,25,44,48,67,70,89};
	
    TFile* file=new TFile(Form("%s/%s",pardir,filename),"update");
    TTree *tree_in=(TTree*)file->Get("PSD");

    int nevents=tree_in->GetEntries();
    int maxentry;
    if(nevents < max){
        maxentry = nevents;
    }
    else{
        maxentry = max;
    }
    //printf("event_no = %d\n",tree_in->GetEntries());

    TCanvas *canvans=new TCanvas("canvas","canvas",600,500);
    //TH2F *h8=new TH2F("h8","h8",500,0,5000,500,0,5000);
    TH2F *h5_xpos=new TH2F("h5_xpos","h5_xpos",200,0,1000,3200,0,16000);
    h5_xpos->GetXaxis()->SetTitle("Dy5");h5_xpos->GetYaxis()->SetTitle("Dy8");
    TH2F *h5_xneg=new TH2F("h5_xneg","h5_xneg",200,0,1000,3200,0,16000);
    h5_xneg->GetXaxis()->SetTitle("Dy5");h5_xneg->GetYaxis()->SetTitle("Dy8");
    TH2F *h5_ypos=new TH2F("h5_ypos","h5_ypos",200,0,1000,3200,0,16000);
    h5_ypos->GetXaxis()->SetTitle("Dy5");h5_ypos->GetYaxis()->SetTitle("Dy8");
    TH2F *h5_yneg=new TH2F("h5_yneg","h5_yneg",200,0,1000,3200,0,16000);
    h5_yneg->GetXaxis()->SetTitle("Dy5");h5_yneg->GetYaxis()->SetTitle("Dy8");
/*
    canvans->Print(Form("%s/dy8_match.pdf[",pardir));

    for(int i=0;i<41;i++){
        h8->SetTitle(Form("X%d",i+1));
        tree_in->Draw(Form("xpos[%d]:xneg[%d] >> h8",id8_pos[i],id8_neg[40-i]),"","",maxentry);
        canvans->Print(Form("%s/dy8_match.pdf",pardir));
        h8->SetTitle(Form("Y%d",41-i));
        tree_in->Draw(Form("ypos[%d]:yneg[%d] >> h8",id8_pos[i],id8_neg[40-i]),"","",maxentry);
        canvans->Print(Form("%s/dy8_match.pdf",pardir));
    }
    canvans->Print(Form("%s/dy8_match.pdf]",pardir));
*/
    canvans->Print(Form("%s/dy5_match.pdf[",outDir));
    for(int i=0;i<41;i++){
        h5_xpos->SetName(Form("hxpos_dy58_%d",i+1));
        h5_xpos->SetTitle(Form("X%d_pos",i+1));
        tree_in->Draw(Form("xpos[%d]:xpos[%d] >> hxpos_dy58_%d",id8_pos[i],id5_pos[i],i+1),"","",maxentry);
        canvans->Print(Form("%s/dy5_match.pdf",outDir));
        h5_xpos->Write(0,TObject::kOverwrite);

        h5_xneg->SetName(Form("hxneg_dy58_%d",i+1));
        h5_xneg->SetTitle(Form("X%d_neg",i+1));
        tree_in->Draw(Form("xneg[%d]:xneg[%d] >> hxneg_dy58_%d",id8_pos[40-i],id5_pos[40-i],i+1),"","",maxentry);
        canvans->Print(Form("%s/dy5_match.pdf",outDir));
        h5_xneg->Write(0,TObject::kOverwrite);

        h5_ypos->SetName(Form("hypos_dy58_%d",i+1));
        h5_ypos->SetTitle(Form("Y%d_pos",i+1));
        tree_in->Draw(Form("ypos[%d]:ypos[%d] >> hypos_dy58_%d",id8_pos[40-i],id5_pos[40-i],i+1),"","",maxentry);
        canvans->Print(Form("%s/dy5_match.pdf",outDir));
        h5_ypos->Write(0,TObject::kOverwrite);

        h5_yneg->SetName(Form("hyneg_dy58_%d",i+1));
        h5_yneg->SetTitle(Form("Y%d_neg",i+1));
        tree_in->Draw(Form("yneg[%d]:yneg[%d] >> hyneg_dy58_%d",id8_pos[i],id5_pos[i],i+1),"","",maxentry);
        canvans->Print(Form("%s/dy5_match.pdf",outDir));
        h5_yneg->Write(0,TObject::kOverwrite);
    }
    canvans->Print(Form("%s/dy5_match.pdf]",outDir));


    delete canvans;
    delete file;

    return 0;
}

TF1* linear_fit(TProfile* hprofile,Double_t *fitrange)
{
    char FunName[100];
    sprintf(FunName,"Fitfcn_%s",hprofile->GetName());
    TF1* ffitold=(TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
    if(ffitold) delete ffitold;

    TF1* ffit=new TF1(FunName,"pol1",fitrange[0],fitrange[1]);
    hprofile->Fit(FunName,"R0");

    return ffit;
}

int fit_dy58(const char* infile,const char* pedfile,const char* outdir,const char* outfile,int pedcut=5,float range=200.0)
{
    //--- pedestal ----------------
    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};

    int channel;
    float mean,sigma;
    float xpedmean_dy8_pos[41],xpedmean_dy5_pos[41],xpedmean_dy8_neg[41],xpedmean_dy5_neg[41];
    float ypedmean_dy8_pos[41],ypedmean_dy5_pos[41],ypedmean_dy8_neg[41],ypedmean_dy5_neg[41];
    float xpedsigma_dy8_pos[41],xpedsigma_dy5_pos[41],xpedsigma_dy8_neg[41],xpedsigma_dy5_neg[41];
    float ypedsigma_dy8_pos[41],ypedsigma_dy5_pos[41],ypedsigma_dy8_neg[41],ypedsigma_dy5_neg[41];
    TFile *f_ped=new TFile(pedfile);

    TTree *tree_ped=(TTree*)f_ped->Get("xpos_ped");
    tree_ped->SetBranchAddress("channel",&channel);
    tree_ped->SetBranchAddress("mean",&mean);
    tree_ped->SetBranchAddress("sigma",&sigma);
    tree_ped->BuildIndex("channel");
    for(int i=0;i<41;i++){
        tree_ped->GetEntryWithIndex(id8[i]+1);
        xpedmean_dy8_pos[i]=mean;
        xpedsigma_dy8_pos[i]=sigma;

        tree_ped->GetEntryWithIndex(id5[i]+1);
        xpedmean_dy5_pos[i]=mean;
        xpedsigma_dy5_pos[i]=sigma;
    }

    tree_ped=(TTree*)f_ped->Get("xneg_ped");
    tree_ped->SetBranchAddress("channel",&channel);
    tree_ped->SetBranchAddress("mean",&mean);
    tree_ped->SetBranchAddress("sigma",&sigma);
    tree_ped->BuildIndex("channel");
    for(int i=0;i<41;i++){
        tree_ped->GetEntryWithIndex(id8[40-i]+1);
        xpedmean_dy8_neg[i]=mean;
        xpedsigma_dy8_neg[i]=sigma;

        tree_ped->GetEntryWithIndex(id5[40-i]+1);
        xpedmean_dy5_neg[i]=mean;
        xpedsigma_dy5_neg[i]=sigma;
    }

    tree_ped=(TTree*)f_ped->Get("ypos_ped");
    tree_ped->SetBranchAddress("channel",&channel);
    tree_ped->SetBranchAddress("mean",&mean);
    tree_ped->SetBranchAddress("sigma",&sigma);
    tree_ped->BuildIndex("channel");
    for(int i=0;i<41;i++){
        tree_ped->GetEntryWithIndex(id8[40-i]+1);
        ypedmean_dy8_pos[i]=mean;
        ypedsigma_dy8_pos[i]=sigma;

        tree_ped->GetEntryWithIndex(id5[40-i]+1);
        ypedmean_dy5_pos[i]=mean;
        ypedsigma_dy5_pos[i]=sigma;
    }

    tree_ped=(TTree*)f_ped->Get("yneg_ped");
    tree_ped->SetBranchAddress("channel",&channel);
    tree_ped->SetBranchAddress("mean",&mean);
    tree_ped->SetBranchAddress("sigma",&sigma);
    tree_ped->BuildIndex("channel");
    for(int i=0;i<41;i++){
        tree_ped->GetEntryWithIndex(id8[i]+1);
        ypedmean_dy8_neg[i]=mean;
        ypedsigma_dy8_neg[i]=sigma;

        tree_ped->GetEntryWithIndex(id5[i]+1);
        ypedmean_dy5_neg[i]=mean;
        ypedsigma_dy5_neg[i]=sigma;
    }
    delete f_ped;

    //----get profile of dy58,and fit---------------------------
    Double_t fitrange[2];
    TFile *f_in=new TFile(infile,"update");
    TH1F* hdy58_dist=new TH1F("hdy58_dist","Dy58 Ratio Distribution",200,40.0,60.0);
    TH2F* hdy58;
    TProfile* hdy58_pfx;
    TF1* lffit;
    Float_t dy58_ratio;

    FILE* fp=fopen(Form("%s/%s.txt",outdir,outfile),"w");
    fprintf(fp,"dy58_ratio:\n");
    fprintf(fp,"index\txpos\txneg\typos\tyneg\n");

    TCanvas* can=new TCanvas("can_dy58_fit","can_dy58_fit",900,900);
    can->Divide(2,2);
    can->Print(Form("%s/%s.pdf[",outdir,outfile));

    for(int i=0;i<41;i++){
        fprintf(fp,"%d\t",i+1);

        can->cd(1);
        hdy58=(TH2F*)f_in->Get(Form("hxpos_dy58_%d",i+1));
        hdy58_pfx=hdy58->ProfileX();
        fitrange[0]=xpedmean_dy5_pos[i]+pedcut*xpedsigma_dy5_pos[i];
        fitrange[1]=xpedmean_dy5_pos[i]+pedcut*xpedsigma_dy5_pos[i]+range;
        lffit=linear_fit(hdy58_pfx,fitrange);
        dy58_ratio=lffit->GetParameter(1);
        fprintf(fp,"%.3f\t",dy58_ratio);
        hdy58_dist->Fill(dy58_ratio);
        lffit->SetLineColor(kRed);
        hdy58->Draw();
        lffit->Draw("same");

        can->cd(2);
        hdy58=(TH2F*)f_in->Get(Form("hxneg_dy58_%d",i+1));
        hdy58_pfx=hdy58->ProfileX();
        fitrange[0]=xpedmean_dy5_neg[i]+pedcut*xpedsigma_dy5_neg[i];
        fitrange[1]=xpedmean_dy5_neg[i]+pedcut*xpedsigma_dy5_neg[i]+range;
        lffit=linear_fit(hdy58_pfx,fitrange);
        dy58_ratio=lffit->GetParameter(1);
        fprintf(fp,"%.3f\t",dy58_ratio);
        hdy58_dist->Fill(dy58_ratio);
        lffit->SetLineColor(kRed);
        hdy58->Draw();
        lffit->Draw("same");

        can->cd(3);
        hdy58=(TH2F*)f_in->Get(Form("hypos_dy58_%d",i+1));
        hdy58_pfx=hdy58->ProfileX();
        fitrange[0]=ypedmean_dy5_pos[i]+pedcut*ypedsigma_dy5_pos[i];
        fitrange[1]=ypedmean_dy5_pos[i]+pedcut*ypedsigma_dy5_pos[i]+range;
        lffit=linear_fit(hdy58_pfx,fitrange);
        dy58_ratio=lffit->GetParameter(1);
        fprintf(fp,"%.3f\t",dy58_ratio);
        hdy58_dist->Fill(dy58_ratio);
        lffit->SetLineColor(kRed);
        hdy58->Draw();
        lffit->Draw("same");

        can->cd(4);
        hdy58=(TH2F*)f_in->Get(Form("hyneg_dy58_%d",i+1));
        hdy58_pfx=hdy58->ProfileX();
        fitrange[0]=ypedmean_dy5_neg[i]+pedcut*ypedsigma_dy5_neg[i];
        fitrange[1]=ypedmean_dy5_neg[i]+pedcut*ypedsigma_dy5_neg[i]+range;
        lffit=linear_fit(hdy58_pfx,fitrange);
        dy58_ratio=lffit->GetParameter(1);
        fprintf(fp,"%.3f\n",dy58_ratio);
        hdy58_dist->Fill(dy58_ratio);
        lffit->SetLineColor(kRed);
        hdy58->Draw();
        lffit->Draw("same");

        can->Print(Form("%s/%s.pdf",outdir,outfile));
    }
    can->Print(Form("%s/%s.pdf]",outdir,outfile));

    hdy58_dist->Write(0,TObject::kOverwrite);
    fclose(fp);
    delete f_in;
    delete can;

    return 0;
}

int draw_channels(const char* pardir,const char* filename,const char* outDir,const char* outName)
{
    float xpos_mean[90]={351.177,311.881,281.001,336.517,258.712,313.785,367.483,306.686,365.061,400.932,389.930,414.645,381.898,411.290,457.176,472.376,463.821,539.654,495.049,496.161,441.246,514.411,469.810,480.274,360.767,386.358,352.114,449.317,322.460,344.191,420.976,343.341,240.283,314.475,384.396,319.270,298.015,267.687,222.540,274.143,163.888,253.895,150.799,200.748,175.790,145.786,524.328,515.608,526.098,480.832,485.753,494.254,464.564,435.246,343.468,312.436,341.823,331.957,229.788,180.386,167.866,143.693,157.214,148.327,113.189,226.599,165.258,195.581,410.142,399.905,425.209,444.228,467.265,398.202,483.074,502.183,482.189,552.946,492.655,465.126,451.151,434.841,447.516,462.986,427.054,431.159,441.765,404.810,370.627,377.385};
    float xneg_mean[90]={446.400,520.242,514.340,474.583,496.813,511.273,451.954,431.075,434.599,440.726,489.202,435.936,499.171,509.421,467.255,535.061,496.474,563.510,574.983,541.864,537.884,496.657,596.976,366.517,391.959,375.442,351.044,405.371,442.040,453.049,475.211,478.180,435.166,421.951,420.496,458.119,457.642,411.283,510.404,547.936,484.632,467.074,450.268,434.619,498.571,393.499,420.707,422.288,433.151,393.430,429.871,415.825,436.522,474.949,460.671,484.791,485.606,453.007,471.585,500.918,497.301,502.769,396.582,485.579,465.391,514.925,499.670,426.066,424.271,463.562,454.416,464.023,409.735,455.579,354.684,375.137,407.772,453.645,461.221,462.112,465.888,447.161,534.270,467.812,467.011,472.202,423.727,475.979,441.521,414.418};
    float ypos_mean[90]={452.297,396.719,465.625,479.740,528.990,451.267,480.754,451.002,424.205,458.843,438.643,442.601,436.487,414.472,413.737,369.849,492.507,466.578,337.399,441.559,374.001,292.447,310.534,376.861,384.398,464.878,458.101,409.130,432.157,468.180,525.836,426.869,496.900,449.711,540.653,521.539,493.202,489.634,462.673,486.457,399.784,497.665,358.482,479.102,343.498,421.276,399.754,497.368,455.586,454.538,413.200,431.225,420.818,345.452,323.384,384.369,422.946,351.598,326.737,311.787,384.763,351.652,284.502,269.942,281.036,241.743,251.200,222.739,462.892,475.289,442.673,449.252,420.579,475.795,479.726,412.265,460.324,410.219,446.704,432.104,474.956,529.460,500.751,453.041,417.644,420.573,377.888,386.825,331.463,392.985};
    float yneg_mean[90]={417.906,472.811,455.049,452.441,453.630,426.665,413.650,458.842,491.942,415.372,438.802,455.469,401.188,413.670,430.686,444.817,489.856,400.641,473.606,387.638,458.840,453.257,439.884,315.324,256.090,265.857,255.013,298.694,286.077,283.444,285.370,303.172,315.953,275.825,234.171,293.651,249.059,276.830,270.094,341.402,282.482,244.047,257.179,300.565,285.312,233.734,427.574,432.091,328.272,419.002,395.342,355.594,297.112,425.499,384.624,395.578,408.127,322.187,396.572,409.241,430.655,393.992,374.721,431.261,373.104,342.076,386.307,329.733,499.101,476.227,497.179,514.796,493.029,511.467,495.702,494.944,478.117,528.368,518.425,441.704,552.574,457.136,484.718,451.821,499.494,413.811,472.878,330.334,307.009,269.479};

    float xpos_sigma[90]={7.186,6.971,1.921,6.608,6.890,6.466,6.866,6.532,6.517,6.399,6.495,6.736,6.680,6.601,6.576,6.403,6.517,6.354,6.584,6.436,6.515,2.087,6.424,6.550,6.334,2.325,6.498,6.272,6.350,6.442,6.311,6.185,6.248,6.107,6.220,6.225,6.167,6.207,6.184,6.088,6.138,6.183,6.122,5.988,2.209,6.111,6.611,6.220,2.247,6.754,6.473,6.740,6.371,6.809,6.513,6.743,6.442,6.705,6.470,7.167,6.596,6.872,6.700,7.551,6.528,7.843,10.055,2.087,6.243,6.128,2.049,6.289,6.143,6.184,6.144,6.031,6.181,6.202,6.191,6.021,6.215,6.181,6.096,6.222,6.132,6.091,6.216,6.188,6.161,2.080};
    float xneg_sigma[90]={7.258,4.987,2.121,4.287,6.122,3.680,4.568,4.003,4.629,3.939,4.070,3.956,4.271,3.629,4.322,3.798,4.123,3.764,4.256,3.939,4.276,2.339,3.906,4.058,3.844,2.200,3.981,4.007,3.865,3.975,4.005,4.152,3.950,3.892,3.830,3.926,3.873,3.916,3.567,4.027,3.688,3.869,3.824,3.871,2.012,3.662,3.840,3.820,1.905,3.924,3.690,4.091,3.875,3.930,3.780,3.888,3.772,4.040,3.761,3.832,4.007,3.941,4.112,4.074,4.092,4.524,4.385,2.165,3.807,3.919,1.954,3.832,3.858,3.813,3.888,3.853,3.759,3.819,3.864,3.873,3.881,3.933,3.821,3.935,3.894,3.724,3.922,4.022,3.869,2.067};
    float ypos_sigma[90]={6.189,6.342,2.350,6.030,6.732,6.012,6.314,5.739,6.390,5.335,5.762,5.607,5.814,5.617,5.837,5.541,5.629,5.484,5.699,5.428,5.683,1.881,5.564,5.553,5.539,2.220,5.417,5.404,5.365,5.449,5.327,5.227,5.271,5.372,5.308,5.241,5.359,5.287,5.248,5.226,5.356,5.271,5.243,5.296,1.908,5.336,5.615,5.409,2.077,5.731,5.560,5.558,5.503,5.528,5.338,6.021,5.660,6.206,5.443,6.132,5.899,7.639,6.225,6.751,6.093,6.532,6.155,1.884,5.385,5.328,1.763,5.360,5.378,5.348,5.390,5.330,5.265,5.261,5.454,5.268,5.323,5.400,5.180,5.424,5.330,5.304,5.402,5.329,5.376,1.715};
    float yneg_sigma[90]={6.251,6.709,2.311,6.073,7.180,6.005,7.595,5.919,6.332,5.667,6.053,5.759,5.986,5.689,5.805,5.514,5.795,5.675,5.825,5.724,5.772,2.028,5.733,5.656,5.737,2.277,5.668,5.813,5.579,5.704,5.587,5.525,5.518,5.692,5.598,5.619,5.543,5.573,5.521,5.534,5.430,5.531,5.503,5.619,2.133,5.534,5.749,5.618,2.027,5.963,5.713,5.758,5.704,5.854,5.561,6.067,5.680,5.911,5.828,6.276,5.836,6.449,5.981,6.281,5.897,6.339,6.067,1.879,5.761,5.650,2.164,5.546,5.851,5.519,5.993,5.689,5.701,5.764,5.561,5.649,5.621,5.635,5.709,5.833,5.681,5.842,5.721,5.698,5.878,1.956};

    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};

    TFile* file=new TFile(Form("%s/%s",pardir,filename));
    int channel;
    float ped,sigma;
    char label[4][5]={"xpos","xneg","ypos","yneg"};
    TFile* file_ped=new TFile(Form("%s/new_ped.root",outDir),"recreate");
    TTree* tree_ped[4];
    for(int i=0;i<4;i++){
        tree_ped[i]=new TTree(Form("%s_ped",label[i]),Form("%s_ped",label[i]));
        tree_ped[i]->Branch("channel",&channel,"channel/I");
        tree_ped[i]->Branch("mean",&ped,"mean/F");
        tree_ped[i]->Branch("sigma",&sigma,"sigma/F");
    }

    TCanvas* xcan=new TCanvas("xcan","xcan",800,900);
    TCanvas* ycan=new TCanvas("ycan","ycan",800,900);
    xcan->Divide(2,2);
    ycan->Divide(2,2);
    for(int i=0;i<4;i++){
        xcan->cd(i+1);gPad->SetLogy();
        ycan->cd(i+1);gPad->SetLogy();
    }
    xcan->Print(Form("%s/xmips.pdf[",outDir));
    ycan->Print(Form("%s/ymips.pdf[",outDir));
    /*
    FILE* fp_xpos=fopen(Form("%s/xpos_mip.csv",pardir),"w");
    FILE* fp_xneg=fopen(Form("%s/xneg_mip.csv",pardir),"w");
    FILE* fp_ypos=fopen(Form("%s/ypos_mip.csv",pardir),"w");
    FILE* fp_yneg=fopen(Form("%s/yneg_mip.csv",pardir),"w");
    */
    //FILE* fp_mip=fopen(Form("%s/mips.csv",outDir),"w");
    //fprintf(fp_mip,"channel,xpos,xneg,ypos,yneg\n");
    FILE* fp_mip2=fopen(Form("%s/mips_2.csv",outDir),"w");
    fprintf(fp_mip2,"channel,xpos,xneg,ypos,yneg\n");

    TH1F *hxpos8[41],*hxpos5[41],*hxneg8[41],*hxneg5[41];
    TH1F *hypos8[41],*hypos5[41],*hyneg8[41],*hyneg5[41];
    /*
    TTree* tree_in=(TTree*)file->Get("PSD");
    Int_t begin=130001;
    for(int ch_id=0;ch_id<41;ch_id++){
        hxpos8[41]=new TH1F(Form("xpos_%d_dy8",ch_id+1),Form("xpos_%d_dy8",ch_id+1),3000,100,3100);
        tree_in->Project(Form("xpos_%d_dy8",ch_id+1),Form("xpos[%d]",id8[ch_id]),"","",1000000000,begin);
        hxpos5[41]=new TH1F(Form("xpos_%d_dy5",ch_id+1),Form("xpos_%d_dy5",ch_id+1),800,100,900);
        tree_in->Project(Form("xpos_%d_dy5",ch_id+1),Form("xpos[%d]",id5[ch_id]),"","",1000000000,begin);
        hxneg8[41]=new TH1F(Form("xneg_%d_dy8",ch_id+1),Form("xneg_%d_dy8",ch_id+1),3000,100,3100);
        tree_in->Project(Form("xneg_%d_dy8",ch_id+1),Form("xneg[%d]",id8[40-ch_id]),"","",1000000000,begin);
        hxneg5[41]=new TH1F(Form("xneg_%d_dy5",ch_id+1),Form("xneg_%d_dy5",ch_id+1),800,100,900);
        tree_in->Project(Form("xneg_%d_dy5",ch_id+1),Form("xneg[%d]",id5[40-ch_id]),"","",1000000000,begin);
        hypos8[41]=new TH1F(Form("ypos_%d_dy8",ch_id+1),Form("ypos_%d_dy8",ch_id+1),3000,100,3100);
        tree_in->Project(Form("ypos_%d_dy8",ch_id+1),Form("ypos[%d]",id8[40-ch_id]),"","",1000000000,begin);
        hypos5[41]=new TH1F(Form("ypos_%d_dy5",ch_id+1),Form("ypos_%d_dy5",ch_id+1),800,100,900);
        tree_in->Project(Form("ypos_%d_dy5",ch_id+1),Form("ypos[%d]",id5[40-ch_id]),"","",1000000000,begin);
        hyneg8[41]=new TH1F(Form("yneg_%d_dy8",ch_id+1),Form("yneg_%d_dy8",ch_id+1),3000,100,3100);
        tree_in->Project(Form("yneg_%d_dy8",ch_id+1),Form("yneg[%d]",id8[ch_id]),"","",1000000000,begin);
        hyneg5[41]=new TH1F(Form("yneg_%d_dy5",ch_id+1),Form("yneg_%d_dy5",ch_id+1),800,100,900);
        tree_in->Project(Form("yneg_%d_dy5",ch_id+1),Form("yneg[%d]",id5[ch_id]),"","",1000000000,begin);
    }
    */
    TF1* ffit;
    float mpv,fwhm;
    for(int ch_id=0;ch_id<41;ch_id++){
        //fprintf(fp_mip,"%d,",ch_id+1);
        fprintf(fp_mip2,"%d,",ch_id+1);

        hxpos8[ch_id]=(TH1F*)file->Get(Form("xpos_%d",id8[ch_id]+1));
        hxpos8[ch_id]->Rebin(10);
        hxpos8[ch_id]->SetTitle(Form("X_layer,Strip_%d,pos_dy8",ch_id+1));
        ffit=langaus(hxpos8[ch_id],mpv,fwhm,xpos_mean[id8[ch_id]],xpos_sigma[id8[ch_id]]);
        xcan->cd(1);
        hxpos8[ch_id]->Draw();
        ffit->Draw("lsame");
        //fprintf(fp_mip,"%.2f,",(mpv-xpos_mean[id8[ch_id]]));
        ffit=pedfit(hxpos8[ch_id],xpos_mean[id8[ch_id]],xpos_sigma[id8[ch_id]],ped,sigma);
        ffit->Draw("lsame");
         fprintf(fp_mip2,"%.2f,",(mpv-ped));
         channel=id8[ch_id]+1;
         tree_ped[0]->Fill();

        hxpos5[ch_id]=(TH1F*)file->Get(Form("xpos_%d",id5[ch_id]+1));
        hxpos5[ch_id]->Rebin(10);
        hxpos5[ch_id]->SetTitle(Form("X_layer,Strip_%d,pos_dy5",ch_id+1));
        ffit=pedfit(hxpos5[ch_id],xpos_mean[id5[ch_id]],xpos_sigma[id5[ch_id]],ped,sigma);
        channel=id5[ch_id]+1;
        tree_ped[0]->Fill();
        xcan->cd(2);
        hxpos5[ch_id]->Draw();
        ffit->Draw("lsame");

        hxneg8[ch_id]=(TH1F*)file->Get(Form("xneg_%d",id8[40-ch_id]+1));
        hxneg8[ch_id]->Rebin(10);
        hxneg8[ch_id]->SetTitle(Form("X_layer,Strip_%d,neg_dy8",ch_id+1));
        ffit=langaus(hxneg8[ch_id],mpv,fwhm,xneg_mean[id8[40-ch_id]],xneg_sigma[id8[40-ch_id]]);
        xcan->cd(3);
        hxneg8[ch_id]->Draw();
        ffit->Draw("lsame");
        //fprintf(fp_mip,"%.2f,",(mpv-xneg_mean[id8[40-ch_id]]));
        ffit=pedfit(hxneg8[ch_id],xneg_mean[id8[40-ch_id]],xneg_sigma[id8[40-ch_id]],ped,sigma);
        ffit->Draw("lsame");
         fprintf(fp_mip2,"%.2f,",(mpv-ped));
         channel=id8[40-ch_id]+1;
         tree_ped[1]->Fill();

        hxneg5[ch_id]=(TH1F*)file->Get(Form("xneg_%d",id5[40-ch_id]+1));
        hxneg5[ch_id]->Rebin(10);
        hxneg5[ch_id]->SetTitle(Form("X_layer,Strip_%d,neg_dy5",ch_id+1));
        ffit=pedfit(hxneg5[ch_id],xneg_mean[id5[40-ch_id]],xneg_sigma[id5[40-ch_id]],ped,sigma);
        channel=id5[40-ch_id]+1;
        tree_ped[1]->Fill();
        xcan->cd(4);
        hxneg5[ch_id]->Draw();
        ffit->Draw("lsame");
        xcan->Print(Form("%s/xmips.pdf",outDir));

        hypos8[ch_id]=(TH1F*)file->Get(Form("ypos_%d",id8[40-ch_id]+1));
        hypos8[ch_id]->Rebin(10);
        hypos8[ch_id]->SetTitle(Form("Y_layer,Strip_%d,pos_dy8",ch_id+1));
        ffit=langaus(hypos8[ch_id],mpv,fwhm,ypos_mean[id8[40-ch_id]],ypos_sigma[id8[40-ch_id]]);
        ycan->cd(1);
        hypos8[ch_id]->Draw();
        ffit->Draw("lsame");
        //fprintf(fp_mip,"%.2f,",(mpv-ypos_mean[id8[40-ch_id]]));
        ffit=pedfit(hypos8[ch_id],ypos_mean[id8[40-ch_id]],ypos_sigma[id8[40-ch_id]],ped,sigma);
        ffit->Draw("lsame");
         fprintf(fp_mip2,"%.2f,",(mpv-ped));
         channel=id8[40-ch_id]+1;
         tree_ped[2]->Fill();

        hypos5[ch_id]=(TH1F*)file->Get(Form("ypos_%d",id5[40-ch_id]+1));
        hypos5[ch_id]->Rebin(10);
        hypos5[ch_id]->SetTitle(Form("Y_layer,Strip_%d,pos_dy5",ch_id+1));
       ffit= pedfit(hypos5[ch_id],ypos_mean[id5[40-ch_id]],ypos_sigma[id5[40-ch_id]],ped,sigma);
        channel=id5[40-ch_id]+1;
        tree_ped[2]->Fill();
        ycan->cd(2);
        hypos5[ch_id]->Draw();
        ffit->Draw("lsame");

        hyneg8[ch_id]=(TH1F*)file->Get(Form("yneg_%d",id8[ch_id]+1));
        hyneg8[ch_id]->Rebin(10);
        hyneg8[ch_id]->SetTitle(Form("Y_layer,Strip_%d,neg_dy8",ch_id+1));
        ffit=langaus(hyneg8[ch_id],mpv,fwhm,yneg_mean[id8[ch_id]],yneg_sigma[id8[ch_id]]);
        ycan->cd(3);
        hyneg8[ch_id]->Draw();
        ffit->Draw("lsame");
        //fprintf(fp_mip,"%.2f\n",(mpv-yneg_mean[id8[ch_id]]));
        ffit=pedfit(hyneg8[ch_id],yneg_mean[id8[ch_id]],yneg_sigma[id8[ch_id]],ped,sigma);
        ffit->Draw("lsame");
         fprintf(fp_mip2,"%.2f\n",(mpv-ped));
         channel=id8[ch_id]+1;
         tree_ped[3]->Fill();

        hyneg5[ch_id]=(TH1F*)file->Get(Form("yneg_%d",id5[ch_id]+1));
        hyneg5[ch_id]->Rebin(10);
        hyneg5[ch_id]->SetTitle(Form("Y_layer,Strip_%d,neg_dy5",ch_id+1));
        ffit=pedfit(hyneg5[ch_id],yneg_mean[id5[ch_id]],yneg_sigma[id5[ch_id]],ped,sigma);
        channel=id5[ch_id]+1;
        tree_ped[3]->Fill();
        ycan->cd(4);
        hyneg5[ch_id]->Draw();
        ffit->Draw("lsame");
        ycan->Print(Form("%s/ymips.pdf",outDir));
    }

    xcan->Print(Form("%s/xmips.pdf]",outDir));
    ycan->Print(Form("%s/ymips.pdf]",outDir));

    TFile* file_out=new TFile(Form("%s/%s",outDir,outName),"update");
    for(int i=0;i<41;i++){
        hxpos8[i]->Write(0,TObject::kOverwrite);
        hxneg5[i]->Write(0,TObject::kOverwrite);
        hxpos8[i]->Write(0,TObject::kOverwrite);
        hxpos8[i]->Write(0,TObject::kOverwrite);

        hypos8[i]->Write(0,TObject::kOverwrite);
        hyneg5[i]->Write(0,TObject::kOverwrite);
        hypos8[i]->Write(0,TObject::kOverwrite);
        hypos8[i]->Write(0,TObject::kOverwrite);
    }

    file_ped->cd();
    for(int i=0;i<4;i++){
        tree_ped[i]->Write();
    }
    delete file_ped;
    delete file_out;
    delete xcan;
    delete ycan;
    delete file;

    //fclose(fp_mip);
    fclose(fp_mip2);

    return 0;
}

//-----end---------------------------------
/****************END**********************************/

int select_file(const char* dir,char *filename)
{
    Int_t selected_id;
    char* dirs[100];
    Int_t couter=0;
    TString full_pathname;
    full_pathname=gSystem->ExpandPathName(dir);
    TSystemDirectory pedestal_dir(full_pathname.Data(),full_pathname.Data());
    TList *files=pedestal_dir.GetListOfFiles();
    if(files){
        if(files->GetSize()!=3){
            printf("\t#There are more than 1 raw data files.\n\tPlease select the one you wanna to analyze:\n");
            files->Sort();
            TSystemFile *file;
            TIter next(files);

            while((file=(TSystemFile*)next())){
                if(!file->IsDirectory()){
                    couter++;
                    printf("\t%d) %s\n",couter,file->GetName());
                    char* pdir=new char[200];
                    sprintf(pdir,"%s",file->GetName());
                    dirs[couter-1]=pdir;
                }
            }
            /*
            for(int i=0;i<couter;i++){
                printf("%d. %s\n",i+1,dirs[i]);
            }
            */

            scanf("%d",&selected_id);
            printf("\t#You have selected file: %s\n",dirs[selected_id-1]);
            sprintf(filename,"%s",dirs[selected_id-1]);
        }
        else{
            files->Sort();
            TSystemFile *file;
            TIter next(files);
            while((file=(TSystemFile*)next())){
                if(!file->IsDirectory()){
                    sprintf(filename,"%s",file->GetName());
                }
            }

        }
    }

    for(int i=0;i<couter;i++){
        delete [] (dirs[i]);
    }

    return 0;
}

int getconfig_lowthresh(const char* infile,const char* outdir,int level=3)
{
    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};
    char label[][10]={"xpos","xneg","ypos","yneg"};
    const char fee_id[4]={0x20,0x24,0x28,0x2c};
    const int limit=8*0xFF-0x200;

    char outfile[200];
    char outbuffer[2];
    TFile* file_in=new TFile(infile);
    TTree* tree_in=0;
    std::ofstream* out=0;

    Float_t mean,sigma;
    Int_t channel;
    for(int i=0;i<4;i++){
        tree_in=(TTree*)file_in->Get(Form("%s_ped",label[i]));
        tree_in->BuildIndex("channel");

        tree_in->SetBranchAddress("channel",&channel);
        tree_in->SetBranchAddress("sigma",&sigma);
        tree_in->SetBranchAddress("mean",&mean);

        sprintf(outfile,"%s/%s_lowthresh.bin",outdir,label[i]);
        out=new std::ofstream(outfile,std::ios_base::binary);
        out->write(fee_id+i,1);
        for(int j=0;j<90;j++){
            tree_in->GetEntryWithIndex(j+1);

            outbuffer[0]=j&0xFF;
            if((mean+level*sigma) > (limit+7)){
                outbuffer[1]=0xFF;
                printf("warning: %s(fee_id=%x),channel_%d low threshold larger than the limit,it will be set to 0xFF\n",label[i],fee_id[i],i+1);
            }
            else{
                int tmp=(mean+level*sigma+0x200)/8;
                outbuffer[1]=tmp&0xFF;
            }
            out->write(outbuffer,2);
        }
        out->close();
        delete out;
    }
    delete file_in;

    return 0;
}

int getconfig_highthresh(const char* outdir,char thresh=0x0)
{
    char label[][10]={"xpos","xneg","ypos","yneg"};
    const char fee_id[4]={0x20,0x24,0x28,0x2c};

    char outfile[200];
    char outbuffer[2];
    std::ofstream* out=0;

    for(int i=0;i<4;i++){
        sprintf(outfile,"%s/%s_highthresh.bin",outdir,label[i]);
        out=new std::ofstream(outfile,std::ios_base::binary);
        out->write(fee_id+i,1);
        for(int j=0;j<90;j++){
            outbuffer[0]=j&0xFF;
            outbuffer[1]=thresh;
            out->write(outbuffer,2);
        }
        out->close();
        delete out;
    }

    return 0;
}

int print_pedestal(const char* infile,const char* outdir)
{
    int id8[41]={0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66};
    int id5[41]={23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88};
    char label[][10]={"xpos","xneg","ypos","yneg"};

    char outfile[200];
    TFile* file_in=new TFile(infile);
    TTree* tree_in=0;
    FILE* fp_mean;
    FILE* fp_sigma;

    Float_t mean,sigma;
    Int_t channel;
    for(int i=0;i<4;i++){
        tree_in=(TTree*)file_in->Get(Form("%s_ped",label[i]));
        tree_in->BuildIndex("channel");

        tree_in->SetBranchAddress("channel",&channel);
        tree_in->SetBranchAddress("sigma",&sigma);
        tree_in->SetBranchAddress("mean",&mean);

        sprintf(outfile,"%s/%s_mean.txt",outdir,label[i]);
        fp_mean=fopen(outfile,"w");
        sprintf(outfile,"%s/%s_sigma.txt",outdir,label[i]);
        fp_sigma=fopen(outfile,"w");
        for(int j=0;j<90;j++){
            tree_in->GetEntryWithIndex(j+1);

            if(j==89){
                fprintf(fp_mean,"%.3f",mean);
                fprintf(fp_sigma,"%.3f",sigma);
            }
            else{
                fprintf(fp_mean,"%.3f,",mean);
                fprintf(fp_sigma,"%.3f,",sigma);
            }
        }
        fclose(fp_mean);
        fclose(fp_sigma);
    }

    delete file_in;

    return 0;
}
/*
int analyze_mip_normal(const char* date,const char* testdir,const char* out="analyze")
{
    char rawdir[200];
    char analyzedir[200];
    char dir_analyze_mip[200];
    char ped_analyze_dir[200];

    sprintf(rawdir,"%s/raw",testdir);
    sprintf(analyzedir,"%s/%s",testdir,out);
    sprintf(dir_analyze_mip,"%s/mip",analyzedir);
    sprintf(ped_analyze_dir,"%s/analyze",testdir);

    if(!(gSystem->OpenDirectory(rawdir))){
        printf("Error! : No raw data directory under '%s'\n",testdir);
        return -1;
    }
    printf("WORKING DIRECTORY: %s\n",gSystem->WorkingDirectory());
    if(!(gSystem->OpenDirectory(analyzedir))){
        gSystem->MakeDirectory(analyzedir);
    }

    sprintf(dir_analyze_mip,"%s",Form("%s/mip",analyzedir));
    if(!(gSystem->OpenDirectory(dir_analyze_mip))){
        gSystem->MakeDirectory(dir_analyze_mip);
    }

    //--------------------------convert raw data to root files--------------------------------------
    char filename[200];
    printf("Start converting raw data files...\n");

    select_file(Form("%s/mip",rawdir),filename);
    printf("\tconvert %s: \n",Form("%s/mip/%s",rawdir,filename));
    //convert_event(Form("%s/mip",rawdir),filename,dir_analyze_mip,"raw.root");
    convert_psd_scidata(0x20,Form("%s/mip",rawdir),filename,dir_analyze_mip,"raw.root",date);
    //--------------------------------------------------------------------------------------------

    draw_channels(dir_analyze_mip,"raw.root");
    draw_mapping(dir_analyze_mip,"raw.root");
    draw_mip(Form("%s/raw.root",dir_analyze_mip),Form("%s/hv/hv.root",ped_analyze_dir));

    return 0;
}

int analyze_mip_compressed(const char* date,const char* testdir,const char* out="analyze")
{
    char rawdir[200];
    char analyzedir[200];
    char dir_analyze_mip[200];
    char ped_analyze_dir[200];

    sprintf(rawdir,"%s/raw",testdir);
    sprintf(analyzedir,"%s/%s",testdir,out);
    sprintf(dir_analyze_mip,"%s/mip",analyzedir);
    sprintf(ped_analyze_dir,"%s/analyze",testdir);

    if(!(gSystem->OpenDirectory(rawdir))){
        printf("Error! : No raw data directory under '%s'\n",testdir);
        return -1;
    }
    printf("WORKING DIRECTORY: %s\n",gSystem->WorkingDirectory());
    if(!(gSystem->OpenDirectory(analyzedir))){
        gSystem->MakeDirectory(analyzedir);
    }

    sprintf(dir_analyze_mip,"%s",Form("%s/mip",analyzedir));
    if(!(gSystem->OpenDirectory(dir_analyze_mip))){
        gSystem->MakeDirectory(dir_analyze_mip);
    }

    //--------------------------convert raw data to root files--------------------------------------
    char filename[200];
    printf("Start converting raw data files...\n");

    select_file(Form("%s/mip",rawdir),filename);
    printf("\tconvert %s: \n",Form("%s/mip/%s",rawdir,filename));
    //convert_event(Form("%s/mip",rawdir),filename,dir_analyze_mip,"raw.root");
    convert_psd_scidata(0x60,Form("%s/mip",rawdir),filename,dir_analyze_mip,"raw.root",date);
    //--------------------------------------------------------------------------------------------

    draw_channels(dir_analyze_mip,"raw.root");
    draw_mapping(dir_analyze_mip,"raw.root");
    draw_mip(Form("%s/raw.root",dir_analyze_mip),Form("%s/hv/hv.root",ped_analyze_dir));

    return 0;
}

int analyze_mip_smaller(const char* date,const char* testdir,const char* out="analyze")
{
    char rawdir[400];
    char analyzedir[400];
    char dir_analyze_mip[400];
    char ped_analyze_dir[400];

    sprintf(rawdir,"%s/raw",testdir);
    sprintf(analyzedir,"%s/%s",testdir,out);
    sprintf(dir_analyze_mip,"%s/mip",analyzedir);
    sprintf(ped_analyze_dir,"%s/analyze",testdir);

    if(!(gSystem->OpenDirectory(rawdir))){
        printf("Error! : No raw data directory under '%s'\n",testdir);
        return -1;
    }
    printf("WORKING DIRECTORY: %s\n",gSystem->WorkingDirectory());
    if(!(gSystem->OpenDirectory(analyzedir))){
        gSystem->MakeDirectory(analyzedir);
    }

    sprintf(dir_analyze_mip,"%s",Form("%s/mip",analyzedir));
    if(!(gSystem->OpenDirectory(dir_analyze_mip))){
        gSystem->MakeDirectory(dir_analyze_mip);
    }

    //--------------------------convert raw data to root files--------------------------------------
    char filename[200];
    printf("Start converting raw data files...\n");

    select_file(Form("%s/mip",rawdir),filename);
    printf("\tconvert %s: \n",Form("%s/mip/%s",rawdir,filename));
    //convert_event(Form("%s/mip",rawdir),filename,dir_analyze_mip,"raw.root");
    convert_psd_scidata(0x30,Form("%s/mip",rawdir),filename,dir_analyze_mip,"raw.root",date);
    //--------------------------------------------------------------------------------------------

    draw_channels(dir_analyze_mip,"raw.root");
    draw_mapping(dir_analyze_mip,"raw.root");
    draw_mip(Form("%s/raw.root",dir_analyze_mip),Form("%s/hv/hv.root",ped_analyze_dir));

    return 0;
}

int analyze_ped(const char* date,const char* testdir,const char* out="analyze",const char* refdir="standard_ped")
{
    char rawdir[200];
    char analyzedir[200];
    char dir_analyze_hv[200];
    char dir_analyze_nohv[200];
    char log_file[200];
    FILE* fp=0;
    sprintf(rawdir,"%s/raw",testdir);
    sprintf(analyzedir,"%s/%s",testdir,out);

    if(!(gSystem->OpenDirectory(rawdir))){
        printf("Error! : No raw data directory under '%s'\n",testdir);
        fprintf(fp,"\tError! : No raw data directory under '%s'\n",testdir);
        return -1;
    }
    printf("WORKING DIRECTORY: %s\n",gSystem->WorkingDirectory());
    if(!(gSystem->OpenDirectory(analyzedir))){
        printf("\t#making new directory: %s\n",analyzedir);
        gSystem->MakeDirectory(analyzedir);

        sprintf(dir_analyze_nohv,"%s",Form("%s/nohv",analyzedir));
        gSystem->MakeDirectory(dir_analyze_nohv);

        sprintf(dir_analyze_hv,"%s",Form("%s/hv",analyzedir));
        gSystem->MakeDirectory(dir_analyze_hv);
    }
    else{
        printf("\tWarning: %s already exists,really wanna save result into this dir?(Yes:y,No:n)\n",analyzedir);
        char flag;
        scanf("%c",&flag);
        if(flag == 'y' || flag == 'Y'){
            printf("overwriting %s\n",analyzedir);
            sprintf(dir_analyze_nohv,"%s",Form("%s/nohv",analyzedir));
            sprintf(dir_analyze_hv,"%s",Form("%s/hv",analyzedir));
        }
        else if(flag == 'n' || flag == 'N'){
            printf("Enter new output directory name:\n");
            char new_out[200];
            scanf("%s",new_out);

            sprintf(analyzedir,"%s/%s",testdir,new_out);
            printf("Making new directory: %s\n",analyzedir);
            gSystem->MakeDirectory(analyzedir);

            sprintf(dir_analyze_nohv,"%s",Form("%s/nohv",analyzedir));
            gSystem->MakeDirectory(dir_analyze_nohv);

            sprintf(dir_analyze_hv,"%s",Form("%s/hv",analyzedir));
            gSystem->MakeDirectory(dir_analyze_hv);
        }
        else{
            printf("Error: Please enter Y or N\n");
            return -1;
        }
    }
    printf("OUTPUT DIRECTORY: %s\n\n",analyzedir);
    sprintf(log_file,"%s/log_ped.txt",analyzedir);
    fp=fopen(log_file,"w");
    fprintf(fp,"Log Info of pedestal analysis:\n");
    //--------------------------convert raw data to root files--------------------------------------
    char filename[200];
    printf("Start converting raw data files...\n");
    fprintf(fp,"\tRaw data converting:\n");

    select_file(Form("%s/nohv",rawdir),filename);
    printf("\tconvert %s: \n",Form("%s/nohv/%s",rawdir,filename));
    //convert_event(Form("%s/nohv",rawdir),filename,Form("%s/nohv",analyzedir),"raw.root");
    convert_psd_scidata(0x20,Form("%s/nohv",rawdir),filename,Form("%s/nohv",analyzedir),"raw.root",date);
    fprintf(fp,"\nconvert %s to %s\n",Form("%s/nohv/%s",rawdir,filename),Form("%s/nohv/raw.root",analyzedir));

    select_file(Form("%s/hv",rawdir),filename);
    printf("\tconvert %s: \n",Form("%s/hv/%s",rawdir,filename));
    //convert_event(Form("%s/hv",rawdir),filename,Form("%s/hv",analyzedir),"raw.root");
    convert_psd_scidata(0x20,Form("%s/hv",rawdir),filename,Form("%s/hv",analyzedir),"raw.root",date);
    fprintf(fp,"\nconvert %s to %s\n",Form("%s/hv/%s",rawdir,filename),Form("%s/hv/raw.root",analyzedir));
    fprintf(fp,"\n\n");
    //---------------------------get ped fitting values----------------------------------------------
    printf("\nStart fitting pedestals...\n");
    getpedseed_event(dir_analyze_nohv,"raw.root",dir_analyze_nohv,"nohv");
    getpedseed_event(dir_analyze_hv,"raw.root",dir_analyze_hv,"hv");

    //---------------------------draw graphs---------------------------------------------------------
    fprintf(fp,"ped reference directory:\n");
    fprintf(fp,"%s\n",refdir);
    //draw_rel(Form("%s/hv.root",dir_analyze_hv),analyzedir,Form("%s/nohv.root",dir_analyze_nohv));
    draw_rel(Form("%s/nohv.root",dir_analyze_nohv),dir_analyze_nohv,Form("%s/analyze/nohv/nohv.root",refdir));
    draw_rel(Form("%s/hv.root",dir_analyze_hv),dir_analyze_hv,Form("%s/analyze/hv/hv.root",refdir));

    //draw_relp(Form("%s/hv.root",dir_analyze_hv),analyzedir,Form("%s/nohv.root",dir_analyze_nohv));
    draw_relp(Form("%s/nohv.root",dir_analyze_nohv),dir_analyze_nohv,Form("%s/analyze/nohv/nohv.root",refdir));
    draw_relp(Form("%s/hv.root",dir_analyze_hv),dir_analyze_hv,Form("%s/analyze/hv/hv.root",refdir));

    draw_ped(dir_analyze_nohv,"raw.root");
    draw_ped(dir_analyze_hv,"raw.root");

    //---------------------------
    //draw_mapping(dir_analyze_hv,"raw.root");
    //draw_mapping(dir_analyze_nohv,"raw.root");
    //------------------------------
    fclose(fp);

    return 0;
}

int analyze_calib(const char* date,const Char_t* parentDir,const Char_t* infile,const Char_t* outDir)
{
    convert_psd_scidata(0xa0,parentDir,infile,outDir,"calib_raw.root",date);
    fit_calibration(outDir,"calib_raw.root",outDir,"calib_result");

    return 0;
}
*/


//for 2014 ps beam_test

void ReverseXAxis (TH1 *h)
{
   // Remove the current axis
   h->GetXaxis()->SetLabelOffset(999);
   h->GetXaxis()->SetTickLength(0);

   // Redraw the new axis
   gPad->Update();
   TGaxis *newaxis = new TGaxis(gPad->GetUxmax(),
                                gPad->GetUymin(),
                                gPad->GetUxmin(),
                                gPad->GetUymin(),
                                h->GetXaxis()->GetXmin(),
                                h->GetXaxis()->GetXmax(),
                                510,"-");
   newaxis->SetLabelOffset(-0.03);
   newaxis->Draw();
}

void ReverseYAxis (TH1 *h)
{
   // Remove the current axis
   h->GetYaxis()->SetLabelOffset(999);
   h->GetYaxis()->SetTickLength(0);

   // Redraw the new axis
   gPad->Update();
   TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),
                                gPad->GetUymax(),
                                gPad->GetUxmin()-0.001,
                                gPad->GetUymin(),
                                h->GetYaxis()->GetXmin(),
                                h->GetYaxis()->GetXmax(),
                                510,"+");
   newaxis->SetLabelOffset(-0.03);
   newaxis->Draw();
}

void extract_psd(const char* mipfile,const char* pedfile,const char* outDir,const char* outName)
{
    //gStyle->SetOptStat(0);

    int id8[4][41]={{0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66},
                   {66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,47,46,22,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,1,0},
                   {66,65,64,63,62,61,60,59,58,57,56,55,54,53,52,51,50,49,47,46,22,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,1,0},
                   {0,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,46,47,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66}};
    int id5[4][41]={{23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88},
                   {88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,69,68,45,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,24,23},
                   {88,87,86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,69,68,45,43,42,41,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,24,23},
                   {23,24,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88}};



    TString label[4]={"xpos","xneg","ypos","yneg"};//fee_id:0-xpos,1-xneg,2-ypos,3-yneg
    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];//[fee_id][strip_id][dynode_id]
                                                   //dynode_id: 0-dy5,1-dy8

    Float_t ped_meanbuffer,ped_sigmabuffer;
    TFile* fped=new TFile(pedfile);
    TTree* tped=0;
    TString tmpstr;
    for(int fee_id=0;fee_id<4;fee_id++){
        tmpstr=label[fee_id]+"_ped";
        tped=(TTree*)fped->Get(tmpstr.Data());
        tped->SetBranchAddress("mean",&ped_meanbuffer);
        tped->SetBranchAddress("sigma",&ped_sigmabuffer);
        tped->BuildIndex("channel");
        for(int i=0;i<41;i++){
            tped->GetEntryWithIndex(id8[fee_id][i]+1);
            ped_mean[fee_id][i][1]=ped_meanbuffer;
            ped_sigma[fee_id][i][1]=ped_sigmabuffer;

            tped->GetEntryWithIndex(id5[fee_id][i]+1);
            ped_mean[fee_id][i][0]=ped_meanbuffer;
            ped_sigma[fee_id][i][0]=ped_sigmabuffer;
        }
        delete tped;
    }
    delete fped;

    TFile *file_out=new TFile(Form("%s/%s",outDir,outName),"recreate");
    TTree *tree_ped=new TTree("ped_data","Pedestal from MIPs distribution");
    tree_ped->Branch("ped_mean",ped_mean,"ped_mean[4][41][2]/F");
    tree_ped->Branch("ped_sigma",ped_sigma,"ped_sigma[4][41][2]/F");
    tree_ped->Fill();
    file_out->mkdir("Pedestal");
    file_out->cd("Pedestal");
    tree_ped->Write();
    ///----------------------------

    Int_t xpos[90],xneg[90],ypos[90],yneg[90];
    TFile* fmip=new TFile(mipfile);
    TTree* tmip=(TTree*)fmip->Get("PSD");
    tmip->SetBranchAddress("xpos",xpos);
    tmip->SetBranchAddress("xneg",xneg);
    tmip->SetBranchAddress("ypos",ypos);
    tmip->SetBranchAddress("yneg",yneg);



    Float_t xmean_geo[41],ymean_geo[41];
    Float_t xmean_alg[41],ymean_alg[41];
    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8
    Float_t nopedmipbuffer[4][41][2];

    file_out->cd();
    TTree* tree_raw=new TTree("raw_data","Raw ADC counts of each valid channel");
    tree_raw->SetAutoSave(50000000);
    tree_raw->Branch("raw_data",mipbuffer,"raw_data[4][41][2]/I");
    TTree* tree_noped=new TTree("noped_data","ADC counts without pedestal of each valid channel");
    tree_noped->Branch("noped_data",nopedmipbuffer,"noped_data[4][41][2]/F");
    tree_noped->SetAutoSave(50000000);
    TTree* tree_calib=new TTree("calib_data","ADC counts of each PSD strip");
    tree_calib->Branch("xmean_geo",xmean_geo,"xmean_geo[41]/F");
    tree_calib->Branch("xmean_alg",xmean_alg,"xmean_alg[41]/F");
    tree_calib->Branch("ymean_geo",ymean_geo,"ymean_geo[41]/F");
    tree_calib->Branch("ymean_alg",ymean_alg,"ymean_alg[41]/F");
    tree_calib->SetAutoSave(50000000);

    Int_t entries=tmip->GetEntries();
    for(int i=0;i<entries;i++){
        tmip->GetEntry(i);

        //x
        for(int j=0;j<41;j++){
            mipbuffer[0][j][0]=xpos[id5[0][j]];
            mipbuffer[0][j][1]=xpos[id8[0][j]];
            nopedmipbuffer[0][j][0]=mipbuffer[0][j][0]-ped_mean[0][j][0];
            nopedmipbuffer[0][j][1]=mipbuffer[0][j][1]-ped_mean[0][j][1];

            mipbuffer[1][j][0]=xneg[id5[1][j]];
            mipbuffer[1][j][1]=xneg[id8[1][j]];
            nopedmipbuffer[1][j][0]=mipbuffer[1][j][0]-ped_mean[1][j][0];
            nopedmipbuffer[1][j][1]=mipbuffer[1][j][1]-ped_mean[1][j][1];
            //
            if(nopedmipbuffer[0][j][1]<0) nopedmipbuffer[0][j][1]=0;
            if(nopedmipbuffer[1][j][1]<0) nopedmipbuffer[1][j][1]=0;
            xmean_alg[j]=(nopedmipbuffer[0][j][1]+nopedmipbuffer[1][j][1])/2;
            xmean_geo[j]=TMath::Sqrt(nopedmipbuffer[0][j][1]*nopedmipbuffer[1][j][1]);

        }

        //y
        for(int j=0;j<41;j++){
            mipbuffer[2][j][0]=ypos[id5[2][j]];
            mipbuffer[2][j][1]=ypos[id8[2][j]];
            nopedmipbuffer[2][j][0]=mipbuffer[2][j][0]-ped_mean[2][j][0];
            nopedmipbuffer[2][j][1]=mipbuffer[2][j][1]-ped_mean[2][j][1];

            mipbuffer[3][j][0]=yneg[id5[3][j]];
            mipbuffer[3][j][1]=yneg[id8[3][j]];
            nopedmipbuffer[3][j][0]=mipbuffer[3][j][0]-ped_mean[3][j][0];
            nopedmipbuffer[3][j][1]=mipbuffer[3][j][1]-ped_mean[3][j][1];
            //
            if(nopedmipbuffer[2][j][1]<0) nopedmipbuffer[2][j][1]=0;
            if(nopedmipbuffer[3][j][1]<0) nopedmipbuffer[3][j][1]=0;
            ymean_alg[j]=(nopedmipbuffer[2][j][1]+nopedmipbuffer[3][j][1])/2;
            ymean_geo[j]=TMath::Sqrt(nopedmipbuffer[2][j][1]*nopedmipbuffer[3][j][1]);

        }

        tree_raw->Fill();
        tree_noped->Fill();
        tree_calib->Fill();
    }

    file_out->mkdir("Raw");
    file_out->cd("Raw");
    tree_raw->Write();
    tree_noped->Write();
    file_out->mkdir("Analysis");
    file_out->cd("Analysis");
    tree_calib->Write();


    delete file_out;
    delete fmip;

}

void draw_hitnum(const char* psdfile,const char* basename,const Int_t limit=5)
{
    TFile* file_in=new TFile(psdfile,"update");
    TString outDir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //--------
    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];//[fee_id][strip_id][dynode_id]
                                                   //dynode_id: 0-dy5,1-dy8
    Float_t thresh[4][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();
    delete tree_ped;
    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=ped_mean[0][ch_id][1]+ped_sigma[0][ch_id][1]*limit;
        thresh[1][ch_id]=ped_mean[1][ch_id][1]+ped_sigma[1][ch_id][1]*limit;
        thresh[2][ch_id]=ped_mean[2][ch_id][1]+ped_sigma[2][ch_id][1]*limit;
        thresh[3][ch_id]=ped_mean[3][ch_id][1]+ped_sigma[3][ch_id][1]*limit;
    }

    //-------------
    Float_t xstrip_pos[41]={-400,-380,-360,-340,-320,-300,-280,-260,-240,-220,-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,
                            0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400};
    Float_t ystrip_pos[41]={-400,-380,-360,-340,-320,-300,-280,-260,-240,-220,-200,-180,-160,-140,-120,-100,-80,-60,-40,-20,
                            0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400};

    file_in->cd();
    TH1F *hhit_x=new TH1F(Form("hhit_x_%dsigma",limit),Form("Multi Hit Num Distribution in X plane(%dsigma)",limit),12,-0.5,11.5);
    hhit_x->GetXaxis()->SetTitle("Hit_Number in X plane");
    hhit_x->GetXaxis()->CenterTitle();
    hhit_x->GetYaxis()->SetTitle("Event Num");
    hhit_x->GetYaxis()->CenterTitle();
    TH1F *hhit_y=new TH1F(Form("hhit_y_%dsigma",limit),Form("Multi Hit Num Distribution in Y plane(%dsigma)",limit),12,-0.5,11.5);
    hhit_y->GetXaxis()->SetTitle("Hit_Number in Y plane");
    hhit_y->GetXaxis()->CenterTitle();
    hhit_y->GetYaxis()->SetTitle("Event Num");
    hhit_y->GetYaxis()->CenterTitle();

    TH2F *hpos_xy=new TH2F(Form("hpos_xy_%dsigma",limit),Form("Hit Position Distribution in PSD(%dsigma)",limit),41,-410,410,41,-410,410);
    hpos_xy->GetXaxis()->SetTitle("Y(mm)");
    hpos_xy->GetXaxis()->CenterTitle();
    hpos_xy->GetYaxis()->SetTitle("X(mm)");
    hpos_xy->GetYaxis()->CenterTitle();

    TH2F *hhit_xy=new TH2F(Form("hhit_xy_%dsigma",limit),Form("MultiHit_X VS MultiHit_Y(%dsigma)",limit),12,-0.5,11.5,12,-0.5,11.5);
    hhit_xy->GetXaxis()->SetTitle("X Hit_num");
    hhit_xy->GetXaxis()->CenterTitle();
    hhit_xy->GetYaxis()->SetTitle("Y Hit_num");
    hhit_xy->GetYaxis()->CenterTitle();

    Int_t hitnum_x,hitnum_y,zerohit_x,zerohit_y,zerohit_xy;
    Int_t hitid_x[41],hitid_y[41];
    //Float_t strippos_x[41],strippos_y[41];
    Int_t hit_x,hit_y,hit_xy;
    Int_t hitstrip_x,hitstrip_y;
    Float_t hitpos_x,hitpos_y;

    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8
    //Float_t nopedmipbuffer[4][41][2];
    Float_t xmean_geo[41],ymean_geo[41];

    TTree* tree_raw=(TTree*)file_in->Get("Raw/raw_data");
    tree_raw->SetBranchAddress("raw_data",mipbuffer);
    //TTree* tree_noped=(TTree*)file_in->Get("Raw/noped_data");
    //tree_noped->SetBranchAddress("noped_data",nopedmipbuffer);
    TTree* tree_calib=(TTree*)file_in->Get("Analysis/calib_data");
    tree_calib->SetBranchAddress("xmean_geo",xmean_geo);
    tree_calib->SetBranchAddress("ymean_geo",ymean_geo);
    TTree* tree_hitnum=new TTree(Form("hitnum_%dsigma",limit),Form("hitnum_%dsigma",limit));
    tree_hitnum->Branch("hitnum_x",&hitnum_x,"hitnum_x/I");
    tree_hitnum->Branch("hitnum_y",&hitnum_y,"hitnum_y/I");
    tree_hitnum->Branch("hitid_x",hitid_x,"hitid_x[hitnum_x]/I");
    tree_hitnum->Branch("hitid_y",hitid_y,"hitid_y[hitnum_y]/I");

    zerohit_x=0;zerohit_y=0;zerohit_xy=0;
    hit_x=0;hit_y=0;hit_xy=0;
    Int_t entries=tree_raw->GetEntries();
    for(int i=0;i<entries;i++){
        tree_raw->GetEntry(i);
        tree_calib->GetEntry(i);
        //tree_noped->GetEntry(i);
        hitnum_x=0;hitnum_y=0;
        hitstrip_x=-1;hitstrip_y=-1;

        //x
        for(int j=0;j<41;j++){
            if(mipbuffer[0][j][1]>thresh[0][j] || mipbuffer[1][j][1]>thresh[1][j]){
                hitid_x[hitnum_x]=j;
                hitnum_x++;
            }
        }
        hhit_x->Fill(hitnum_x);

        hitstrip_x=TMath::LocMax(41,xmean_geo);
        hitpos_x=xstrip_pos[hitstrip_x];

        //y
        for(int j=0;j<41;j++){
            if(mipbuffer[2][j][1]>thresh[2][j] || mipbuffer[3][j][1]>thresh[3][j]){
                hitid_y[hitnum_y]=j;
                hitnum_y++;
            }
        }
        hhit_y->Fill(hitnum_y);
        //
        hitstrip_y=TMath::LocMax(41,ymean_geo);
        hitpos_y=ystrip_pos[hitstrip_y];

        ///////////////
        if(hitnum_x==0){
            zerohit_x++;
        }
        if(hitnum_y==0){
            zerohit_y++;
        }
        if(hitnum_x==0 && hitnum_y==0){
            zerohit_xy++;
        }

        if(hitnum_x>0){
            hit_x++;
        }
        if(hitnum_y>0){
            hit_y++;
        }
        if(hitnum_x>0 && hitnum_y>0){
            hit_xy++;
        }
        //pos
        if(hitnum_x>0 && hitnum_y>0){
            hpos_xy->Fill(-hitpos_x,-hitpos_y);
        }
        hhit_xy->Fill(hitnum_x,hitnum_y);
        tree_hitnum->Fill();

    }

    TCanvas* can=new TCanvas("can","can",800,800);
    can->Divide(2,2);
    can->cd(1);
    gPad->SetLogy();
    hhit_x->Draw();
    can->cd(2);
    gPad->SetLogy();
    hhit_y->Draw();
    can->cd(3);
    hpos_xy->Draw();
    gPad->Update();
    ReverseXAxis(hpos_xy);
    ReverseYAxis(hpos_xy);
    can->cd(4);
    hhit_xy->Draw("text");
    can->Print(Form("%s/%s_%dsigma.eps",outDir.Data(),basename,limit));

    FILE* fp=fopen(Form("%s/%s_%dsigma.txt",outDir.Data(),basename,limit),"w");
    fprintf(fp,"zero hit x: %d\n",zerohit_x);
    fprintf(fp,"zero hit y: %d\n",zerohit_y);
    fprintf(fp,"zero hit xy: %d\n",zerohit_xy);
    fprintf(fp,"================\n");
    fprintf(fp,"hit x: %d\t%f\n",hit_x,100.*hit_xy/hit_x);
    fprintf(fp,"hit y: %d\t%f\n",hit_y,100.*hit_xy/hit_y);
    fprintf(fp,"hit xy: %d\n",hit_xy);
    fprintf(fp,"================\n");
    fprintf(fp,"total events: %d\n",entries);
    fclose(fp);

    file_in->cd("Analysis");
    tree_hitnum->Write(0,TObject::kOverwrite);
    if(!file_in->cd("Result")){
        file_in->mkdir("Result");
        file_in->cd("Result");
    }
    hhit_x->Write(0,TObject::kOverwrite);
    hhit_y->Write(0,TObject::kOverwrite);
    hpos_xy->Write(0,TObject::kOverwrite);
    hhit_xy->Write(0,TObject::kOverwrite);

    delete can;
    delete file_in;

}


void bt_draw_channel(const char* psdfile,Int_t limit=5)
{
    TString dir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //----------

    TFile* file_in=new TFile(psdfile);
    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];
    Float_t thresh[4][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();

    delete tree_ped;

    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=ped_mean[0][ch_id][1]+ped_sigma[0][ch_id][1]*limit;
        thresh[1][ch_id]=ped_mean[1][ch_id][1]+ped_sigma[1][ch_id][1]*limit;
        thresh[2][ch_id]=ped_mean[2][ch_id][1]+ped_sigma[2][ch_id][1]*limit;
        thresh[3][ch_id]=ped_mean[3][ch_id][1]+ped_sigma[3][ch_id][1]*limit;
    }
    //------------------------------------------
    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8

    TTree* tree_raw=(TTree*)file_in->Get("Raw/raw_data");
    tree_raw->SetBranchAddress("raw_data",mipbuffer);
    TH1F* hdy8_xpos[41];
    TH1F* hdy8_xneg[41];
    TH1F* hdy8_ypos[41];
    TH1F* hdy8_yneg[41];
    for(int i=0;i<41;i++){
        hdy8_xpos[i]=new TH1F(Form("hdy8_xpos_%d",i+1),Form("hdy8_xpos_%d",i+1),400,-0.5,3999.5);
        hdy8_xneg[i]=new TH1F(Form("hdy8_xneg_%d",i+1),Form("hdy8_xneg_%d",i+1),400,-0.5,3999.5);
        hdy8_ypos[i]=new TH1F(Form("hdy8_ypos_%d",i+1),Form("hdy8_ypos_%d",i+1),400,-0.5,3999.5);
        hdy8_yneg[i]=new TH1F(Form("hdy8_yneg_%d",i+1),Form("hdy8_yneg_%d",i+1),400,-0.5,3999.5);
    }


    Int_t entries=tree_raw->GetEntries();
    for(int i=0;i<entries;i++){
        tree_raw->GetEntry(i);

        for(int j=0;j<41;j++){
            hdy8_xpos[j]->Fill(mipbuffer[0][j][1]);
            hdy8_xneg[j]->Fill(mipbuffer[1][j][1]);
            hdy8_ypos[j]->Fill(mipbuffer[2][j][1]);
            hdy8_yneg[j]->Fill(mipbuffer[3][j][1]);
        }

    }

    TCanvas* can=new TCanvas("can","can",800,800);
    can->Divide(2,2);
    can->Print(Form("%s/mips_anlysis_%dsigma.pdf[",dir.Data(),limit));
    TLine *line_xpos=new TLine;
    line_xpos->SetLineColor(kRed);
    line_xpos->SetLineWidth(1);
    TLine *line_xneg=new TLine;
    line_xneg->SetLineColor(kRed);
    line_xneg->SetLineWidth(1);
    TLine *line_ypos=new TLine;
    line_ypos->SetLineColor(kRed);
    line_ypos->SetLineWidth(1);
    TLine *line_yneg=new TLine;
    line_yneg->SetLineColor(kRed);
    line_yneg->SetLineWidth(1);
    for(int i=0;i<41;i++){
        can->cd(1);
        gPad->SetLogy();
        hdy8_xpos[i]->Draw();
        line_xpos->DrawLine(thresh[0][i],0,thresh[0][i],10000);
        can->cd(2);
        gPad->SetLogy();
        hdy8_xneg[i]->Draw();
        line_xneg->DrawLine(thresh[1][i],0,thresh[1][i],10000);
        can->cd(3);
        gPad->SetLogy();
        hdy8_ypos[i]->Draw();
        line_ypos->DrawLine(thresh[2][i],0,thresh[2][i],10000);
        can->cd(4);
        gPad->SetLogy();
        hdy8_yneg[i]->Draw();
        line_yneg->DrawLine(thresh[3][i],0,thresh[3][i],10000);
        can->Print(Form("%s/mips_anlysis_%dsigma.pdf",dir.Data(),limit));
    }

    can->Print(Form("%s/mips_anlysis_%dsigma.pdf]",dir.Data(),limit));
    delete line_xpos;delete line_xneg;
    delete line_ypos;delete line_yneg;
    delete can;
    delete file_in;

}


void bt_draw_meangeo(const char* psdfile)
{
    TString dir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //----------

    TFile* file_in=new TFile(psdfile);

    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];
    Float_t thresh[2][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();

    delete tree_ped;

    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=TMath::Sqrt(ped_sigma[0][ch_id][1]*ped_sigma[1][ch_id][1]);
        thresh[1][ch_id]=TMath::Sqrt(ped_sigma[2][ch_id][1]*ped_sigma[3][ch_id][1]);
    }

    //------------------------------------------
    Float_t xmean_geo[41],ymean_geo[41];
    Int_t hitnum_x,hitnum_y;
    TTree* tree_calib=(TTree*)file_in->Get("Analysis/calib_data");
    tree_calib->SetBranchAddress("xmean_geo",xmean_geo);
    tree_calib->SetBranchAddress("ymean_geo",ymean_geo);
    TTree* tree_hitnum=(TTree*)file_in->Get("Analysis/hitnum_5sigma");
    tree_hitnum->SetBranchAddress("hitnum_x",&hitnum_x);
    tree_hitnum->SetBranchAddress("hitnum_y",&hitnum_y);
    tree_calib->AddFriend(tree_hitnum);

    TH1F* hdy8_xmean[41];
    TH1F* hdy8_ymean[41];
    TH1F* hdy8_xmean_singlehit[41];
    TH1F* hdy8_ymean_singlehit[41];
    for(int i=0;i<41;i++){
        hdy8_xmean[i]=new TH1F(Form("hdy8_xmean_%d",i+1),Form("hdy8_xmean_%d",i+1),200,-0.5,1999.5);
        hdy8_ymean[i]=new TH1F(Form("hdy8_ymean_%d",i+1),Form("hdy8_ymean_%d",i+1),200,-0.5,1999.5);
        hdy8_xmean_singlehit[i]=new TH1F(Form("hdy8_xmean_singlehit_%d",i+1),Form("hdy8_xmean_singlehit_%d",i+1),200,-0.5,1999.5);
        hdy8_ymean_singlehit[i]=new TH1F(Form("hdy8_ymean_singlehit_%d",i+1),Form("hdy8_ymean_singlehit_%d",i+1),200,-0.5,1999.5);
        hdy8_xmean_singlehit[i]->SetLineColor(kBlack);
        hdy8_ymean_singlehit[i]->SetLineColor(kBlack);
    }


    Int_t entries=tree_calib->GetEntries();
    for(int i=0;i<entries;i++){
        tree_calib->GetEntry(i);

        for(int j=0;j<41;j++){
            hdy8_xmean[j]->Fill(xmean_geo[j]);
            hdy8_ymean[j]->Fill(ymean_geo[j]);
            if(hitnum_x==1){
                hdy8_xmean_singlehit[j]->Fill(xmean_geo[j]);
            }
            if(hitnum_y==1){
                hdy8_ymean_singlehit[j]->Fill(ymean_geo[j]);
            }

        }

    }

    TCanvas* can=new TCanvas("can","can",800,400);
    can->Divide(2,1);
    can->Print(Form("%s/mips_geomean.pdf[",dir.Data()));
    TLine *line_x=new TLine;
    line_x->SetLineColor(kRed);
    line_x->SetLineWidth(1);
    TLine *line_y=new TLine;
    line_y->SetLineColor(kRed);
    line_y->SetLineWidth(1);

    TF1* ffit;
    Float_t mpv,fwhm;
    for(int i=0;i<41;i++){
        can->cd(1);
        gPad->SetLogy();
        //ffit=langaus(hdy8_xmean_singlehit[i],mpv,fwhm,0,thresh[0][i]);
        hdy8_xmean[i]->Draw();
        hdy8_xmean_singlehit[i]->Draw("same");
        //ffit->Draw("lsame");
        line_x->DrawLine(5*thresh[0][i],0,thresh[0][i],10000);
        can->cd(2);
        gPad->SetLogy();
        //ffit=langaus(hdy8_ymean_singlehit[i],mpv,fwhm,0,thresh[1][i]);
        hdy8_ymean[i]->Draw();
        hdy8_ymean_singlehit[i]->Draw("same");
        //ffit->Draw("lsame");
        line_y->DrawLine(5*thresh[1][i],0,thresh[1][i],10000);

        can->Print(Form("%s/mips_geomean.pdf",dir.Data()));
    }

    can->Print(Form("%s/mips_geomean.pdf]",dir.Data()));
    delete line_x;delete line_y;
    delete can;
    delete file_in;
}


void bt_draw_singlehit(const char* psdfile,Int_t limit=5)
{
    TString dir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //----------

    TFile* file_in=new TFile(psdfile);

    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];
    Float_t thresh[4][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();

    delete tree_ped;

    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=ped_mean[0][ch_id][1]+ped_sigma[0][ch_id][1]*limit;
        thresh[1][ch_id]=ped_mean[1][ch_id][1]+ped_sigma[1][ch_id][1]*limit;
        thresh[2][ch_id]=ped_mean[2][ch_id][1]+ped_sigma[2][ch_id][1]*limit;
        thresh[3][ch_id]=ped_mean[3][ch_id][1]+ped_sigma[3][ch_id][1]*limit;
    }
    //------------------------------------------
    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8
    Int_t hitnum_x,hitnum_y;

    TTree* tree_raw=(TTree*)file_in->Get("Raw/raw_data");
    tree_raw->SetBranchAddress("raw_data",mipbuffer);
    TTree* tree_hitnum=(TTree*)file_in->Get(Form("Analysis/hitnum_%dsigma",limit));
    tree_hitnum->SetBranchAddress("hitnum_x",&hitnum_x);
    tree_hitnum->SetBranchAddress("hitnum_y",&hitnum_y);
    tree_raw->AddFriend(tree_hitnum);

    TH1F* hdy8_xpos[41];
    TH1F* hdy8_xneg[41];
    TH1F* hdy8_ypos[41];
    TH1F* hdy8_yneg[41];
    for(int i=0;i<41;i++){
        hdy8_xpos[i]=new TH1F(Form("hdy8_xpos_%d",i+1),Form("hdy8_xpos_%d",i+1),400,-0.5,3999.5);
        hdy8_xneg[i]=new TH1F(Form("hdy8_xneg_%d",i+1),Form("hdy8_xneg_%d",i+1),400,-0.5,3999.5);
        hdy8_ypos[i]=new TH1F(Form("hdy8_ypos_%d",i+1),Form("hdy8_ypos_%d",i+1),400,-0.5,3999.5);
        hdy8_yneg[i]=new TH1F(Form("hdy8_yneg_%d",i+1),Form("hdy8_yneg_%d",i+1),400,-0.5,3999.5);
    }


    Int_t entries=tree_raw->GetEntries();
    for(int i=0;i<entries;i++){
        tree_raw->GetEntry(i);

        if(hitnum_x==1){
            for(int j=0;j<41;j++){
                hdy8_xpos[j]->Fill(mipbuffer[0][j][1]);
                hdy8_xneg[j]->Fill(mipbuffer[1][j][1]);

            }
        }

        if(hitnum_y==1){
            for(int j=0;j<41;j++){
                hdy8_ypos[j]->Fill(mipbuffer[2][j][1]);
                hdy8_yneg[j]->Fill(mipbuffer[3][j][1]);

            }
        }


    }

    TCanvas* can=new TCanvas("can","can",800,800);
    can->Divide(2,2);
    can->Print(Form("%s/mips_onehit_%dsigma.pdf[",dir.Data(),limit));
    TLine *line_xpos=new TLine;
    line_xpos->SetLineColor(kRed);
    line_xpos->SetLineWidth(1);
    TLine *line_xneg=new TLine;
    line_xneg->SetLineColor(kRed);
    line_xneg->SetLineWidth(1);
    TLine *line_ypos=new TLine;
    line_ypos->SetLineColor(kRed);
    line_ypos->SetLineWidth(1);
    TLine *line_yneg=new TLine;
    line_yneg->SetLineColor(kRed);
    line_yneg->SetLineWidth(1);
    for(int i=0;i<41;i++){
        can->cd(1);
        gPad->SetLogy();
        hdy8_xpos[i]->Draw();
        line_xpos->DrawLine(thresh[0][i],0,thresh[0][i],10000);
        can->cd(2);
        gPad->SetLogy();
        hdy8_xneg[i]->Draw();
        line_xneg->DrawLine(thresh[1][i],0,thresh[1][i],10000);
        can->cd(3);
        gPad->SetLogy();
        hdy8_ypos[i]->Draw();
        line_ypos->DrawLine(thresh[2][i],0,thresh[2][i],10000);
        can->cd(4);
        gPad->SetLogy();
        hdy8_yneg[i]->Draw();
        line_yneg->DrawLine(thresh[3][i],0,thresh[3][i],10000);
        can->Print(Form("%s/mips_onehit_%dsigma.pdf",dir.Data(),limit));
    }

    can->Print(Form("%s/mips_onehit_%dsigma.pdf]",dir.Data(),limit));
    delete line_xpos;delete line_xneg;
    delete line_ypos;delete line_yneg;
    delete can;
    delete file_in;

}

void bt_efficiency(const char* psdfile,Int_t sample_strip_x=25,Int_t sample_strip_y=25,Int_t limit=5)
{
    TString dir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //----------

    TFile* file_in=new TFile(psdfile);

    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];
    Float_t thresh[4][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();

    delete tree_ped;

    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=ped_mean[0][ch_id][1]+ped_sigma[0][ch_id][1]*limit;
        thresh[1][ch_id]=ped_mean[1][ch_id][1]+ped_sigma[1][ch_id][1]*limit;
        thresh[2][ch_id]=ped_mean[2][ch_id][1]+ped_sigma[2][ch_id][1]*limit;
        thresh[3][ch_id]=ped_mean[3][ch_id][1]+ped_sigma[3][ch_id][1]*limit;
    }
    //------------------------------------------
    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8
    Int_t hitnum_x,hitnum_y;
    Int_t hitid_x[41],hitid_y[41];

    TTree* tree_raw=(TTree*)file_in->Get("Raw/raw_data");
    tree_raw->SetBranchAddress("raw_data",mipbuffer);
    TTree* tree_hitnum=(TTree*)file_in->Get(Form("Analysis/hitnum_5sigma"));
    tree_hitnum->SetBranchAddress("hitnum_x",&hitnum_x);
    tree_hitnum->SetBranchAddress("hitnum_y",&hitnum_y);
    tree_hitnum->SetBranchAddress("hitid_x",hitid_x);
    tree_hitnum->SetBranchAddress("hitid_y",hitid_y);
    tree_raw->AddFriend(tree_hitnum);


    //Int_t totalnum_x=0;Int_t totalnum_y=0;
    Int_t samplenum_x=0;Int_t samplenum_y=0;
    Int_t firednum_x=0;Int_t firednum_y=0;
    Int_t temp_x,temp_y;
    Int_t entries=tree_raw->GetEntries();
    //totalnum_x=entries;totalnum_y=entries;
    for(int i=0;i<entries;i++){
        tree_raw->GetEntry(i);
        temp_x=0;temp_y=0;
        /*
        if(hitnum_x==1){
            if(hitid_x[0]==sample_strip_x){
                samplenum_x++;
                for(int j=0;j<41;j++){
                    if(mipbuffer[2][j][1]>thresh[2][j] || mipbuffer[3][j][1]>thresh[3][j]){
                        temp_y++;
                    }

                }
                if(temp_y>0)    firednum_y++;
            }
        }

        if(hitnum_y==1){
            if(hitid_y[0]==sample_strip_y){
                samplenum_y++;
                for(int j=0;j<41;j++){
                    if(mipbuffer[0][j][1]>thresh[0][j] || mipbuffer[1][j][1]>thresh[1][j]){
                        temp_x++;
                    }

                }
                if(temp_x>0)    firednum_x++;
            }
        }
        */
        if(hitnum_x==1&&hitnum_y==1){
            if(hitid_x[0]==sample_strip_x){
                samplenum_x++;
                for(int j=0;j<41;j++){
                    if(mipbuffer[2][j][1]>thresh[2][j] || mipbuffer[3][j][1]>thresh[3][j]){
                        temp_y++;
                    }

                }
                if(temp_y>0)    firednum_y++;
            }

            if(hitid_y[0]==sample_strip_y){
                samplenum_y++;
                for(int j=0;j<41;j++){
                    if(mipbuffer[0][j][1]>thresh[0][j] || mipbuffer[1][j][1]>thresh[1][j]){
                        temp_x++;
                    }

                }
                if(temp_x>0)    firednum_x++;
            }
        }

    }

    FILE* fp=fopen(Form("%s/eff_%dsigma.csv",dir.Data(),limit),"w");
    fprintf(fp,"entries\t%d\n",entries);
    fprintf(fp,"X:\n");
    fprintf(fp,"x_sample:\t%d\ty_fired:\t%d\ty_efficiency\t%f\n",samplenum_x,firednum_y,100.*firednum_y/samplenum_x);
    fprintf(fp,"Y:\n");
    fprintf(fp,"y_sample:\t%d\tx_fired:\t%d\tx_efficiency\t%f\n",samplenum_y,firednum_x,100.*firednum_x/samplenum_y);
    fclose(fp);

    delete file_in;
}

void bt_draw_backsca(const char* psdfile)
{
    gStyle->SetOptStat(0);
    TFile* file_in=new TFile(psdfile,"update");
    TString outDir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //
    Int_t hitnum_x,hitnum_y;
    Int_t hitid_x[41],hitid_y[41];

    TTree* tree_hitnum=(TTree*)file_in->Get("Analysis/hitnum_5sigma");
    TH1F* hall_x=new TH1F("hall_x","Distribution of Fired Strips for All evnets(X)",41,-0.5,40.5);
    hall_x->GetXaxis()->SetTitle("Strip_ID");
    hall_x->GetYaxis()->SetTitle("Event_NUM");
    TH1F* hsingle_x=new TH1F("hsingle_x","Distribution of Fired Strips for Single_hit evnets(X)",41,-0.5,40.5);
    TH1F* hdiff_x;
    hall_x->SetLineColor(kBlue);hall_x->SetLineWidth(2);
    hsingle_x->SetLineColor(kRed);hsingle_x->SetLineWidth(2);


    TH1F* hall_y=new TH1F("hall_y","Distribution of Fired Strips for All evnets(Y)",41,-0.5,40.5);
    hall_y->GetXaxis()->SetTitle("Strip_ID");
    hall_y->GetYaxis()->SetTitle("Event_NUM");
    TH1F* hsingle_y=new TH1F("hsingle_y","Distribution of Fired Strips for Single_hit evnets(Y)",41,-0.5,40.5);
    TH1F* hdiff_y;
    hall_y->SetLineColor(kBlue);hall_y->SetLineWidth(2);
    hsingle_y->SetLineColor(kRed);hsingle_y->SetLineWidth(2);


    Double_t ymin = 0;
    Double_t ymax = 0.5;
    Double_t dy = (ymax-ymin)/0.8; //10 per cent margins top and bottom
    Double_t xmin = -0.5;
    Double_t xmax = 40.5;
    Double_t dx = (xmax-xmin)/0.8;

    TCanvas* can=new TCanvas("can","Distribution of Fired Strips",800,400);
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogy();

    tree_hitnum->Draw("hitid_x>>hall_x");
    tree_hitnum->Draw("hitid_x[0]>>hsingle_x","hitnum_x==1","same");
    //gPad->BuildLegend();
    hdiff_x=(TH1F*)hsingle_x->Clone("hdiff_x");
    hdiff_x->SetLineColor(kGreen);hdiff_x->SetLineWidth(2);hdiff_x->SetLineStyle(5);hdiff_x->SetFillStyle(3005);hdiff_x->SetFillColor(kGreen);
    hdiff_x->Divide(hall_x);
    TPad* pad1=new TPad("pad1","",0,0,1,1);
    pad1->SetFillStyle(4000);
    pad1->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad1->Draw();
    pad1->cd();
    hdiff_x->Draw("][sames");
    TGaxis *axis1 = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
    axis1->SetLabelColor(kGreen);
    axis1->Draw();
    TLegend *leg_x=new TLegend(0.2,0.6,0.5,0.9);
    //leg_x->SetNColumns(3);
    leg_x->AddEntry(hall_x,"All","l");
    leg_x->AddEntry(hsingle_x,"Single","l");
    leg_x->AddEntry(hdiff_x,"Single/All","l");
    leg_x->Draw();

    can->cd(2);
    gPad->SetLogy();
    tree_hitnum->Draw("hitid_y>>hall_y");
    tree_hitnum->Draw("hitid_y[0]>>hsingle_y","hitnum_y==1","same");
    //gPad->BuildLegend();
    hdiff_y=(TH1F*)hsingle_y->Clone("hdiff_y");
    hdiff_y->SetLineColor(kGreen);hdiff_y->SetLineWidth(2);hdiff_y->SetLineStyle(5);hdiff_y->SetFillStyle(3005);hdiff_y->SetFillColor(kGreen);
    hdiff_y->Divide(hall_y);
    TPad* pad2=new TPad("pad2","",0,0,1,1);
    pad2->SetFillStyle(4000);
    pad2->Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy);
    pad2->Draw();
    pad2->cd();
    hdiff_y->Draw("][sames");
    TGaxis *axis2 = new TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L");
    axis2->SetLabelColor(kGreen);
    axis2->Draw();
    TLegend *leg_y=new TLegend(0.2,0.6,0.5,0.9);
    //leg_x->SetNColumns(3);
    leg_y->AddEntry(hall_y,"All","l");
    leg_y->AddEntry(hsingle_y,"Single","l");
    leg_y->AddEntry(hdiff_y,"Single/All","l");
    leg_y->Draw();

    can->Print(Form("%s/back_sca.eps",outDir.Data()));
    delete can;
    delete file_in;
}

void bt_draw_singlechannel(const char* psdfile,Int_t chid_x=25,Int_t chid_y=25,Int_t limit=5)
{
    gStyle->SetOptStat(0);

    TString dir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //----------

    TFile* file_in=new TFile(psdfile);

    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];
    Float_t thresh[4][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();

    delete tree_ped;

    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=ped_mean[0][ch_id][1]+ped_sigma[0][ch_id][1]*limit;
        thresh[1][ch_id]=ped_mean[1][ch_id][1]+ped_sigma[1][ch_id][1]*limit;
        thresh[2][ch_id]=ped_mean[2][ch_id][1]+ped_sigma[2][ch_id][1]*limit;
        thresh[3][ch_id]=ped_mean[3][ch_id][1]+ped_sigma[3][ch_id][1]*limit;
    }
    //------------------------------------------
    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8
    Int_t hitnum_x,hitnum_y;

    TTree* tree_raw=(TTree*)file_in->Get("Raw/raw_data");
    tree_raw->SetBranchAddress("raw_data",mipbuffer);
    TTree* tree_hitnum=(TTree*)file_in->Get(Form("Analysis/hitnum_%dsigma",limit));
    tree_hitnum->SetBranchAddress("hitnum_x",&hitnum_x);
    tree_hitnum->SetBranchAddress("hitnum_y",&hitnum_y);
    tree_raw->AddFriend(tree_hitnum);

    TH1F* hdy8_xpos_all=new TH1F("hdy8_xpos_all","hdy8_xpos_all",200,-0.5,1999.5);
    TH1F* hdy8_xneg_all=new TH1F("hdy8_xneg_all","hdy8_xneg_all",200,-0.5,1999.5);
    TH1F* hdy8_ypos_all=new TH1F("hdy8_ypos_all","hdy8_ypos_all",200,-0.5,1999.5);
    TH1F* hdy8_yneg_all=new TH1F("hdy8_yneg_all","hdy8_yneg_all",200,-0.5,1999.5);
    hdy8_xpos_all->SetLineColor(kBlue);
    hdy8_xneg_all->SetLineColor(kBlue);
    hdy8_ypos_all->SetLineColor(kBlue);
    hdy8_yneg_all->SetLineColor(kBlue);

    TH1F* hdy8_xpos_single=new TH1F("hdy8_xpos_single","hdy8_xpos_single",200,-0.5,1999.5);
    TH1F* hdy8_xneg_single=new TH1F("hdy8_xneg_single","hdy8_xneg_single",200,-0.5,1999.5);
    TH1F* hdy8_ypos_single=new TH1F("hdy8_ypos_single","hdy8_ypos_single",200,-0.5,1999.5);
    TH1F* hdy8_yneg_single=new TH1F("hdy8_yneg_single","hdy8_yneg_single",200,-0.5,1999.5);
    hdy8_xpos_single->SetLineColor(kRed);
    hdy8_xneg_single->SetLineColor(kRed);
    hdy8_ypos_single->SetLineColor(kRed);
    hdy8_yneg_single->SetLineColor(kRed);

    TH1F* hdy8_xpos_zero=new TH1F("hdy8_xpos_zero","hdy8_xpos_zero",200,-0.5,1999.5);
    TH1F* hdy8_xneg_zero=new TH1F("hdy8_xneg_zero","hdy8_xneg_zero",200,-0.5,1999.5);
    TH1F* hdy8_ypos_zero=new TH1F("hdy8_ypos_zero","hdy8_ypos_zero",200,-0.5,1999.5);
    TH1F* hdy8_yneg_zero=new TH1F("hdy8_yneg_zero","hdy8_yneg_zero",200,-0.5,1999.5);
    hdy8_xpos_zero->SetLineColor(kGreen);
    hdy8_xneg_zero->SetLineColor(kGreen);
    hdy8_ypos_zero->SetLineColor(kGreen);
    hdy8_yneg_zero->SetLineColor(kGreen);

    Int_t entries=tree_raw->GetEntries();
    for(int i=0;i<entries;i++){
        tree_raw->GetEntry(i);

        hdy8_xpos_all->Fill(mipbuffer[0][chid_x-1][1]);
        hdy8_xneg_all->Fill(mipbuffer[1][chid_x-1][1]);
        hdy8_ypos_all->Fill(mipbuffer[2][chid_y-1][1]);
        hdy8_yneg_all->Fill(mipbuffer[3][chid_y-1][1]);
        if(hitnum_x==1){
            hdy8_xpos_single->Fill(mipbuffer[0][chid_x-1][1]);
            hdy8_xneg_single->Fill(mipbuffer[1][chid_x-1][1]);
        }

        if(hitnum_y==1){
            hdy8_ypos_single->Fill(mipbuffer[2][chid_y-1][1]);
            hdy8_yneg_single->Fill(mipbuffer[3][chid_y-1][1]);
        }

        if(hitnum_x==0 && hitnum_y==0){
            hdy8_xpos_zero->Fill(mipbuffer[0][chid_x-1][1]);
            hdy8_xneg_zero->Fill(mipbuffer[1][chid_x-1][1]);
            hdy8_ypos_zero->Fill(mipbuffer[2][chid_y-1][1]);
            hdy8_yneg_zero->Fill(mipbuffer[3][chid_y-1][1]);
        }

    }

    TLine *line_xpos=new TLine;
    line_xpos->SetLineColor(kBlack);
    line_xpos->SetLineWidth(1);
    TLine *line_xneg=new TLine;
    line_xneg->SetLineColor(kBlack);
    line_xneg->SetLineWidth(1);
    TLine *line_ypos=new TLine;
    line_ypos->SetLineColor(kBlack);
    line_ypos->SetLineWidth(1);
    TLine *line_yneg=new TLine;
    line_yneg->SetLineColor(kBlack);
    line_yneg->SetLineWidth(1);

    TCanvas* can_x=new TCanvas("can_x","can_x",800,400);
    can_x->Divide(2,1);
    can_x->cd(1);
    gPad->SetLogy();
    hdy8_xpos_all->Draw();
    hdy8_xpos_single->Draw("same");
    hdy8_xpos_zero->Draw("same");
    line_xpos->DrawLine(thresh[0][chid_x-1],0,thresh[0][chid_x-1],20000);
    TLegend* legend_xpos=new TLegend(0.6,0.5,0.9,0.9);
    legend_xpos->AddEntry(hdy8_xpos_all,"All","l");
    legend_xpos->AddEntry(hdy8_xpos_single,"Single","l");
    legend_xpos->AddEntry(hdy8_xpos_zero,"Zero","l");
    legend_xpos->AddEntry(line_xpos,Form("%d*Sigma",limit),"l");
    legend_xpos->Draw();
    hdy8_xpos_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_xpos_all->SetTitle(Form("X%d_Pos",chid_x));
    //gPad->BuildLegend();
    can_x->cd(2);
    gPad->SetLogy();
    hdy8_xneg_all->Draw();
    hdy8_xneg_single->Draw("same");
    hdy8_xneg_zero->Draw("same");
    line_xneg->DrawLine(thresh[1][chid_x-1],0,thresh[1][chid_x-1],10000);
    TLegend* legend_xneg=new TLegend(0.6,0.5,0.9,0.9);
    legend_xneg->AddEntry(hdy8_xneg_all,"All","l");
    legend_xneg->AddEntry(hdy8_xneg_single,"Single","l");
    legend_xneg->AddEntry(hdy8_xneg_zero,"Zero","l");
    legend_xneg->AddEntry(line_xneg,Form("%d*Sigma",limit),"l");
    legend_xneg->Draw();
    hdy8_xneg_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_xneg_all->SetTitle(Form("X%d_Neg",chid_x));
    //gPad->BuildLegend();
    can_x->Print(Form("%s/mips_compare_%dsigma_x.eps",dir.Data(),limit));

    TCanvas* can_y=new TCanvas("can_y","can_y",800,400);
    can_y->Divide(2,1);
    can_y->cd(1);
    gPad->SetLogy();
    hdy8_ypos_all->Draw();
    hdy8_ypos_single->Draw("same");
    hdy8_ypos_zero->Draw("same");
    line_ypos->DrawLine(thresh[2][chid_y-1],0,thresh[2][chid_y-1],10000);
    TLegend* legend_ypos=new TLegend(0.6,0.5,0.9,0.9);
    legend_ypos->AddEntry(hdy8_ypos_all,"All","l");
    legend_ypos->AddEntry(hdy8_ypos_single,"Single","l");
    legend_ypos->AddEntry(hdy8_ypos_zero,"Zero","l");
    legend_ypos->AddEntry(line_ypos,Form("%d*Sigma",limit),"l");
    legend_ypos->Draw();
    hdy8_ypos_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_ypos_all->SetTitle(Form("Y%d_Pos",chid_x));
    //gPad->BuildLegend();
    can_y->cd(2);
    gPad->SetLogy();
    hdy8_yneg_all->Draw();
    hdy8_yneg_single->Draw("same");
    hdy8_yneg_zero->Draw("same");
    line_yneg->DrawLine(thresh[3][chid_y-1],0,thresh[3][chid_y-1],10000);
    TLegend* legend_yneg=new TLegend(0.6,0.5,0.9,0.9);
    legend_yneg->AddEntry(hdy8_yneg_all,"All","l");
    legend_yneg->AddEntry(hdy8_yneg_single,"Single","l");
    legend_yneg->AddEntry(hdy8_yneg_zero,"Zero","l");
    legend_yneg->AddEntry(line_yneg,Form("%d*Sigma",limit),"l");
    legend_yneg->Draw();
    hdy8_yneg_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_yneg_all->SetTitle(Form("Y%d_Neg",chid_x));
    //gPad->BuildLegend();
    can_y->Print(Form("%s/mips_compare_%dsigma_y.eps",dir.Data(),limit));

    delete line_xpos;delete line_xneg;
    delete line_ypos;delete line_yneg;
    delete can_x;delete can_y;
    delete file_in;

}

void bt_draw_twohits(const char* psdfile,Int_t chid_x=25,Int_t chid_y=25,Int_t limit=5)
{
    gStyle->SetOptStat(0);

    TString dir=gSystem->DirName(gSystem->ExpandPathName(psdfile));
    //----------

    TFile* file_in=new TFile(psdfile);

    Float_t ped_mean[4][41][2],ped_sigma[4][41][2];
    Float_t thresh[4][41];
    TTree* tree_ped=(TTree*)file_in->Get("Pedestal/ped_data");
    tree_ped->SetBranchAddress("ped_mean",ped_mean);
    tree_ped->SetBranchAddress("ped_sigma",ped_sigma);
    tree_ped->GetEntry();

    delete tree_ped;

    for(int ch_id=0;ch_id<41;ch_id++){
        thresh[0][ch_id]=ped_mean[0][ch_id][1]+ped_sigma[0][ch_id][1]*limit;
        thresh[1][ch_id]=ped_mean[1][ch_id][1]+ped_sigma[1][ch_id][1]*limit;
        thresh[2][ch_id]=ped_mean[2][ch_id][1]+ped_sigma[2][ch_id][1]*limit;
        thresh[3][ch_id]=ped_mean[3][ch_id][1]+ped_sigma[3][ch_id][1]*limit;
    }
    //------------------------------------------
    Int_t mipbuffer[4][41][2];//[fee_id][strip_id][dynode_id]
                                //dynode_id: 0-dy5,1-dy8
    Int_t hitnum_x,hitnum_y;
    Int_t hitid_x[41],hitid_y[41];

    TTree* tree_raw=(TTree*)file_in->Get("Raw/raw_data");
    tree_raw->SetBranchAddress("raw_data",mipbuffer);
    TTree* tree_hitnum=(TTree*)file_in->Get(Form("Analysis/hitnum_%dsigma",limit));
    tree_hitnum->SetBranchAddress("hitnum_x",&hitnum_x);
    tree_hitnum->SetBranchAddress("hitnum_y",&hitnum_y);
    tree_hitnum->SetBranchAddress("hitid_x",hitid_x);
    tree_hitnum->SetBranchAddress("hitid_y",hitid_y);
    tree_raw->AddFriend(tree_hitnum);

    TH1F* hdy8_xpos_all=new TH1F("hdy8_xpos_all","hdy8_xpos_all",300,-0.5,2999.5);
    TH1F* hdy8_xneg_all=new TH1F("hdy8_xneg_all","hdy8_xneg_all",300,-0.5,2999.5);
    TH1F* hdy8_ypos_all=new TH1F("hdy8_ypos_all","hdy8_ypos_all",300,-0.5,2999.5);
    TH1F* hdy8_yneg_all=new TH1F("hdy8_yneg_all","hdy8_yneg_all",300,-0.5,2999.5);
    hdy8_xpos_all->SetLineColor(kBlue);
    hdy8_xneg_all->SetLineColor(kBlue);
    hdy8_ypos_all->SetLineColor(kBlue);
    hdy8_yneg_all->SetLineColor(kBlue);

    TH1F* hdy8_xpos_single=new TH1F("hdy8_xpos_single","hdy8_xpos_single",300,-0.5,2999.5);
    TH1F* hdy8_xneg_single=new TH1F("hdy8_xneg_single","hdy8_xneg_single",300,-0.5,2999.5);
    TH1F* hdy8_ypos_single=new TH1F("hdy8_ypos_single","hdy8_ypos_single",300,-0.5,2999.5);
    TH1F* hdy8_yneg_single=new TH1F("hdy8_yneg_single","hdy8_yneg_single",300,-0.5,2999.5);
    hdy8_xpos_single->SetLineColor(kRed);
    hdy8_xneg_single->SetLineColor(kRed);
    hdy8_ypos_single->SetLineColor(kRed);
    hdy8_yneg_single->SetLineColor(kRed);

    TH1F* hdy8_xpos_double=new TH1F("hdy8_xpos_double","hdy8_xpos_double",300,-0.5,2999.5);
    TH1F* hdy8_xneg_double=new TH1F("hdy8_xneg_double","hdy8_xneg_double",300,-0.5,2999.5);
    TH1F* hdy8_ypos_double=new TH1F("hdy8_ypos_double","hdy8_ypos_double",300,-0.5,2999.5);
    TH1F* hdy8_yneg_double=new TH1F("hdy8_yneg_double","hdy8_yneg_double",300,-0.5,2999.5);
    hdy8_xpos_double->SetLineColor(kBlack);
    hdy8_xneg_double->SetLineColor(kBlack);
    hdy8_ypos_double->SetLineColor(kBlack);
    hdy8_yneg_double->SetLineColor(kBlack);

    Int_t entries=tree_raw->GetEntries();
    for(int i=0;i<entries;i++){
        tree_raw->GetEntry(i);

        hdy8_xpos_all->Fill(mipbuffer[0][chid_x-1][1]);
        hdy8_xneg_all->Fill(mipbuffer[1][chid_x-1][1]);
        hdy8_ypos_all->Fill(mipbuffer[2][chid_y-1][1]);
        hdy8_yneg_all->Fill(mipbuffer[3][chid_y-1][1]);
        if(hitnum_x==1){
            hdy8_xpos_single->Fill(mipbuffer[0][chid_x-1][1]);
            hdy8_xneg_single->Fill(mipbuffer[1][chid_x-1][1]);
            //hdy8_xpos_double->Fill(mipbuffer[0][chid_x-1][1]);
            //hdy8_xneg_double->Fill(mipbuffer[1][chid_x-1][1]);
        }

        if(hitnum_y==1){
            hdy8_ypos_single->Fill(mipbuffer[2][chid_y-1][1]);
            hdy8_yneg_single->Fill(mipbuffer[3][chid_y-1][1]);
            //hdy8_ypos_double->Fill(mipbuffer[2][chid_x-1][1]);
            //hdy8_yneg_double->Fill(mipbuffer[3][chid_x-1][1]);
        }

        if(hitnum_x==2 && TMath::Abs(hitid_x[0]-hitid_x[1])==1){
            if(hitid_x[0]==(chid_x-1) || hitid_x[1]==(chid_x-1)){
                hdy8_xpos_double->Fill(mipbuffer[0][chid_x-1][1]);
                hdy8_xneg_double->Fill(mipbuffer[1][chid_x-1][1]);

                //hdy8_xpos_single->Fill(mipbuffer[0][chid_x-1][1]);
                //hdy8_xneg_single->Fill(mipbuffer[1][chid_x-1][1]);
            }
        }

        if(hitnum_y==2 && TMath::Abs(hitid_y[0]-hitid_y[1])==1){
            if(hitid_y[0]==(chid_y-1) || hitid_y[1]==(chid_y-1)){
                hdy8_ypos_double->Fill(mipbuffer[2][chid_y-1][1]);
                hdy8_yneg_double->Fill(mipbuffer[3][chid_y-1][1]);

                //hdy8_ypos_single->Fill(mipbuffer[2][chid_y-1][1]);
                //hdy8_yneg_single->Fill(mipbuffer[3][chid_y-1][1]);
            }
        }

    }

    TLine *line_xpos=new TLine;
    line_xpos->SetLineColor(kGreen);
    line_xpos->SetLineWidth(1);
    TLine *line_xneg=new TLine;
    line_xneg->SetLineColor(kGreen);
    line_xneg->SetLineWidth(1);
    TLine *line_ypos=new TLine;
    line_ypos->SetLineColor(kGreen);
    line_ypos->SetLineWidth(1);
    TLine *line_yneg=new TLine;
    line_yneg->SetLineColor(kGreen);
    line_yneg->SetLineWidth(1);

    TCanvas* can_x=new TCanvas("can_x","can_x",800,400);
    can_x->Divide(2,1);
    can_x->cd(1);
    gPad->SetLogy();
    hdy8_xpos_all->Draw();
    hdy8_xpos_single->Draw("same");
    hdy8_xpos_double->Draw("same");
    //line_xpos->DrawLine(thresh[0][chid_x-1],0,thresh[0][chid_x-1],20000);
    TLegend* legend_xpos=new TLegend(0.6,0.5,0.9,0.9);
    legend_xpos->AddEntry(hdy8_xpos_all,"All","l");
    legend_xpos->AddEntry(hdy8_xpos_single,"Single","l");
    legend_xpos->AddEntry(hdy8_xpos_double,"Double","l");
    legend_xpos->AddEntry(line_xpos,Form("%d*Sigma",limit),"l");
    legend_xpos->Draw();
    hdy8_xpos_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_xpos_all->SetTitle(Form("X%d_Pos",chid_x));
    //gPad->BuildLegend();
    can_x->cd(2);
    gPad->SetLogy();
    hdy8_xneg_all->Draw();
    hdy8_xneg_single->Draw("same");
    hdy8_xneg_double->Draw("same");
    //line_xneg->DrawLine(thresh[1][chid_x-1],0,thresh[1][chid_x-1],10000);
    TLegend* legend_xneg=new TLegend(0.6,0.5,0.9,0.9);
    legend_xneg->AddEntry(hdy8_xneg_all,"All","l");
    legend_xneg->AddEntry(hdy8_xneg_single,"Single","l");
    legend_xneg->AddEntry(hdy8_xneg_double,"Double","l");
    legend_xneg->AddEntry(line_xneg,Form("%d*Sigma",limit),"l");
    legend_xneg->Draw();
    hdy8_xneg_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_xneg_all->SetTitle(Form("X%d_Neg",chid_x));
    //gPad->BuildLegend();
    can_x->Print(Form("%s/mips_twohits_%dsigma_x.eps",dir.Data(),limit));

    TCanvas* can_y=new TCanvas("can_y","can_y",800,400);
    can_y->Divide(2,1);
    can_y->cd(1);
    gPad->SetLogy();
    hdy8_ypos_all->Draw();
    hdy8_ypos_single->Draw("same");
    hdy8_ypos_double->Draw("same");
    //line_ypos->DrawLine(thresh[2][chid_y-1],0,thresh[2][chid_y-1],10000);
    TLegend* legend_ypos=new TLegend(0.6,0.5,0.9,0.9);
    legend_ypos->AddEntry(hdy8_ypos_all,"All","l");
    legend_ypos->AddEntry(hdy8_ypos_single,"Single","l");
    legend_ypos->AddEntry(hdy8_ypos_double,"Double","l");
    legend_ypos->AddEntry(line_ypos,Form("%d*Sigma",limit),"l");
    legend_ypos->Draw();
    hdy8_ypos_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_ypos_all->SetTitle(Form("Y%d_Pos",chid_x));
    //gPad->BuildLegend();
    can_y->cd(2);
    gPad->SetLogy();
    hdy8_yneg_all->Draw();
    hdy8_yneg_single->Draw("same");
    hdy8_yneg_double->Draw("same");
    //line_yneg->DrawLine(thresh[3][chid_y-1],0,thresh[3][chid_y-1],10000);
    TLegend* legend_yneg=new TLegend(0.6,0.5,0.9,0.9);
    legend_yneg->AddEntry(hdy8_yneg_all,"All","l");
    legend_yneg->AddEntry(hdy8_yneg_single,"Single","l");
    legend_yneg->AddEntry(hdy8_yneg_double,"Double","l");
    legend_yneg->AddEntry(line_yneg,Form("%d*Sigma",limit),"l");
    legend_yneg->Draw();
    hdy8_yneg_all->GetXaxis()->SetTitle("ADC counts");
    hdy8_yneg_all->SetTitle(Form("Y%d_Neg",chid_x));
    //gPad->BuildLegend();
    can_y->Print(Form("%s/mips_twohits_%dsigma_y.eps",dir.Data(),limit));

    delete line_xpos;delete line_xneg;
    delete line_ypos;delete line_yneg;
    delete can_x;delete can_y;
    delete file_in;

}
