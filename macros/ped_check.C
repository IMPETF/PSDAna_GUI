#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"

int get_average(const char* filename,Float_t *sigma_avg)
{
    Int_t ch_id[82]={1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,
                     46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89};
    char label[][10]={"xpos","xneg","ypos","yneg"};
    TFile* file_in=new TFile(filename);
    TTree* tree_in=0;
    Float_t sigma,sum;
    Int_t channel;

    for(int  i=0;i<4;i++){
        tree_in=(TTree*)file_in->Get(Form("%s_ped",label[i]));
        tree_in->BuildIndex("channel");

        tree_in->SetBranchAddress("channel",&channel);
        tree_in->SetBranchAddress("sigma",&sigma);
        //tree_in->SetBranchAddress("mean",&mean);
        sum=0;
        for(int j=0;j<82;j++){
            tree_in->GetEntryWithIndex(ch_id[j]);
            sum+=sigma;
        }
        sigma_avg[i]=sum/82.0;

        delete tree_in;
    }

    delete file_in;

    return 0;
}

int ped_check(const char* dirname,const char* outname,const int test_round_start,const int test_round_stop)
{
    gROOT->Reset();

    Float_t sigma_avg[4];
    Int_t round_id;

    TFile *fileout=new TFile(Form("%s/%s",dirname,outname),"RECREATE");
    TTree *tree_out=new TTree("ped_sigma_avg","Pedestal Sigma Average");
    tree_out->Branch("round_id",&round_id,"round_id/I");
    tree_out->Branch("xpos",sigma_avg,"xpos/F");
    tree_out->Branch("xneg",sigma_avg+1,"xneg/F");
    tree_out->Branch("ypos",sigma_avg+2,"ypos/F");
    tree_out->Branch("yneg",sigma_avg+3,"yneg/F");

    char ped_filename[]="hv.root";
    char dir[256];
    char filename[256];
    for(round_id=test_round_start;round_id<=test_round_stop;round_id++){
        sprintf(dir,"%s/test%d",dirname,round_id);
        if(!(gSystem->OpenDirectory(dir))){
            printf("Warning! : dir not exists '%s'\n",dir);
            return -1;
        }

        sprintf(dir,"%s/test%d/analyze",dirname,round_id);
        if(!(gSystem->OpenDirectory(dir))){
            printf("Warning! : dir not exists '%s'\n",dir);
            return -1;
        }

        sprintf(filename,"%s/test%d/analyze/hv/%s",dirname,round_id,ped_filename);
        get_average(filename,sigma_avg);
        tree_out->Fill();
    }

    fileout->cd();
    tree_out->Write();
    delete fileout;

    return 0;
}

int draw_graph(const char* inname,const char *outname)
{
    char label[][10]={"xpos","xneg","ypos","yneg"};
    Float_t sigma_avg[4];
    Int_t round_id;

    TFile *file_in=new TFile(inname,"update");
    TTree* tree_in=(TTree* )file_in->Get("ped_sigma_avg");
    tree_in->SetBranchAddress("round_id",&round_id);
    tree_in->SetBranchAddress("xpos",sigma_avg);
    tree_in->SetBranchAddress("xneg",sigma_avg+1);
    tree_in->SetBranchAddress("ypos",sigma_avg+2);
    tree_in->SetBranchAddress("yneg",sigma_avg+3);

    TGraph* gr[4];
    TMultiGraph* mg=new TMultiGraph("mg","Ped Sigma");
    for(int i=0;i<4;i++){
        gr[i]=new TGraph();
        gr[i]->SetName(label[i]);
        gr[i]->SetMarkerStyle(20+i);
        gr[i]->SetMarkerColor(2+i);
        gr[i]->SetLineColor(2+i);
        mg->Add(gr[i]);
    }

    Int_t entries=tree_in->GetEntries();
    for(int i=0;i<entries;i++){
        tree_in->GetEntry(i);
        for(int j=0;j<4;j++){
            gr[j]->SetPoint(i,round_id,sigma_avg[j]);
        }
    }

    TCanvas* can=new TCanvas("can","can",500,500);
    mg->Draw("APL");
    can->BuildLegend();

    file_in->cd();
    mg->Write(0,TObject::kOverwrite);

    delete file_in;

    return 0;
}
