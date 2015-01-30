#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TGClient.h"
#include "TGWindow.h"
#include "TGSplitter.h"
#include "TGLayout.h"
#include "TVirtualX.h"
#include "TG3DLine.h"
#include "TGText.h"
#include "TGLabel.h"
#include "TGComboBox.h"
#include "TGListView.h"
#include "TGFSContainer.h"
#include "TGButton.h"
#include "TGTextView.h"
#include "TCanvasImp.h"
#include "TRootCanvas.h"
#include "TStyle.h"
#include "TVirtualX.h"
#include "TString.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TKey.h"
#include "TGFSContainer.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TIterator.h"
#include "TCollection.h"
#include "TMarker.h"
#include "TAxis.h"
#include "TContextMenu.h"
#include "TClassMenuItem.h"
#include "TGTextEntry.h"
#include "TGMsgBox.h"
#include "TGProgressBar.h"
#include "TGTextBuffer.h"
#include "TGSplitter.h"
#include "TGNumberEntry.h"
#include "TVirtualPadEditor.h"
#include "TGFrame.h"
#include "TRootEmbeddedCanvas.h"
#include "TGFileDialog.h"
#include "TGButtonGroup.h"
#include "TTree.h"
#include <ctime>
#include <string.h>
#include <iostream>
#include <fstream>
#include "GuiFrame.h"
using namespace std;

ClassImp(GuiFrame)

//------globals----------////
int triggercounttotal1=0; int triggercounttotal2=0; int triggercounttotal3=0; int triggercounttotal4=0;

GuiFrame *gGuiFrame =0;
//-------------------------------
void HBOOK1(const char *hname,const char *title, int nxbin, float xlow, float xup, float vmx)
{
    if(!(gGuiFrame->datafile)) return;
    new TH1F(hname, title, nxbin, xlow, xup);
}

void HBOOK2(const char *hname,const char *title, int nxbin, float xlow, float xup, int nybin, float ylow, float yup, float vmx)
{
    if(!(gGuiFrame->datafile)) return;
    new TH2F(hname, title, nxbin, xlow, xup, nybin, ylow, yup);
}

void HF1(const char *hname, float x, float weight)
{
    if(!(gGuiFrame->datafile)) return;
    if(gGuiFrame->datafile->Get(hname))
    ((TH1F *)(gGuiFrame->datafile->Get(hname)))->Fill(x, weight);

}

void HF1cont(const char *hname, int x, float y)
{
    if(!(gGuiFrame->datafile)) return;
    if(gGuiFrame->datafile->Get(hname))
    {
        if(x==0)
        ((TH1F *)(gGuiFrame->datafile->Get(hname)))->Reset();
        ((TH1F *)(gGuiFrame->datafile->Get(hname)))->GetXaxis()->SetRange(0,x+30);
        ((TH1F *)(gGuiFrame->datafile->Get(hname)))->SetBinContent(x, y);
    }
        //((TH1F *)(gDirectory->Get(hname)))->Fill(x, y);
}

void HF2(const char *hname, float x, float y, float weight)
{
    if(!(gGuiFrame->datafile)) return;
    if(gGuiFrame->datafile->Get(hname))
    ((TH2F *)(gGuiFrame->datafile->Get(hname)))->Fill(x, y, weight);
}

#include "TH_init.icc"
#include "TH_analyze.icc"
#include "test_for_gui.C"

//--------GuiFrame Definition-----------------
GuiFrame::GuiFrame(const TGWindow *p, UInt_t w, UInt_t h):
TGMainFrame(p, w, h, kHorizontalFrame)

{
    if(gGuiFrame!=0)
    {
        cout<<"Only one instance of 'GuiFrame' permitted"<<endl;
        return;
    }
    gGuiFrame = this;
    bin_dir=gSystem->WorkingDirectory();
    /////////////////////////////////
    SetCleanup(kDeepCleanup);
    Connect("CloseWindow()", "GuiFrame", this, "CloseWindow()");
    DontCallClose();
    ObjList = 0;
    objcurr = 0;
    fClose = kTRUE;
    //----------Left Frame---------
    fFFileList = new TGGroupFrame(this,"Engineer Data Decoding",kVerticalFrame|kFixedWidth);
    fFFileList->SetTitlePos(TGGroupFrame::kCenter);
    fFFileList->Resize(350,500);
    AddFrame(fFFileList, new TGLayoutHints(kLHintsExpandY, 1, 1, 1, 1));
    //----------LeftFrame: open dialog--------
        TGHorizontalFrame *butt1 = new TGHorizontalFrame(fFFileList, 380, 25);
            fLFilename_Engineer = new TGLabel(butt1, new TGString("FileName:"));
            fFileEngineer = new TGTextEntry(butt1, fTbEngineer = new TGTextBuffer(200));
            //fFileEngineer->Resize(380, fFileEngineer->GetDefaultHeight());
            fBdisplayfile  =  new TGTextButton(butt1, " &Open..",kB_openEngineer);
            fBdisplayfile->Connect("Clicked()", "GuiFrame", this, "OpenEngineerData()");

            butt1->AddFrame(fLFilename_Engineer, new TGLayoutHints(kLHintsCenterY|kLHintsLeft,1,1,1,1));
            butt1->AddFrame(fFileEngineer, new TGLayoutHints(kLHintsCenterY|kLHintsCenterX|kLHintsExpandX,1,1,1,1));
            butt1->AddFrame(fBdisplayfile, new TGLayoutHints(kLHintsCenterY|kLHintsRight,1,1,1,1));

        fFFileList->AddFrame(butt1,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
        //----------LeftFrame: seperate line 1--
        TGHorizontal3DLine *hsep_1 = new TGHorizontal3DLine(fFFileList, 380, 3);
        fFFileList->AddFrame(hsep_1, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

        //----------LeftFrame: control panel----
        TGLayoutHints *flayoutButt = new TGLayoutHints(kLHintsCenterY|kLHintsCenterX|kLHintsExpandX|kLHintsExpandY,1,1,1,1);
        TGHorizontalFrame *butt2 = new TGHorizontalFrame(fFFileList, 600, 50, kFitWidth);
            fBonefile  =  new TGTextButton(butt2, " &Run ",kB_process);
            fBonefile->Connect("Clicked()", "GuiFrame", this, "ProcessFiles()");
            fBstop  =  new TGTextButton(butt2, " &Stop ",kB_stop);
            fBstop->Connect("Clicked()", "GuiFrame", this, "Stop()");
            fBresetall =  new TGTextButton(butt2, " RESETALL ", kB_resetall);
            fBresetall->Connect("Clicked()", "GuiFrame", this, "ImB_ResetAllTH()");
            fBresetcurr = new TGTextButton(butt2, " RESETCUR ", kB_resetcurr);
            fBresetcurr->Connect("Clicked()", "GuiFrame", this, "ImB_ResetCurrTH()");

            butt2->AddFrame(fBonefile, flayoutButt);
            butt2->AddFrame(fBstop, flayoutButt);
            butt2->AddFrame(fBresetcurr, flayoutButt);
            butt2->AddFrame(fBresetall,  flayoutButt);

        fFFileList->AddFrame(butt2,new TGLayoutHints(kLHintsExpandX,1,1,1,1));

        TGHorizontalFrame *butt3 = new TGHorizontalFrame(fFFileList, 600, 50, kFitWidth);
            fBintegral =  new TGTextButton(butt3, " INTEGRAL ", kB_integral);
            fBintegral->Connect("Clicked()", "GuiFrame", this, "ImB_integral()");
            fBprevious = new TGTextButton(butt3, " PREVIOUS ", kB_privious);
            fBprevious->Connect("Clicked()", "GuiFrame", this, "ImB_previous()");
            fBnext     = new TGTextButton(butt3, "   NEXT   ", kB_next);
            fBnext->Connect("Clicked()", "GuiFrame", this, "ImB_next()");
            fBupdate   = new TGTextButton(butt3, "  UPDATE  ", kB_update);
            fBupdate->Connect("Clicked()", "GuiFrame", this, "ImB_update()");

            butt3->AddFrame(fBintegral,  flayoutButt);
            butt3->AddFrame(fBprevious, flayoutButt);
            butt3->AddFrame(fBnext, flayoutButt);
            butt3->AddFrame(fBupdate, flayoutButt);

            fBonefile->SetToolTipText("Process file");
            fBstop->SetToolTipText("Stop processing");
            fBresetall->SetToolTipText("Reset all histograms");
            fBresetcurr->SetToolTipText("Reset current histogram");
            fBintegral->SetToolTipText("Integral given bins");
            fBprevious->SetToolTipText("Draw previous histogram");
            fBnext->SetToolTipText("Draw next histogram");
            fBupdate->SetToolTipText("Update current histogram");

        fFFileList->AddFrame(butt3,new TGLayoutHints(kLHintsExpandX,1,1,1,1));

        //---LeftFrame: seperate line 2----------
        TGHorizontal3DLine *hsep_2 = new TGHorizontal3DLine(fFFileList, 380, 3);
        fFFileList->AddFrame(hsep_2, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));

        //---LeftFrame: List View--------------
        flvFile = new TGListView(fFFileList, 380, 480);
        fFFileList->AddFrame(flvFile, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY|kLHintsCenterX, 1, 1, 1, 1));
        Pixel_t white;
        gClient->GetColorByName("white", white);
        fFileCont = new TGFileContainer(flvFile, kSunkenFrame, white);
        fFileCont->Connect("DoubleClicked(TGFrame* , Int_t)", "GuiFrame", this, "OnDoubleClick(TGLVEntry* ,Int_t)");
        //fFileCont->Associate(this);

        //---LeftFrame: Progress Bar-----------
        TGLayoutHints *fLHPr = new TGLayoutHints(kLHintsBottom |kLHintsExpandX, 5, 5,  5, 10);
        fHProgH = new TGHProgressBar(fFFileList,TGProgressBar::kStandard, 200);
        fHProgH->ShowPosition();
        fHProgH->SetFillType(TGProgressBar::kBlockFill);
        fFFileList->AddFrame(fHProgH, fLHPr);

    //---Splitter Left------------------
    TGVSplitter *vsp_left = new TGVSplitter(this, 4, 4);
    vsp_left->SetFrame(fFFileList, kTRUE);
    AddFrame(vsp_left, new TGLayoutHints(kLHintsLeft|kLHintsExpandY));

    /********************************************************************/
    //---Center Frame----------------------------
    fFDrawPanel = new TGVerticalFrame(this, 600, 500);
    AddFrame(fFDrawPanel,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
        //----------CenterFrame: seperate line 3--
        TGHorizontal3DLine *hsep_3 = new TGHorizontal3DLine(fFDrawPanel, 380, 3);
        fFDrawPanel->AddFrame(hsep_3, new TGLayoutHints(kLHintsExpandX, 1, 1, 1, 1));
        //---CenterFrame: draw frame------------
        TGVerticalFrame* fFCanvasFrame=new TGVerticalFrame(fFDrawPanel,600,600,kFixedHeight);
            TGHorizontalFrame* fFDrawFrame=new TGHorizontalFrame(fFCanvasFrame,600,25);
                //--draw frame: draw option------
                TGHorizontalFrame *fFdrawopt = new TGHorizontalFrame(fFDrawFrame, 180, 25);
                flabDrawOpt = new TGLabel(fFdrawopt, " Draw Option: ");
                fFdrawopt->AddFrame(flabDrawOpt, new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsCenterY, 1, 1, 1, 1));
                fcomDrawOpt = new TGComboBox(fFdrawopt);
                fcomDrawOpt->Resize(90, 25);
                fcomDrawOpt->AddEntry("",    0);
                fcomDrawOpt->AddEntry("COLZ", 1);
                fcomDrawOpt->AddEntry("LEGO", 2);
                fcomDrawOpt->AddEntry("CONT", 3);
                fcomDrawOpt->AddEntry("SURF", 4);
                fFdrawopt->AddFrame(fcomDrawOpt, new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsCenterY, 1, 1, 1,1));
            fFDrawFrame->AddFrame(fFdrawopt, new TGLayoutHints(kLHintsLeft|kLHintsExpandX, 1, 1, 1, 1));

            //---draw frame: Divide Pad------------
            TGHorizontalFrame *fFDivPad = new TGHorizontalFrame(fFDrawFrame, 180, 25);
                flabDivPad = new TGLabel(fFDivPad, " Divide Pad: ");
                fFDivPad->AddFrame(flabDivPad, new TGLayoutHints(kLHintsLeft|kLHintsExpandX|kLHintsCenterY, 10, 1, 1, 1));
                fcomDivPad = new TGComboBox(fFDivPad);
                fcomDivPad->Resize(90, 25);
                fcomDivPad->AddEntry("1 X 1", 0);
                fcomDivPad->AddEntry("2 X 2", 1);
                fcomDivPad->AddEntry("2 X 3", 2);
                fcomDivPad->AddEntry("3 X 3", 3);
                fcomDivPad->Select(0, false);
                fFDivPad->AddFrame(fcomDivPad, new TGLayoutHints(kLHintsRight|kLHintsExpandX|kLHintsCenterY, 1, 1, 1,1));
            fFDrawFrame->AddFrame(fFDivPad, new TGLayoutHints(kLHintsRight|kLHintsExpandX, 1, 1, 1, 1));
        fFCanvasFrame->AddFrame(fFDrawFrame,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));

        //---CenterFrame: seperate line 4----------
        TGHorizontal3DLine *hsep_4 = new TGHorizontal3DLine(fFCanvasFrame, 380, 3);
        fFCanvasFrame->AddFrame(hsep_4, new TGLayoutHints(kLHintsTop|kLHintsExpandX, 1, 1, 1, 1));

        //---CenterFrame: Embedded Canvas--------------------------
        fEmbeddedCanvas = new TRootEmbeddedCanvas("canpf",fFCanvasFrame,600,600);
        canpf=fEmbeddedCanvas->GetCanvas();
        fFCanvasFrame->AddFrame(fEmbeddedCanvas,new TGLayoutHints(kLHintsTop|kLHintsExpandX|kLHintsExpandY,1,1,1,1));

        fFDrawPanel->AddFrame(fFCanvasFrame,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));
        //---CenterFrame: Split line-----------------
        TGHSplitter* hsp=new TGHSplitter(fFDrawPanel,600,2);
        hsp->SetFrame(fFCanvasFrame,kTRUE);
        fFDrawPanel->AddFrame(hsp,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));
        //---CenterFrame: Text View-------------------------
        TGHorizontalFrame *fFvtext = new TGHorizontalFrame(fFDrawPanel, 600, 100, kRaisedFrame);
            fviewText = new TGTextView(fFvtext, 600, 190);
            fFvtext->AddFrame(fviewText, new TGLayoutHints(kLHintsExpandY|kLHintsExpandX, 2, 2, 6, 2));//
        fFDrawPanel->AddFrame(fFvtext, new TGLayoutHints(kLHintsBottom|kLHintsExpandY|kLHintsExpandX, 2, 2, 6, 2));//
    /****************************************************************/
    //---Splitter Right------------------
    TGVSplitter *vsp_right = new TGVSplitter(this, 4, 4);
    AddFrame(vsp_right, new TGLayoutHints(kLHintsLeft|kLHintsExpandY));
    //---Right Frame--------------------
    fFScidataPanel = new TGGroupFrame(this,"Science Data Decoding",kVerticalFrame|kFixedWidth);
    fFScidataPanel->SetTitlePos(TGGroupFrame::kCenter);
    fFScidataPanel->Resize(350,500);
    AddFrame(fFScidataPanel,new TGLayoutHints(kLHintsRight|kLHintsExpandY,1,1,1,1));
        //---RightFrame: Config --------------
            //--Process--
        fHBG_process = new TGHButtonGroup(fFScidataPanel,"Process");
            new TGRadioButton(fHBG_process,"RawDataDecoding",kB_rawdatadecoding);
            new TGRadioButton(fHBG_process,"RootFileAnalyzing",kB_rootfileanalyze);
        fHBG_process->SetButton(kB_rawdatadecoding);
        fHBG_process->Connect("Pressed(Int_t)","GuiFrame",this,"OnProcessConfig(Int_t)");
        fFScidataPanel->AddFrame(fHBG_process,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
            //--Type--
        fHBG_type = new TGHButtonGroup(fFScidataPanel,"Type");
            new TGRadioButton(fHBG_type,"Pedestal",kB_pedestal);
            new TGRadioButton(fHBG_type,"Calibration",kB_calibration);
            new TGRadioButton(fHBG_type,"MIPs",kB_mips);
        fHBG_type->SetButton(kB_pedestal);
        analyze_type=kB_pedestal;
        rawdata_type=kNormal;
        fHBG_type->Connect("Pressed(Int_t)","GuiFrame",this,"OnTypeConfig(Int_t)");
        fFScidataPanel->AddFrame(fHBG_type,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
            //--Mode--
        fHBG_mode = new TGHButtonGroup(fFScidataPanel,"Mode");
            new TGRadioButton(fHBG_mode,"Normal",kB_normal);
            new TGRadioButton(fHBG_mode,"Compressed",kB_compressed);
            new TGRadioButton(fHBG_mode,"Smaller",kB_smaller);
        fHBG_mode->SetButton(kB_normal);
        fHBG_mode->Connect("Pressed(Int_t)","GuiFrame",this,"OnModeConfig(Int_t)");
        //fHBG_mode->SetState(kFALSE);
        fFScidataPanel->AddFrame(fHBG_mode,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
            //--Time Configutation--
        TGGroupFrame *fGF_time = new TGGroupFrame(fFScidataPanel,"Time Info",kVerticalFrame);
        fGF_time->Resize(350,100);
        fFScidataPanel->AddFrame(fGF_time,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
            TGHorizontalFrame *fFH_ymd=new TGHorizontalFrame(fGF_time,350,25);
                fNE_year=new TGNumberEntry(fFH_ymd,2013,4,kNE_year,TGNumberFormat::kNESInteger,
                                                                TGNumberFormat::kNEANonNegative,
                                                                TGNumberFormat::kNELLimitMinMax,2013,2100);
                TGLabel *label_year=new TGLabel(fFH_ymd,"Year:");
                fFH_ymd->AddFrame(label_year,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fFH_ymd->AddFrame(fNE_year,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                fNE_month=new TGNumberEntry(fFH_ymd,1,2,kNE_month,TGNumberFormat::kNESInteger,
                                                                TGNumberFormat::kNEANonNegative,
                                                                TGNumberFormat::kNELLimitMinMax,1,12);
                TGLabel *label_month=new TGLabel(fFH_ymd,"Month:");
                fFH_ymd->AddFrame(label_month,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fFH_ymd->AddFrame(fNE_month,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                fNE_day=new TGNumberEntry(fFH_ymd,1,2,kNE_day,TGNumberFormat::kNESInteger,
                                                                TGNumberFormat::kNEANonNegative,
                                                                TGNumberFormat::kNELLimitMinMax,1,31);
                TGLabel *label_day=new TGLabel(fFH_ymd,"Day:");
                fFH_ymd->AddFrame(label_day,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fFH_ymd->AddFrame(fNE_day,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
            fGF_time->AddFrame(fFH_ymd,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));

            TGHorizontalFrame *fFH_hms=new TGHorizontalFrame(fGF_time,350,25);
                fNE_hour=new TGNumberEntry(fFH_hms,0,2,kNE_hour,TGNumberFormat::kNESInteger,
                                                                TGNumberFormat::kNEANonNegative,
                                                                TGNumberFormat::kNELLimitMinMax,0,23);
                TGLabel *label_hour=new TGLabel(fFH_hms,"Hour:");
                fFH_hms->AddFrame(label_hour,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fFH_hms->AddFrame(fNE_hour,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                fNE_minute=new TGNumberEntry(fFH_hms,0,2,kNE_minute,TGNumberFormat::kNESInteger,
                                                                TGNumberFormat::kNEANonNegative,
                                                                TGNumberFormat::kNELLimitMinMax,0,59);
                TGLabel *label_minute=new TGLabel(fFH_hms,"Minute:");
                fFH_hms->AddFrame(label_minute,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fFH_hms->AddFrame(fNE_minute,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                fNE_second=new TGNumberEntry(fFH_hms,0,2,kNE_second,TGNumberFormat::kNESInteger,
                                                                TGNumberFormat::kNEANonNegative,
                                                                TGNumberFormat::kNELLimitMinMax,0,59);
                TGLabel *label_second=new TGLabel(fFH_hms,"Second:");
                fFH_hms->AddFrame(label_second,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fFH_hms->AddFrame(fNE_second,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
            fGF_time->AddFrame(fFH_hms,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));

             TGHorizontalFrame *fFH_gettime=new TGHorizontalFrame(fGF_time,350,25);
                fCBN_usetimecode=new TGCheckButton(fFH_gettime,"UseTimecode",kB_useTimecode);
                fCBN_usetimecode->SetToolTipText("Enable/Disable start/stop time configuration");
                fCBN_usetimecode->Connect("Toggled(Bool_t)","GuiFrame",this,"OnUseTimecode(Bool_t)");
                fCBN_usetimecode->SetOn();
                fFH_gettime->AddFrame(fCBN_usetimecode,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                fTBN_stoptime = new TGTextButton(fFH_gettime,"STOP",kB_getStoptime);
                fTBN_stoptime->Connect("Clicked()","GuiFrame",this,"OnGetStopTime()");
                fTBN_stoptime->SetToolTipText("Configure Stop Time");
                fFH_gettime->AddFrame(fTBN_stoptime,new TGLayoutHints(kLHintsRight|kLHintsExpandX,1,1,1,1));
                fTBN_starttime = new TGTextButton(fFH_gettime,"START",kB_getStarttime);
                fTBN_starttime->SetToolTipText("Configure Start Time");
                fTBN_starttime->Connect("Clicked()","GuiFrame",this,"OnGetStartTime()");
                fFH_gettime->AddFrame(fTBN_starttime,new TGLayoutHints(kLHintsRight|kLHintsExpandX,1,1,1,1));
             fGF_time->AddFrame(fFH_gettime,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));
        //---RightFrame: Analyze-------------
        fGF_analyze = new TGGroupFrame(fFScidataPanel,"Decode",kVerticalFrame);
        fGF_analyze->SetTitlePos(TGGroupFrame::kCenter);
        fGF_analyze->Resize(350,300);
        fFScidataPanel->AddFrame(fGF_analyze,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
        //--input--
        TGHorizontalFrame *input = new TGHorizontalFrame(fGF_analyze, 350, 25);
            fFile_Sci_input = new TGTextEntry(input, fTb_Sci_input = new TGTextBuffer(200));
            fBopen_Sci_input  =  new TGTextButton(input, " &Input..",kB_inputSciData);
            fBopen_Sci_input->Connect("Clicked()", "GuiFrame", this, "OnInputSciData()");
            input->AddFrame(fFile_Sci_input, new TGLayoutHints(kLHintsCenterY|kLHintsCenterX|kLHintsExpandX,1,1,1,1));
            input->AddFrame(fBopen_Sci_input, new TGLayoutHints(kLHintsCenterY|kLHintsRight,1,1,1,1));
        fGF_analyze->AddFrame(input,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
        //--output--
        TGHorizontalFrame *output = new TGHorizontalFrame(fGF_analyze, 350, 25);
            fFile_Sci_output = new TGTextEntry(output, fTb_Sci_output = new TGTextBuffer(200));
            fBopen_Sci_output  =  new TGTextButton(output, " &Output..",kB_outputSciData);
            fBopen_Sci_output->Connect("Clicked()", "GuiFrame", this, "OnOutputSciData()");
            output->AddFrame(fFile_Sci_output, new TGLayoutHints(kLHintsCenterY|kLHintsCenterX|kLHintsExpandX,1,1,1,1));
            output->AddFrame(fBopen_Sci_output, new TGLayoutHints(kLHintsCenterY|kLHintsRight,1,1,1,1));
        fGF_analyze->AddFrame(output,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
        //--ref---
        TGHorizontalFrame *ref = new TGHorizontalFrame(fGF_analyze, 350, 25);
            fFile_Sci_ref = new TGTextEntry(ref, fTb_Sci_ref = new TGTextBuffer(200));
            fBopen_Sci_ref = new TGTextButton(ref, " &Reference..",kB_refSciData);
            fBopen_Sci_ref->Connect("Clicked()", "GuiFrame", this, "OnRefSciData()");
            fFile_Sci_ref->SetState(kFALSE);
            fBopen_Sci_ref->SetState(kButtonDisabled);
            //fBopen_Sci_ref->SetState(kButtonUp);
            ref->AddFrame(fFile_Sci_ref, new TGLayoutHints(kLHintsCenterY|kLHintsCenterX|kLHintsExpandX,1,1,1,1));
            ref->AddFrame(fBopen_Sci_ref, new TGLayoutHints(kLHintsCenterY|kLHintsRight,1,1,1,1));
        fGF_analyze->AddFrame(ref,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
        //--command--
        TGGroupFrame *fGF_command=new TGGroupFrame(fGF_analyze,"Command",kHorizontalFrame);
            fTBN_analyzeSci = new TGTextButton(fGF_command,"Decode",kB_analyzeSci);
            fTBN_analyzeSci->Connect("Clicked()","GuiFrame",this,"OnSciDecode()");
            fTBN_showResult = new TGTextButton(fGF_command,"Analyze",kB_showResult);
            fTBN_showResult->Connect("Clicked()","GuiFrame",this,"OnSciAnalyze()");
            fTBN_showResult->SetState(kButtonDisabled);
            fGF_command->AddFrame(fTBN_analyzeSci,new TGLayoutHints(kLHintsExpandX,2,2,10,1));
            fGF_command->AddFrame(fTBN_showResult,new TGLayoutHints(kLHintsExpandX,2,2,10,1));
        fGF_analyze->AddFrame(fGF_command,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
        //--result--
        TGGroupFrame *fGF_result=new TGGroupFrame(fGF_analyze,"Result",kVerticalFrame);
        fGF_result->Resize(350,300);
            //--result: channels---
            TGHorizontalFrame *fHF_channels=new TGHorizontalFrame(fGF_result,350,100);
                TGVerticalFrame *fVF_channels=new TGVerticalFrame(fHF_channels,300,100);
                    TGHorizontalFrame *fHF_channels_par1=new TGHorizontalFrame(fVF_channels,300,50);
                        TGLabel* fLResult_layer=new TGLabel(fHF_channels_par1,"Layer:");
                        fHF_channels_par1->AddFrame(fLResult_layer,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                        fCB_layer=new TGComboBox(fHF_channels_par1);
                        fCB_layer->AddEntry("X",0);
                        fCB_layer->AddEntry("Y",1);
                        fCB_layer->Select(0,kFALSE);
                        fHF_channels_par1->AddFrame(fCB_layer,new TGLayoutHints(kLHintsExpandX|kLHintsLeft|kLHintsExpandY,1,1,1,1));

                        TGLabel* fLResult_strip=new TGLabel(fHF_channels_par1,"Strip:");
                        fHF_channels_par1->AddFrame(fLResult_strip,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                        fCB_strip=new TGComboBox(fHF_channels_par1);
                        for(int i=0;i<41;i++){
                            fCB_strip->AddEntry(Form("%d",i+1),i);
                        }
                        fCB_strip->Select(0,kFALSE);
                        fHF_channels_par1->AddFrame(fCB_strip,new TGLayoutHints(kLHintsExpandX|kLHintsLeft|kLHintsExpandY,1,1,1,1));
                    fVF_channels->AddFrame(fHF_channels_par1,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));

                    TGHorizontalFrame *fHF_channels_par2=new TGHorizontalFrame(fVF_channels,300,50);
                        TGLabel* fLResult_side=new TGLabel(fHF_channels_par2,"Side:");
                        fHF_channels_par2->AddFrame(fLResult_side,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                        fCB_side=new TGComboBox(fHF_channels_par2);
                        fCB_side->AddEntry("Positive",0);
                        fCB_side->AddEntry("Negative",1);
                        fCB_side->Select(0,kFALSE);
                        fHF_channels_par2->AddFrame(fCB_side,new TGLayoutHints(kLHintsExpandX|kLHintsLeft|kLHintsExpandY,1,1,1,1));

                        TGLabel* fLResult_dynode=new TGLabel(fHF_channels_par2,"Dynode:");
                        fHF_channels_par2->AddFrame(fLResult_dynode,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                        fCB_dynode=new TGComboBox(fHF_channels_par2);
                        fCB_dynode->AddEntry("Dy8",0);
                        fCB_dynode->AddEntry("Dy5",1);
                        fCB_dynode->Select(0,kFALSE);
                        fHF_channels_par2->AddFrame(fCB_dynode,new TGLayoutHints(kLHintsExpandX|kLHintsLeft|kLHintsExpandY,1,1,1,1));
                    fVF_channels->AddFrame(fHF_channels_par2,new TGLayoutHints(kLHintsTop|kLHintsExpandX,1,1,1,1));
                fHF_channels->AddFrame(fVF_channels,new TGLayoutHints(kLHintsExpandX|kLHintsLeft,1,1,1,1));

                fTBN_DrawChannel = new TGTextButton(fHF_channels,"Draw",kB_drawChannel);
                fTBN_DrawChannel->Connect("Clicked()","GuiFrame",this,"OnDrawChannel()");
                fHF_channels->AddFrame(fTBN_DrawChannel,new TGLayoutHints(kLHintsRight|kLHintsCenterY,1,1,1,1));

                fTBN_ScanChannels = new TGTextButton(fHF_channels,"Scan",kB_scanChannel);
                fTBN_ScanChannels->Connect("Clicked()","GuiFrame",this,"OnScanChannels()");
                fHF_channels->AddFrame(fTBN_ScanChannels,new TGLayoutHints(kLHintsRight|kLHintsCenterY,1,1,1,1));
            fGF_result->AddFrame(fHF_channels,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,1,1,1,1));
            //--result: seperate line 5----------
            TGHorizontal3DLine *hsep_5 = new TGHorizontal3DLine(fGF_result, 350, 3);
            fGF_result->AddFrame(hsep_5, new TGLayoutHints(kLHintsTop|kLHintsExpandX, 1, 1, 1, 1));
            //--result: pedestals---
            TGHorizontalFrame *fHF_pedestal=new TGHorizontalFrame(fGF_result,350,25);
                TGLabel* fLResult_ped=new TGLabel(fHF_pedestal,"Pedestal:");
                fHF_pedestal->AddFrame(fLResult_ped,new TGLayoutHints(kLHintsLeft,1,1,1,1));
            fGF_result->AddFrame(fHF_pedestal,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
            //--result: calibration--
            TGHorizontalFrame *fHF_calibration=new TGHorizontalFrame(fGF_result,350,25);
                TGLabel* fLResult_calib=new TGLabel(fHF_calibration,"Calibration:");
                fHF_calibration->AddFrame(fLResult_calib,new TGLayoutHints(kLHintsLeft,1,1,1,1));
            fGF_result->AddFrame(fHF_calibration,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
            //--result: mips--
            TGHorizontalFrame *fHF_mips=new TGHorizontalFrame(fGF_result,350,25);
                TGLabel* fLResult_mips=new TGLabel(fHF_mips,"MIPs:");
                fHF_mips->AddFrame(fLResult_mips,new TGLayoutHints(kLHintsLeft,1,1,1,1));
            fGF_result->AddFrame(fHF_mips,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
        fGF_analyze->AddFrame(fGF_result,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
        //---RightFrame: Utility-------------
        TGGroupFrame *fGF_utility = new TGGroupFrame(fFScidataPanel,"Utility",kVerticalFrame);
        fGF_utility->SetTitlePos(TGGroupFrame::kCenter);
        fGF_utility->Resize(350,300);
            //--utility: low thresh--
            TGGroupFrame* fGF_lowthresh=new TGGroupFrame(fGF_utility,"LowThresh Encoding",kVerticalFrame);
                TGHorizontalFrame *fHF_low_encoding=new TGHorizontalFrame(fGF_lowthresh,350,25);
                    TGLabel *fLLow_value=new TGLabel(fHF_low_encoding,"Value:");
                    fNE_lowthresh=new TGNumberEntry(fHF_low_encoding,5,2,kNE_lowthresh,TGNumberFormat::kNESInteger,
                                                                                       TGNumberFormat::kNEANonNegative,
                                                                                       TGNumberFormat::kNELLimitMinMax,0,99);
                    TGLabel *fLLow_sigma=new TGLabel(fHF_low_encoding,"sigma");
                    fTBN_lowthreshEncoding=new TGTextButton(fHF_low_encoding,"Generate",kB_lowthreshEncoding);
                    fTBN_lowthreshEncoding->Connect("Clicked()","GuiFrame",this,"OnLowthreshEncoding()");
                    fHF_low_encoding->AddFrame(fLLow_value,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                    fHF_low_encoding->AddFrame(fNE_lowthresh,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                    fHF_low_encoding->AddFrame(fLLow_sigma,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                    fHF_low_encoding->AddFrame(fTBN_lowthreshEncoding,new TGLayoutHints(kLHintsRight,1,1,1,1));
                fGF_lowthresh->AddFrame(fHF_low_encoding,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));

                TGHorizontalFrame *fHF_low_input=new TGHorizontalFrame(fGF_lowthresh,350,25);
                    fFile_lowthreshOpen=new TGTextEntry(fHF_low_input,fTb_lowthreshOpen=new TGTextBuffer(200));
                    fTBN_lowthreshOpen=new TGTextButton(fHF_low_input,"Open..",kB_lowthreshOpen);
                    fTBN_lowthreshOpen->Connect("Clicked()","GuiFrame",this,"OnInputLowthresh()");
                    fHF_low_input->AddFrame(fFile_lowthreshOpen,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                    fHF_low_input->AddFrame(fTBN_lowthreshOpen,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fGF_lowthresh->AddFrame(fHF_low_input,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));

                TGHorizontalFrame *fHF_low_output=new TGHorizontalFrame(fGF_lowthresh,350,25);
                    fFile_lowthreshSave=new TGTextEntry(fHF_low_output,fTb_lowthreshSave=new TGTextBuffer(200));
                    fTBN_lowthreshSave=new TGTextButton(fHF_low_output,"Save..",kB_lowthreshSave);
                    fTBN_lowthreshSave->Connect("Clicked()","GuiFrame",this,"OnOutputLowthresh()");
                    fHF_low_output->AddFrame(fFile_lowthreshSave,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                    fHF_low_output->AddFrame(fTBN_lowthreshSave,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fGF_lowthresh->AddFrame(fHF_low_output,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
            fGF_utility->AddFrame(fGF_lowthresh,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
            //--utility: high thresh--
            TGGroupFrame* fGF_highthresh=new TGGroupFrame(fGF_utility,"HighThresh Encoding",kVerticalFrame);
                TGHorizontalFrame *fHF_high_encoding=new TGHorizontalFrame(fGF_highthresh,350,25);
                    TGLabel *fLHigh_value=new TGLabel(fHF_high_encoding,"Value(0x):");
                    fNE_highthresh=new TGNumberEntry(fHF_high_encoding,0,2,kNE_highthresh,TGNumberFormat::kNESHex,
                                                                                             TGNumberFormat::kNEANonNegative,
                                                                                             TGNumberFormat::kNELLimitMinMax,0,255);
                    fNE_highthresh->Resize(50,20);
                    fTBN_highthreshEncoding=new TGTextButton(fHF_high_encoding,"Generate",kB_highthreshEncoding);
                    fTBN_highthreshEncoding->Connect("Clicked()","GuiFrame",this,"OnHighthreshEncoding()");
                    fHF_high_encoding->AddFrame(fLHigh_value,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                    fHF_high_encoding->AddFrame(fNE_highthresh,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                    fHF_high_encoding->AddFrame(fTBN_highthreshEncoding,new TGLayoutHints(kLHintsRight,1,1,1,1));
                fGF_highthresh->AddFrame(fHF_high_encoding,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));

                TGHorizontalFrame *fHF_high_output=new TGHorizontalFrame(fGF_highthresh,350,25);
                    fFile_highthreshSave=new TGTextEntry(fHF_high_output,fTb_highthreshSave=new TGTextBuffer(200));
                    fTBN_highthreshSave=new TGTextButton(fHF_high_output,"Save..",kB_highthreshSave);
                    fTBN_highthreshSave->Connect("Clicked()","GuiFrame",this,"OnOutputHighthresh()");
                    fHF_high_output->AddFrame(fFile_highthreshSave,new TGLayoutHints(kLHintsLeft|kLHintsExpandX,1,1,1,1));
                    fHF_high_output->AddFrame(fTBN_highthreshSave,new TGLayoutHints(kLHintsLeft,1,1,1,1));
                fGF_highthresh->AddFrame(fHF_high_output,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
            fGF_utility->AddFrame(fGF_highthresh,new TGLayoutHints(kLHintsExpandX|kLHintsTop,1,1,1,1));
        fFScidataPanel->AddFrame(fGF_utility,new TGLayoutHints(kLHintsExpandX,1,1,1,1));
    //---Right Frame: SetFrame
    vsp_right->SetFrame(fFScidataPanel, kFALSE);
    /***********MainFrame END*******************************************************/
    SetResultState(false);

    SetWindowName("PSD Data Decoding");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    fFileCont->SetDefaultHeaders();
//		fFileCont->DisplayDirectory();
    //	fFileCont->AddFile("..");        // up level directory
    //	fFileCont->Resize();
    fFileCont->StopRefreshTimer();   // stop refreshing
    fFileCont->SetPageDimension(0, 0);
//		fFileCont->SetPageDimension(20, 10);
    fFileCont->SetViewMode(kLVDetails);
    //fFileCont->SetColHeaders("name");
    fFileCont->SetHeaders(2);
    fFileCont->SetColHeaders("     name                             ", "    title                                                            ");//don't delete the spaces

    Resize();
    CreateCanvas();
//-----------------------------------------------------------------

    fi_engineer=new TGFileInfo();
    fi_Sci_input=new TGFileInfo();
    fi_Sci_output=new TGFileInfo();
    fi_Sci_ref=new TGFileInfo();
    fi_lowthreshOpen=new TGFileInfo();
    fi_lowthreshSave=new TGFileInfo();
    fi_highthreshSave=new TGFileInfo();

    datafile = new TFile("myfile.root", "RECREATE");
    TH_init();
    ObjList = datafile->GetList();
    fFileCont->AddFile(datafile->GetName());
    fFileCont->Resize();
    processing = true;

    MakeTH1MenuList();
    RemoveMenuEntry("Delete", THmlist);
    RemoveMenuEntry("Dump", THmlist);
    RemoveMenuEntry("SetName", THmlist);
    RemoveMenuEntry("SetMaximum", THmlist);
    RemoveMenuEntry("SetMinimum", THmlist);
    RemoveMenuEntry("ShowBackground", THmlist);

    MakeTH2MenuList();
    RemoveMenuEntry("Delete", TH2mlist);
    RemoveMenuEntry("Dump", TH2mlist);
    RemoveMenuEntry("SetName", TH2mlist);
    RemoveMenuEntry("SetMaximum", TH2mlist);
    RemoveMenuEntry("SetMinimum", TH2mlist);
    RemoveMenuEntry("ShowBackground", TH2mlist);
//-----------------------------------------------------------------------------
}
//--------------------------------------------------------------------

GuiFrame::~GuiFrame()
{
    if(datafile) {
        datafile->Write(0,TObject::kOverwrite);
        delete datafile;
    }
    if(canpf)  delete canpf;

    delete fi_engineer;
    delete fi_Sci_input;
    delete fi_Sci_output;
    delete fi_Sci_ref;
    delete fi_lowthreshOpen;
    delete fi_highthreshSave;
    delete fi_lowthreshSave;
//------------------------------------------------------
    delete[] onevent;
    processing = false;
    if (TVirtualPadEditor::GetPadEditor(kFALSE) != 0)
        TVirtualPadEditor::Terminate();
//---------------------------------------------------------
}
void GuiFrame::OpenEngineerData()
{
    new TGFileDialog(gClient->GetRoot(),this,kFDOpen,fi_engineer);

    fFileEngineer->Clear();
    fFileEngineer->SetText(fi_engineer->fFilename);

    //ShowText(fTbEngineer->GetString());
}

void GuiFrame::CloseWindow()
{
     if (fClose)
     {
         DeleteWindow();
         gApplication->Terminate(0);
     }
   else {
      fClose = kTRUE;
      TTimer::SingleShot(150, "GuiFrame", this, "CloseWindow()");
   }

}

const char* GuiFrame::GetDrawOpt()
{
    int num =  fcomDrawOpt->GetSelected();
    switch (num)
    {
    case 0:
        return "";
        break;
    case 1:
        return "COLZ";
        break;
    case 2:
        return "LEGO";
        break;
    case 3:
        return "CONT";
        break;
    case 4:
        return "SURF";
        break;
    default:
        return "";
        break;
    }
}

void GuiFrame::ShowText(TGText *text)
{
    ClearTextView();
    fviewText->AddText(text);
    fviewText->Update();
    fviewText->ShowBottom();
}

void GuiFrame::ShowText(const char *text)
{
    ClearTextView();
    fviewText->AddLineFast(text);
    fviewText->Update();
    fviewText->ShowBottom();
}

void GuiFrame::CreateCanvas()
{
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetFrameFillStyle(0);
    gStyle->SetFrameFillStyle(1001);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFuncColor(kRed);
    gStyle->SetStatColor(kWhite);
    gStyle->SetStatBorderSize(1);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetLegendBorderSize(1);
    gStyle->SetFillColor(1);
    gStyle->SetPalette(1);
    gStyle->SetOptFit(0111);
//	gStyle->SetHistLineWidth(0.2);
    /*
        //canpf = new TCanvas("can_pf", "Process_Files",700,200,600,400);
        if(!canpf->GetShowToolBar())   canpf->ToggleToolBar();
        if(!canpf->GetShowEventStatus()) canpf->ToggleEventStatus();
        canpf->SetBorderMode(0);
        canpf->SetFillColor(0);
        canpf->SetFillStyle(0);
        canpf->SetFrameBorderMode(0);
        ((TRootCanvas *)canpf->GetCanvasImp())->DontCallClose();
        //remove some contextmenu items of TCanvas
        MakeTcMenuList();
        RemoveMenuEntry("SetEditable", TCmlist);
        RemoveMenuEntry("DrawClonePad", TCmlist);
        RemoveMenuEntry("DrawClone", TCmlist);
        RemoveMenuEntry("Dump", TCmlist);
        RemoveMenuEntry("SetName", TCmlist);
*/
        TH1F* hist=new TH1F("hist","hist",100,-50,50);
        hist->FillRandom("gaus",100);
        canpf->cd();
        hist->Draw();

}

void GuiFrame::ClearTextView()
{
    if(fviewText->ReturnLineCount()>=50)
    {
        fviewText->GetText()->DelLine(1);
    }
}

void GuiFrame::OnDoubleClick(TGLVEntry *f, Int_t btn)
{
    if (btn!=kButton1) return;
    f=(TGLVEntry *)fFileCont->GetLastActive();
    gVirtualX->SetCursor(fFileCont->GetId(),gVirtualX->CreateCursor(kWatch));

    TString name(f->GetTitle());
    const char* fname = (const char*)f->GetUserData();

    if (fname) {
        DisplayObject(fname,name);
    } else if (name.EndsWith(".root")) {
        DisplayFile(name);
    } else {
        DisplayDirectory(name);
    }
    gVirtualX->SetCursor(fFileCont->GetId(),gVirtualX->CreateCursor(kPointer));
}

void GuiFrame::DisplayFile(const TString &fname)
{
    // display content of ROOT file
    if(!datafile) return;

    fFileCont->RemoveAll();
    //	fFileCont->AddFile(gSystem->WorkingDirectory());
    fFileCont->SetPagePosition(0,0);
    //	fFileCont->SetHeaders(2);
    //	fFileCont->SetColHeaders("name", "title");

    TIter next(datafile->GetList());
    //TKey *key;
    TObject *key;

    while ((key=(TObject*)next())) {

        //if(!(key->IsA())->InheritsFrom("TH1")) continue;
        TString cname = key->GetTitle();
        TString name = key->GetName();
        TGLVEntry *entry = new TGLVEntry(fFileCont,name,cname);
        entry->SetSubnames(key->GetTitle());
        fFileCont->AddItem(entry);

        // user data is a filename
        entry->SetUserData((void*)StrDup(fname));
    }

    Resize();
}

void GuiFrame::DisplayObject(const TString& fname,const TString& name)
{
    // browse object located in file

    if(!datafile) return;
    TDirectory *sav = gDirectory;

    static TFile *file = 0;
    file = (TFile *) gROOT->GetListOfFiles()->FindObject(fname.Data());
    if(!file) return;

    TObject* obj = file->Get(name);
    if (obj) {
        if (!obj->IsFolder()) {
            canpf->cd(0);
            DrawObj(obj);
        } else obj->Print();
    }
    gDirectory = sav;
}

void GuiFrame::DisplayDirectory(const TString &fname)
{
    // display content of directory

    fFileCont->SetDefaultHeaders();
    gSystem->ChangeDirectory(fname);
    fFileCont->ChangeDirectory(fname);
    fFileCont->DisplayDirectory();
    //fFileCont->AddFile("..");  // up level directory
    Resize();
}

void GuiFrame::ImB_update()
{
    if(!canpf) return;
    if(!ObjListOK()) return;
//	int nx=1, ny=1;
//	GetDivPad(nx, ny);
    if(canpf->GetPad(1))
    {
        for(int i=1; i<30; i++)
        {
            if(canpf->GetPad(i))
            {
                canpf->cd(i);
                gPad->Modified();
                gPad->Update();
            }
        }
    }
    else
    {
        gPad->Modified();
        gPad->Update();
    }
}

void GuiFrame::ImB_next()
{
    DrawObj(kB_next);
    //if(!ObjListOK()) return;
    //do
    //{
    //	DrawObj(ObjList->After(objcurr));
    //}while(!(objcurr->IsA()->InheritsFrom("TH1")));
}

void GuiFrame::ImB_previous()
{
    DrawObj(kB_privious);
    //if(!ObjListOK()) return;
    //do
    //{
    //	DrawObj(ObjList->Before(objcurr));
    //}while(!(objcurr->IsA()->InheritsFrom("TH1")));
}

bool GuiFrame::ObjListOK()
{
    if( objcurr && ObjList )
    {
        return true;
    }
    else
    {
        return false;
    }

}

void GuiFrame::DrawObj(TObject *obj)
{
    if(!obj) return;
    if(!canpf) return;
    fHProgH->Reset();
    obj->Draw(GetDrawOpt());
    fHProgH->Increment(50);
    canpf->Update();
    fHProgH->Increment(50);
    objcurr = obj;
    fHProgH->Reset();
}

void GuiFrame::DrawObj(CMDIdentifiers id)
{
    if(!ObjListOK()) return;
    canpf->Clear();
    int nx=1, ny=1;
    GetDivPad(nx, ny);
    if(nx*ny>2)
    {
        canpf->Divide(nx, ny);
        canpf->SetLogy();
    }

    if(id == kB_privious)
    {
        for(int np=nx*ny; np>=1; np--)
        {
            canpf->cd(np);
            if(nx*ny>2) gPad->SetLogy();
            do
            {
                DrawObj(ObjList->Before(objcurr));
            }while(!(objcurr->IsA()->InheritsFrom("TH1")));

        }

    } else if(id == kB_next)
    {
        for(int np=1; np<=nx*ny; np++)
        {
            canpf->cd(np);
            if(nx*ny>2) gPad->SetLogy();
            do
            {
                DrawObj(ObjList->After(objcurr));
            }while(!(objcurr->IsA()->InheritsFrom("TH1")));

        }

    }

}

void GuiFrame::ImB_ResetAllTH()
{
    if(!ObjListOK()) return;
    TIter next(ObjList);
    TObject *obj;
    while(obj = next())
    {
        ResetTH(obj);
    }
    ImB_update();
}

void GuiFrame::ImB_ResetCurrTH()
{
    if(!ObjListOK()) return;
    ResetTH(objcurr);
    ImB_update();
}

void GuiFrame::ResetTH(TObject *obj)
{
    if (obj && (obj->IsA()->InheritsFrom("TH1")))
    {
        ((TH1 *)obj)->Reset();
    }
}

void GuiFrame::ImB_integral()
{
    if(!objcurr) return;
    if(!canpf) return;
    if(canpf->GetPad(1)) return;
    TString cname = objcurr->IsA()->GetName();
    if(!cname.Contains("TH1F", TString::kIgnoreCase)) return;
    ShowText("Select the left bin with mouse:");
    TMarker *p1 = (TMarker *)gPad->WaitPrimitive("TMarker");
    float left = 0.;
    if (p1)
    {
        p1->SetMarkerStyle(3);
        p1->SetMarkerColor(2);
        p1->SetMarkerSize(2.5);
        p1->Draw();
        left = p1->GetX();
    }

    TString text_xl = "X1= ";
    text_xl += left;
    ShowText(text_xl.Data());
    delete p1;
    ShowText("Select the right bin with mouse:");
    TMarker *p2 = (TMarker *)gPad->WaitPrimitive("TMarker");
    float right = 0.;
    if (p2)
    {
        p2->SetMarkerStyle(3);
        p2->SetMarkerColor(2);
        p2->SetMarkerSize(2.5);
        p2->Draw();
        right = p2->GetX();
    }

    TString text_xr = "X2= ";
    text_xr += right;
    ShowText(text_xr.Data());
    p2->Delete();
    TAxis *xaxis = ((TH1*)objcurr)->GetXaxis();
    if(left>right)
    {
        float xxt = left;
        left = right;
        right = xxt;
    }
    int xlbin = xaxis->FindBin(left);
    int xrbin = xaxis->FindBin(right);

    double numb = ((TH1*)objcurr)->Integral(xlbin, xrbin);
    float per=0.0;
    if(((TH1*)objcurr)->GetEntries()>0)
    {
        per = numb*100./(((TH1*)objcurr)->GetEntries());
    }
    TString text_inte = "Integral: ";
    text_inte += numb;
    ShowText(text_inte);
    TString text_per = "Percent: ";
    //TString data_per = "p";
    char s[32];
    sprintf(s, "%.5f", per);
    TString data_per = s;
    TString data_per1 = data_per(0, 4);
    text_per += data_per1;
    text_per += "%";
    ShowText(text_per.Data());

    TString stime = "Integral time: ";
    time_t nowtime;
    time(&nowtime);
    struct tm * timeinfo;
    timeinfo = localtime(&nowtime);
    TString ltime = asctime(timeinfo);
    stime += ltime(11, 8);
    ShowText(stime.Data());
    ShowText("   ");


//	p1->Delete();
//	p2->Delete();
}

void GuiFrame::GetDivPad(int &nx, int &ny)
{
    nx = 1;
    ny = 1;
    int n = fcomDivPad->GetSelected();
    switch(n)
    {
    case 0:
        nx = 1;
        ny = 1;
        break;
    case 1:
        nx = 2;
        ny = 2;
        break;
    case 2:
        nx = 2;
        ny = 3;
        break;
    case 3:
        nx = 3;
        ny = 3;
        break;
    default:
        break;
    }
}

void GuiFrame::RemoveMenuEntry(const char *menuTitle, TList *mlist)
{
    TClassMenuItem *menuItem = 0;
    TString itemTitle;

    if(!mlist) return;
    TIter next(mlist);
    while( menuItem = (TClassMenuItem *)next() )
    {
        itemTitle = menuItem->GetTitle();
        if( itemTitle == menuTitle)
        {
            delete menuItem;
            break;
        }
    }
}

void GuiFrame::MakeTcMenuList()
{
    TClass *c_TCanvas = gROOT->GetClass("TCanvas");
    if(!c_TCanvas) return;

    c_TCanvas->MakeCustomMenuList();
    TCmlist = c_TCanvas->GetMenuList();
}

void GuiFrame::MakeTH1MenuList()
{
    TClass *c_TH1 = gROOT->GetClass("TH1F");
    if(!c_TH1) return;

    c_TH1->MakeCustomMenuList();
    THmlist = c_TH1->GetMenuList();
}

void GuiFrame::MakeTH2MenuList()
{
    TClass *c_TH2 = gROOT->GetClass("TH2F");
    if(!c_TH2) return;

    c_TH2->MakeCustomMenuList();
    TH2mlist = c_TH2->GetMenuList();
}

//----------------------------------------------------------------------------------------
void GuiFrame::ProcessFiles()
{

    TString filename;
    filename=fi_engineer->fFilename;
    char buffer[400];
    sprintf(buffer,"Process Engineer Data:");
    ShowText(buffer);
    sprintf(buffer,"   %s",filename.Data());
    ShowText(buffer);

    DecodeEngineerData(filename);
}

void GuiFrame::Stop()
{
    processing = false;
    fBonefile->SetEnabled(true);
    //Pixel_t yellow;
    //gClient->GetColorByName("red", yellow);
    //fviewText->SetSelectFore(yellow);
    ShowText("   Decoding Stopped!");
    ShowText("Process End.");
//  triggercounttotal=0;

}

void GuiFrame::DecodeEngineerData(TString filename)
{
    fClose = kFALSE;
    char hname[200];
    char top[2];
//	char datatop[4];
    unsigned short event_num;
    int data_length,initNo;
    int Ic[4];
    int Temp[4];
    int addtemp;
    int FPGA_add[32][4];
    int FPGA_temp[32][4];
    int Ic_count_temp[4];
    int Temp_count_temp[4];
    int FPGA_count_temp[4];
    int time_second[4];
    int time_millisecond[4];
    int add[4];
    int i;
    char outinfo[256];

    FILE * fp = fopen(filename.Data(), "rb");

    if (fp == NULL)
    {
         int  buttons, retval;
         EMsgBoxIcon icontype = kMBIconStop;
         buttons = 0;
         gGuiFrame->Disconnect("CloseWindow()");
         filename.Append("can't be opened.");
         new TGMsgBox(gClient->GetRoot(), gGuiFrame,
                fTbEngineer->GetString(), filename,
                icontype, buttons, &retval);

         gGuiFrame->Connect("CloseWindow()", "GuiFrame", this, "CloseWindow()");

         return ;
    }

    data_length=200;
    onevent = new unsigned short[data_length];
    char onedata[400];

    processing = true;
    fBonefile->SetEnabled(false);
    ShowText("   Decoding Start...");

    event_num=0;
    initNo=0;

    addtemp=0;
  for(int i=0;i<32;i++)
  {
    for(int j=0;j<4;j++)
    {
        FPGA_add[i][j]=0;
        FPGA_temp[i][j]=0;
    }
  }
  for (i=0;i<4;i++)
  {
    add[i]=0;
    Ic[i]=0;
    Temp[i]=0;
    Ic_count_temp[i]=0;
    Temp_count_temp[i]=0;
    FPGA_count_temp[i]=0;
    time_second[i]=0;
    time_millisecond[i]=0;
  }

    int FILEID=0;//to determin scientific data or engineering data
  int FILEIDtmp1=0;
  int FILEIDtmp2=0;
    while(!feof(fp)&&processing)
    {
        if(fread(&top,sizeof(char),2,fp)==NULL) break;
    if((top[0]&0x00ff)==0xe2&&(top[1]&0x00ff)==0x25)//synchronous code
        {
            if(fread(&top,sizeof(char),2,fp)==NULL) break;
                else
                    {
                        FILEID=FILEIDHLAdd(top[0],top[1]);
                        //cout<<"FILEID = "<<FILEID<<endl;
                        switch(FILEID)
                        {
                            case 2067://1:for	scientific data
                                break;//1:end for	scientific data

                            case 2313://2:for telemetering of current of FEE at Positive(+)
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                data_length = HLAdd(top[0],top[1])+1;
                                //cout<<"data_length = "<<data_length<<endl;
                  if(fread(&onedata,sizeof(char),data_length,fp)==NULL) break;

                  time_second[0] =((onedata[2]&0xFF)<<24)+((onedata[3]&0xFF)<<16)+((onedata[4]&0xFF)<<8)+(onedata[5]&0xFF);
                  cout<<"time_second[0] = "<<time_second[0]<<endl;
                  time_millisecond[0] =((onedata[6]&0xFF)<<8)+(onedata[7]&0xFF);

                          add[0]=HLAddCurrent(onedata[8],onedata[9]);//PX
                          add[0]=add[0]*0.0004068*500;
                          sprintf(hname, "FEE_I_1");
                  HF1cont(hname,Ic[0],add[0]);
                          Ic[0]++;
                  if(Ic[0]>ED_xChan_maxNo) Ic[0]=0;
                  Ic_count_temp[0]++;
                          triggercounttotal1=Ic_count_temp[0];

                          add[2]=HLAddCurrent(onedata[36],onedata[37]);//PY
                          add[2]=add[2]*0.0004068*500;
                          sprintf(hname, "FEE_I_3");
                  HF1cont(hname,Ic[2],add[2]);
                          Ic[2]++;
                  if(Ic[2]>ED_xChan_maxNo) Ic[2]=0;
                  Ic_count_temp[2]++;
                          triggercounttotal3=Ic_count_temp[2];

                                break;//2:end for telemetering of current of FEE at Positive(+)

                            case 2330://3:for telemetering of temperature of FEE at Positive(+)
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                data_length = HLAdd(top[0],top[1])+1;
                                //cout<<"data_length = "<<data_length<<endl;
                  if(fread(&onedata,sizeof(char),data_length,fp)==NULL) break;

                  addtemp=HLAdd(onedata[8],onedata[9]);//PX
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_1");
                  HF1cont(hname,Temp[0],addtemp);
                  addtemp=HLAdd(onedata[10],onedata[11]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_2");
                  HF1cont(hname,Temp[0],addtemp);
                  addtemp=HLAdd(onedata[12],onedata[13]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_3");
                  HF1cont(hname,Temp[0],addtemp);
                  addtemp=HLAdd(onedata[14],onedata[15]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_4");
                  HF1cont(hname,Temp[0],addtemp);
                  Temp[0]++;
                  if(Temp[0]>ED_xChan_maxNo) Temp[0]=0;
                  Temp_count_temp[0]++;
                  triggercounttotal1=Temp_count_temp[0];

                  addtemp=HLAdd(onedata[240],onedata[241]);//PY
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_9");
                  HF1cont(hname,Temp[2],addtemp);
                  addtemp=HLAdd(onedata[242],onedata[243]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_10");
                  HF1cont(hname,Temp[2],addtemp);
                  addtemp=HLAdd(onedata[244],onedata[245]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_11");
                  HF1cont(hname,Temp[2],addtemp);
                  addtemp=HLAdd(onedata[246],onedata[247]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_12");
                  HF1cont(hname,Temp[2],addtemp);
                  Temp[2]++;
                  if(Temp[2]>ED_xChan_maxNo) Temp[2]=0;
                  Temp_count_temp[2]++;
                  triggercounttotal3=Temp_count_temp[2];

                                break;//3:for telemetering of temperature of FEE at Positive(+)
                            case 2348://4:for telemetering of Status of FEE at Positive(+)
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                data_length = HLAdd(top[0],top[1])+1;
                                //cout<<"data_length = "<<data_length<<endl;
                  if(fread(&onedata,sizeof(char),data_length,fp)==NULL) break;

 // PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX PX
                      //for status_reg[0]       PX
                      FPGA_add[0][0]=onedata[8]&0x3F;
                  sprintf(hname, "FPGA_1_1");	//FPGA_1_1 is for FEEID
                  HF1cont(hname,FPGA_temp[0][0],FPGA_add[0][0]);
                  FPGA_temp[0][0]++;
                  if(FPGA_temp[0][0]>ED_xChan_maxNo)    FPGA_temp[0][0]=0;

                  //for status_reg[1]
                  FPGA_add[1][0]=	onedata[9]&0x01;//
                  sprintf(hname, "FPGA_2_1_1");//FPGA_2_1 is for Auto-power off function, Bit0 if for trigger signal receive enable
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][0]=	(onedata[9]&0x02)>>1;//Bit1 if for operating mode
                  sprintf(hname, "FPGA_2_1_2");//FPGA_2_2 is for Auto-power off function, Bit1 if for operation mode
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][0]=	(onedata[9]&0x0C)>>2;//Bit2-3 if for data output mode
                  sprintf(hname, "FPGA_2_1_3");//FPGA_2_3 is for Auto-power off function, Bit2-3 if for data output mode
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][0]=	(onedata[9]&0x10)>>4;//Bit4 if for trigger signal input selection
                  sprintf(hname, "FPGA_2_1_4");//FPGA_2_4 is for Auto-power off function, Bit4 if for trigger signal input selection
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][0]=	(onedata[9]&0x40)>>6;//Bit6 if for power status display
                  sprintf(hname, "FPGA_2_1_5");//FPGA_2_5 is for Auto-power off function, Bit6 if for power status display
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][0]=	(onedata[9]&0x80)>>7;//Bit7 if for auto power enable switch
                  sprintf(hname, "FPGA_2_1_6");//FPGA_2_6 is for Auto-power off function, Bit7 if for auto power enable switch
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);

                  FPGA_temp[1][0]++;
                  if(FPGA_temp[1][0]>ED_xChan_maxNo)    FPGA_temp[1][0]=0;

                  //for status_reg[2]
                      FPGA_add[2][0]=onedata[10]&0xFF;
                  sprintf(hname, "FPGA_3_1");	//FPGA_3_1 is for peak delay
                  HF1cont(hname,FPGA_temp[2][0],FPGA_add[2][0]);
                  FPGA_temp[2][0]++;
                  if(FPGA_temp[2][0]>ED_xChan_maxNo)FPGA_temp[2][0]=0;

                  //for status_reg[3] and status_reg[4]
                      FPGA_add[3][0]=(onedata[11]&0xF0)>>4;
                  sprintf(hname, "FPGA_4_1");	//FPGA_4_1 is for Trigger status display
                  HF1cont(hname,FPGA_temp[3][0],FPGA_add[3][0]);
                  FPGA_temp[3][0]++;
                  if(FPGA_temp[3][0]>ED_xChan_maxNo)    FPGA_temp[3][0]=0;

                      FPGA_add[4][0]=FPGA4_HLAdd(onedata[11],onedata[12]);
                  sprintf(hname, "FPGA_5_1");	//FPGA_5_1 is for Trigger Number
                  HF1cont(hname,FPGA_temp[4][0],FPGA_add[4][0]);
                  FPGA_temp[4][0]++;
                  if(FPGA_temp[4][0]>ED_xChan_maxNo)    FPGA_temp[4][0]=0;

                  //for status_reg[5] and status_reg[6]
                      FPGA_add[5][0]=(onedata[13]&0xF0)>>4;
                  sprintf(hname, "FPGA_6_1");	//FPGA_6_1 is for times of current auto power-off
                  HF1cont(hname,FPGA_temp[5][0],FPGA_add[5][0]);
                  FPGA_temp[5][0]++;
                  if(FPGA_temp[5][0]>ED_xChan_maxNo)  FPGA_temp[5][0]=0;

                      FPGA_add[6][0]=FPGA4_HLAdd(onedata[13],onedata[14]);
                  sprintf(hname, "FPGA_7_1");	//FPGA_7_1 is for current detection value display
                  HF1cont(hname,FPGA_temp[6][0],FPGA_add[6][0]);
                  FPGA_temp[6][0]++;
                  if(FPGA_temp[6][0]>ED_xChan_maxNo)    FPGA_temp[6][0]=0;

                  //for status_reg[7]
                  FPGA_add[7][0]=(onedata[15]&0x0F);
                  sprintf(hname, "FPGA_8_1_1");	//FPGA_8_1_1 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[7][0],FPGA_add[7][0]);
                      FPGA_add[7][0]=(onedata[15]&0xF0)>>4;
                  sprintf(hname, "FPGA_8_1_2");	//FPGA_8_1_2 is for data management command injection
                  HF1cont(hname,FPGA_temp[7][0],FPGA_add[7][0]);
                  FPGA_temp[7][0]++;
                  if(FPGA_temp[7][0]>ED_xChan_maxNo)    FPGA_temp[7][0]=0;

                  //for status_reg[8] and status_reg[9]
                      FPGA_add[8][0]=FPGA8_HLAdd(onedata[16],onedata[17]);
                  sprintf(hname, "FPGA_9_1");	//FPGA_9_1 is for timer of data management command injection packages
                  HF1cont(hname,FPGA_temp[8][0],FPGA_add[8][0]);
                  FPGA_temp[8][0]++;
                  if(FPGA_temp[8][0]>ED_xChan_maxNo)   FPGA_temp[8][0]=0;

                  //for status_reg[10] and status_reg[11]
                      FPGA_add[10][0]=FPGA8_HLAdd(onedata[18],onedata[19]);
                  sprintf(hname, "FPGA_11_1");	//FPGA_11_1 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[10][0],FPGA_add[10][0]);
                  FPGA_temp[10][0]++;
                  if(FPGA_temp[10][0]>ED_xChan_maxNo)   FPGA_temp[10][0]=0;

                  //for status_reg[12] and status_reg[13]
                      FPGA_add[12][0]=FPGA4_HLAdd(onedata[20],onedata[21]);
                  sprintf(hname, "FPGA_13_1");	//FPGA_13_1 is for threshold of auto power off
                  HF1cont(hname,FPGA_temp[12][0],FPGA_add[12][0]);
                  FPGA_temp[12][0]++;
                  if(FPGA_temp[12][0]>ED_xChan_maxNo)   FPGA_temp[12][0]=0;

                  //for status_reg[14] and status_reg[15]
                      FPGA_add[14][0]=FPGA8_HLAdd(onedata[22],onedata[23]);
                  sprintf(hname, "FPGA_15_1");	//FPGA_15_1 is for timer of power-off
                  HF1cont(hname,FPGA_temp[14][0],FPGA_add[14][0]);
                  FPGA_temp[14][0]++;
                  if(FPGA_temp[14][0]>ED_xChan_maxNo)  FPGA_temp[14][0]=0;

                      //for status_reg[16]
                      FPGA_add[16][0]=(onedata[24]&0xF0)>>4;
                  sprintf(hname, "FPGA_17_1");	//FPGA_17_1 is for command and trigger channel display
                  HF1cont(hname,FPGA_temp[16][0],FPGA_add[16][0]);
                  FPGA_temp[16][0]++;
                  if(FPGA_temp[16][0]>ED_xChan_maxNo) FPGA_temp[16][0]=0;

                  //for status_reg[17] and status_reg[18]
                      FPGA_add[17][0]=FPGA8_HLAdd(onedata[25],onedata[26]);
                  sprintf(hname, "FPGA_18_1");	//FPGA_18_1 is for timer of odd erro
                  HF1cont(hname,FPGA_temp[17][0],FPGA_add[17][0]);
                  FPGA_temp[17][0]++;
                  if(FPGA_temp[17][0]>ED_xChan_maxNo)   FPGA_temp[17][0]=0;

                  //for status_reg[19] and status_reg[20]
                      FPGA_add[19][0]=FPGA8_HLAdd(onedata[27],onedata[28]);
                  sprintf(hname, "FPGA_20_1");	//FPGA_20_1 is for timer of accumative sum
                  HF1cont(hname,FPGA_temp[19][0],FPGA_add[19][0]);
                  FPGA_temp[19][0]++;
                  if(FPGA_temp[19][0]>ED_xChan_maxNo)   FPGA_temp[19][0]=0;

                  //for status_reg[21] and status_reg[22]
                      FPGA_add[21][0]=FPGA8_HLAdd(onedata[29],onedata[30]);
                  sprintf(hname, "FPGA_22_1");	//FPGA_22_1 is for timer of invalid data management injection command
                  HF1cont(hname,FPGA_temp[21][0],FPGA_add[21][0]);
                  FPGA_temp[21][0]++;
                  if(FPGA_temp[21][0]>ED_xChan_maxNo)   FPGA_temp[21][0]=0;

                  //for status_reg[23]
                      FPGA_add[23][0]=onedata[31]&0xFF;
                  sprintf(hname, "FPGA_24_1");	//FPGA_24_1 is for timer of high lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[23][0],FPGA_add[23][0]);
                  FPGA_temp[23][0]++;
                  if(FPGA_temp[23][0]>ED_xChan_maxNo)   FPGA_temp[23][0]=0;

                  //for status_reg[24]
                      FPGA_add[24][0]=onedata[32]&0xFF;
                  sprintf(hname, "FPGA_25_1");	//FPGA_24_1 is for timer of low lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[24][0],FPGA_add[24][0]);
                  FPGA_temp[24][0]++;
                  if(FPGA_temp[24][0]>ED_xChan_maxNo)   FPGA_temp[24][0]=0;

                  //for status_reg[25] and status_reg[26]
                      FPGA_add[25][0]=FPGA8_HLAdd(onedata[33],onedata[34]);
                  sprintf(hname, "FPGA_26_1");	//FPGA_26_1 is for scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[25][0],FPGA_add[25][0]);
                  FPGA_temp[25][0]++;
                  if(FPGA_temp[25][0]>ED_xChan_maxNo)   FPGA_temp[25][0]=0;

                  //for status_reg[27] and status_reg[28]
                      FPGA_add[27][0]=FPGA8_HLAdd(onedata[35],onedata[36]);
                  sprintf(hname, "FPGA_28_1");	//FPGA_28_1 is for accumulative sum of scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[27][0],FPGA_add[27][0]);
                  FPGA_temp[27][0]++;
                  if(FPGA_temp[27][0]>ED_xChan_maxNo)   FPGA_temp[27][0]=0;

                  //for status_reg[29] and status_reg[30]
                      FPGA_add[29][0]=FPGA8_HLAdd(onedata[37],onedata[38]);
                  sprintf(hname, "FPGA_30_1");	//FPGA_30_1 is for timer of trigger width erro
                  HF1cont(hname,FPGA_temp[29][0],FPGA_add[29][0]);
                  FPGA_temp[29][0]++;
                  if(FPGA_temp[29][0]>ED_xChan_maxNo)   FPGA_temp[29][0]=0;
                  FPGA_count_temp[0]++;
                  triggercounttotal1=FPGA_count_temp[0];

//PYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPYPY

                      //for status_reg[0] PY
                      FPGA_add[0][2]=onedata[296]&0x3F;
                  sprintf(hname, "FPGA_1_3");	//FPGA_1_3 is for FEEID
                  HF1cont(hname,FPGA_temp[0][2],FPGA_add[0][2]);
                  FPGA_temp[0][2]++;
                  if(FPGA_temp[0][2]>ED_xChan_maxNo)    FPGA_temp[0][2]=0;

                  //for status_reg[1]
                  FPGA_add[1][2]=	onedata[297]&0x01;//
                  sprintf(hname, "FPGA_2_3_1");//FPGA_2_1 is for Auto-power off function, Bit0 if for trigger signal receive enable
                  HF1cont(hname,FPGA_temp[1][2],FPGA_add[1][2]);
                  FPGA_add[1][2]=	(onedata[297]&0x02)>>1;//Bit1 if for operating mode
                  sprintf(hname, "FPGA_2_3_2");//FPGA_2_2 is for Auto-power off function, Bit1 if for operation mode
                  HF1cont(hname,FPGA_temp[1][2],FPGA_add[1][2]);
                  FPGA_add[1][2]=	(onedata[297]&0x0C)>>2;//Bit2-3 if for data output mode
                  sprintf(hname, "FPGA_2_3_3");//FPGA_2_3 is for Auto-power off function, Bit2-3 if for data output mode
                  HF1cont(hname,FPGA_temp[1][2],FPGA_add[1][2]);
                  FPGA_add[1][2]=	(onedata[297]&0x10)>>4;//Bit4 if for trigger signal input selection
                  sprintf(hname, "FPGA_2_3_4");//FPGA_2_4 is for Auto-power off function, Bit4 if for trigger signal input selection
                  HF1cont(hname,FPGA_temp[1][2],FPGA_add[1][2]);
                  FPGA_add[1][2]=	(onedata[297]&0x40)>>6;//Bit6 if for power status display
                  sprintf(hname, "FPGA_2_3_5");//FPGA_2_5 is for Auto-power off function, Bit6 if for power status display
                  HF1cont(hname,FPGA_temp[1][2],FPGA_add[1][2]);
                  FPGA_add[1][2]=	(onedata[297]&0x80)>>7;//Bit7 if for auto power enable switch
                  sprintf(hname, "FPGA_2_3_6");//FPGA_2_6 is for Auto-power off function, Bit7 if for auto power enable switch
                  HF1cont(hname,FPGA_temp[1][2],FPGA_add[1][2]);

                  FPGA_temp[1][2]++;
                  if(FPGA_temp[1][2]>ED_xChan_maxNo)    FPGA_temp[1][2]=0;

                  //for status_reg[2]
                      FPGA_add[2][2]=onedata[298]&0xFF;
                  sprintf(hname, "FPGA_3_3");	//FPGA_3_3 is for peak delay
                  HF1cont(hname,FPGA_temp[2][2],FPGA_add[2][2]);
                  FPGA_temp[2][2]++;
                  if(FPGA_temp[2][2]>ED_xChan_maxNo)FPGA_temp[2][2]=0;

                  //for status_reg[3] and status_reg[4]
                      FPGA_add[3][2]=(onedata[299]&0xF0)>>4;
                  sprintf(hname, "FPGA_4_3");	//FPGA_4_3 is for Trigger status display
                  HF1cont(hname,FPGA_temp[3][2],FPGA_add[3][2]);
                  FPGA_temp[3][2]++;
                  if(FPGA_temp[3][2]>ED_xChan_maxNo)    FPGA_temp[3][2]=0;

                      FPGA_add[4][2]=FPGA4_HLAdd(onedata[299],onedata[300]);
                  sprintf(hname, "FPGA_5_3");	//FPGA_5_3 is for Trigger Number
                  HF1cont(hname,FPGA_temp[4][2],FPGA_add[4][2]);
                  FPGA_temp[4][2]++;
                  if(FPGA_temp[4][2]>ED_xChan_maxNo)    FPGA_temp[4][2]=0;

                  //for status_reg[5] and status_reg[6]
                      FPGA_add[5][2]=(onedata[301]&0xF0)>>4;
                  sprintf(hname, "FPGA_6_3");	//FPGA_6_3 is for times of current auto power-off
                  HF1cont(hname,FPGA_temp[5][2],FPGA_add[5][2]);
                  FPGA_temp[5][2]++;
                  if(FPGA_temp[5][2]>ED_xChan_maxNo)  FPGA_temp[5][2]=0;

                      FPGA_add[6][2]=FPGA4_HLAdd(onedata[301],onedata[302]);
                  sprintf(hname, "FPGA_7_3");	//FPGA_7_3 is for current detection value display
                  HF1cont(hname,FPGA_temp[6][2],FPGA_add[6][2]);
                  FPGA_temp[6][2]++;
                  if(FPGA_temp[6][2]>ED_xChan_maxNo)    FPGA_temp[6][2]=0;

                  //for status_reg[7]
                  FPGA_add[7][2]=(onedata[303]&0x0F);
                  sprintf(hname, "FPGA_8_3_1");	//FPGA_8_3_1 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[7][2],FPGA_add[7][2]);
                      FPGA_add[7][2]=(onedata[303]&0xF0)>>4;
                  sprintf(hname, "FPGA_8_3_2");	//FPGA_8_3_2 is for data management command injection
                  HF1cont(hname,FPGA_temp[7][2],FPGA_add[7][2]);
                  FPGA_temp[7][2]++;
                  if(FPGA_temp[7][2]>ED_xChan_maxNo)    FPGA_temp[7][2]=0;

                  //for status_reg[8] and status_reg[9]
                      FPGA_add[8][2]=FPGA8_HLAdd(onedata[304],onedata[305]);
                  sprintf(hname, "FPGA_9_3");	//FPGA_9_3 is for timer of data management command injection packages
                  HF1cont(hname,FPGA_temp[8][2],FPGA_add[8][2]);
                  FPGA_temp[8][2]++;
                  if(FPGA_temp[8][2]>ED_xChan_maxNo)   FPGA_temp[8][2]=0;

                  //for status_reg[10] and status_reg[11]
                      FPGA_add[10][2]=FPGA8_HLAdd(onedata[306],onedata[307]);
                  sprintf(hname, "FPGA_11_3");	//FPGA_11_3 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[10][2],FPGA_add[10][2]);
                  FPGA_temp[10][2]++;
                  if(FPGA_temp[10][2]>ED_xChan_maxNo)   FPGA_temp[10][2]=0;

                  //for status_reg[12] and status_reg[13]
                      FPGA_add[12][2]=FPGA4_HLAdd(onedata[308],onedata[309]);
                  sprintf(hname, "FPGA_13_3");	//FPGA_13_3 is for threshold of auto power off
                  HF1cont(hname,FPGA_temp[12][2],FPGA_add[12][2]);
                  FPGA_temp[12][2]++;
                  if(FPGA_temp[12][2]>ED_xChan_maxNo)   FPGA_temp[12][2]=0;

                  //for status_reg[14] and status_reg[15]
                      FPGA_add[14][2]=FPGA8_HLAdd(onedata[310],onedata[311]);
                  sprintf(hname, "FPGA_15_3");	//FPGA_15_3 is for timer of power-off
                  HF1cont(hname,FPGA_temp[14][2],FPGA_add[14][2]);
                  FPGA_temp[14][2]++;
                  if(FPGA_temp[14][2]>ED_xChan_maxNo)  FPGA_temp[14][2]=0;

                      //for status_reg[16]
                      FPGA_add[16][2]=(onedata[312]&0xF0)>>4;
                  sprintf(hname, "FPGA_17_3");	//FPGA_17_3 is for command and trigger channel display
                  HF1cont(hname,FPGA_temp[16][2],FPGA_add[16][2]);
                  FPGA_temp[16][2]++;
                  if(FPGA_temp[16][2]>ED_xChan_maxNo) FPGA_temp[16][2]=0;

                  //for status_reg[17] and status_reg[18]
                      FPGA_add[17][2]=FPGA8_HLAdd(onedata[313],onedata[314]);
                  sprintf(hname, "FPGA_18_3");	//FPGA_18_3 is for timer of odd erro
                  HF1cont(hname,FPGA_temp[17][2],FPGA_add[17][2]);
                  FPGA_temp[17][2]++;
                  if(FPGA_temp[17][2]>ED_xChan_maxNo)   FPGA_temp[17][2]=0;

                  //for status_reg[19] and status_reg[20]
                      FPGA_add[19][2]=FPGA8_HLAdd(onedata[315],onedata[316]);
                  sprintf(hname, "FPGA_20_3");	//FPGA_20_3 is for timer of accumative sum
                  HF1cont(hname,FPGA_temp[19][2],FPGA_add[19][2]);
                  FPGA_temp[19][2]++;
                  if(FPGA_temp[19][2]>ED_xChan_maxNo)   FPGA_temp[19][2]=0;

                  //for status_reg[21] and status_reg[22]
                      FPGA_add[21][2]=FPGA8_HLAdd(onedata[317],onedata[318]);
                  sprintf(hname, "FPGA_22_3");	//FPGA_22_3 is for timer of invalid data management injection command
                  HF1cont(hname,FPGA_temp[21][2],FPGA_add[21][2]);
                  FPGA_temp[21][2]++;
                  if(FPGA_temp[21][2]>ED_xChan_maxNo)   FPGA_temp[21][2]=0;

                  //for status_reg[23]
                      FPGA_add[23][2]=onedata[319]&0xFF;
                  sprintf(hname, "FPGA_24_3");	//FPGA_24_3 is for timer of high lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[23][2],FPGA_add[23][2]);
                  FPGA_temp[23][2]++;
                  if(FPGA_temp[23][2]>ED_xChan_maxNo)   FPGA_temp[23][2]=0;

                  //for status_reg[24]
                      FPGA_add[24][2]=onedata[320]&0xFF;
                  sprintf(hname, "FPGA_25_3");	//FPGA_24_3 is for timer of low lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[24][2],FPGA_add[24][2]);
                  FPGA_temp[24][2]++;
                  if(FPGA_temp[24][2]>ED_xChan_maxNo)   FPGA_temp[24][2]=0;

                  //for status_reg[25] and status_reg[26]
                      FPGA_add[25][2]=FPGA8_HLAdd(onedata[321],onedata[322]);
                  sprintf(hname, "FPGA_26_3");	//FPGA_26_3 is for scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[25][2],FPGA_add[25][2]);
                  FPGA_temp[25][2]++;
                  if(FPGA_temp[25][2]>ED_xChan_maxNo)   FPGA_temp[25][2]=0;

                  //for status_reg[27] and status_reg[28]
                      FPGA_add[27][2]=FPGA8_HLAdd(onedata[323],onedata[324]);
                  sprintf(hname, "FPGA_28_3");	//FPGA_28_3 is for accumulative sum of scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[27][2],FPGA_add[27][2]);
                  FPGA_temp[27][2]++;
                  if(FPGA_temp[27][2]>ED_xChan_maxNo)   FPGA_temp[27][2]=0;

                  //for status_reg[29] and status_reg[30]
                      FPGA_add[29][2]=FPGA8_HLAdd(onedata[325],onedata[326]);
                  sprintf(hname, "FPGA_30_3");	//FPGA_30_3 is for timer of trigger width erro
                  HF1cont(hname,FPGA_temp[29][2],FPGA_add[29][2]);
                  FPGA_temp[29][2]++;
                  if(FPGA_temp[29][2]>ED_xChan_maxNo)   FPGA_temp[29][2]=0;

                  FPGA_count_temp[2]++;
                  triggercounttotal3=FPGA_count_temp[2];

                                break;//4:end for telemetering of Status of FEE at Positive(+)
                            case 2367://5:for telemetering of current of FEE at Negative(-)
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                data_length = HLAdd(top[0],top[1])+1;
                                //cout<<"data_length = "<<data_length<<endl;
                  if(fread(&onedata,sizeof(char),data_length,fp)==NULL) break;

                          add[1]=HLAddCurrent(onedata[8],onedata[9]);//NX
                          add[1]=add[1]*0.0004068*500;
                          sprintf(hname, "FEE_I_2");
                  HF1cont(hname,Ic[1],add[1]);
                          Ic[1]++;
                          triggercounttotal2=Ic[1];
                  if(Ic[1]>ED_xChan_maxNo) Ic[1]=0;

                          add[3]=HLAddCurrent(onedata[36],onedata[37]);//NY
                          add[3]=add[3]*0.0004068*500;
                          sprintf(hname, "FEE_I_4");
                  HF1cont(hname,Ic[3],add[3]);
                          Ic[3]++;
                          triggercounttotal4=Ic[3];
                  if(Ic[3]>ED_xChan_maxNo) Ic[3]=0;

                                break;//5:end for telemetering of current of FEE at Negative(-)
                            case 2383://6:for telemetering of temperature of FEE at Negative(-)
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                data_length = HLAdd(top[0],top[1])+1;
                                //cout<<"data_length = "<<data_length<<endl;
                  if(fread(&onedata,sizeof(char),data_length,fp)==NULL) break;

                  addtemp=HLAdd(onedata[8],onedata[9]);//NX
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_5");
                  HF1cont(hname,Temp[1],addtemp);
                  addtemp=HLAdd(onedata[10],onedata[11]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_6");
                  HF1cont(hname,Temp[1],addtemp);
                  addtemp=HLAdd(onedata[12],onedata[13]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_7");
                  HF1cont(hname,Temp[1],addtemp);
                  addtemp=HLAdd(onedata[14],onedata[15]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_8");
                  HF1cont(hname,Temp[1],addtemp);
                  Temp[1]++;
                  if(Temp[1]>ED_xChan_maxNo) Temp[1]=0;
                  Temp_count_temp[1]++;
                  triggercounttotal2=Temp_count_temp[1];

                  addtemp=HLAdd(onedata[240],onedata[241]);//NY
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_13");
                  HF1cont(hname,Temp[3],addtemp);
                  addtemp=HLAdd(onedata[242],onedata[243]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_14");
                  HF1cont(hname,Temp[3],addtemp);
                  addtemp=HLAdd(onedata[244],onedata[245]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_15");
                  HF1cont(hname,Temp[3],addtemp);
                  addtemp=HLAdd(onedata[246],onedata[247]);
                  addtemp=TempChange(addtemp);
                  sprintf(hname, "T_16");
                  HF1cont(hname,Temp[3],addtemp);
                  Temp[3]++;
                  if(Temp[3]>ED_xChan_maxNo) Temp[3]=0;
                  Temp_count_temp[3]++;
                  triggercounttotal4=Temp_count_temp[3];
                                break;//6:end for telemetering of temperature of FEE at Negative(-)

                            case 2396://7:for telemetering of Status of FEE at Negative(-)
                            if(fread(&top,sizeof(char),2,fp)==NULL) break; //NX
                                if(fread(&top,sizeof(char),2,fp)==NULL) break;
                                data_length = HLAdd(top[0],top[1])+1;
                                //cout<<"data_length = "<<data_length<<endl;
                  if(fread(&onedata,sizeof(char),data_length,fp)==NULL) break;

// NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX NX


                      //for status_reg[0]       NX
                      FPGA_add[0][1]=onedata[8]&0x3F;
                  sprintf(hname, "FPGA_1_2");	//FPGA_1_2 is for FEEID
                  HF1cont(hname,FPGA_temp[0][1],FPGA_add[0][1]);
                  FPGA_temp[0][1]++;
                  if(FPGA_temp[0][1]>ED_xChan_maxNo)    FPGA_temp[0][1]=0;

                  //for status_reg[1]
                  FPGA_add[1][1]=	onedata[9]&0x01;//
                  sprintf(hname, "FPGA_2_2_1");//FPGA_2_1 is for Auto-power off function, Bit0 if for trigger signal receive enable
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][1]=	(onedata[9]&0x02)>>1;//Bit1 if for operating mode
                  sprintf(hname, "FPGA_2_2_2");//FPGA_2_2 is for Auto-power off function, Bit1 if for operation mode
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][1]=	(onedata[9]&0x0C)>>2;//Bit2-3 if for data output mode
                  sprintf(hname, "FPGA_2_2_3");//FPGA_2_3 is for Auto-power off function, Bit2-3 if for data output mode
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][1]=	(onedata[9]&0x10)>>4;//Bit4 if for trigger signal input selection
                  sprintf(hname, "FPGA_2_2_4");//FPGA_2_4 is for Auto-power off function, Bit4 if for trigger signal input selection
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][1]=	(onedata[9]&0x40)>>6;//Bit6 if for power status display
                  sprintf(hname, "FPGA_2_2_5");//FPGA_2_5 is for Auto-power off function, Bit6 if for power status display
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);
                  FPGA_add[1][1]=	(onedata[9]&0x80)>>7;//Bit7 if for auto power enable switch
                  sprintf(hname, "FPGA_2_2_6");//FPGA_2_6 is for Auto-power off function, Bit7 if for auto power enable switch
                  HF1cont(hname,FPGA_temp[1][0],FPGA_add[1][0]);

                  FPGA_temp[1][1]++;
                  if(FPGA_temp[1][1]>ED_xChan_maxNo)    FPGA_temp[1][1]=0;

                  //for status_reg[2]
                      FPGA_add[2][1]=onedata[10]&0xFF;
                  sprintf(hname, "FPGA_3_2");	//FPGA_3_2 is for peak delay
                  HF1cont(hname,FPGA_temp[2][1],FPGA_add[2][1]);
                  FPGA_temp[2][1]++;
                  if(FPGA_temp[2][1]>ED_xChan_maxNo)FPGA_temp[2][1]=0;

                  //for status_reg[3] and status_reg[4]
                      FPGA_add[3][1]=(onedata[11]&0xF0)>>4;
                  sprintf(hname, "FPGA_4_2");	//FPGA_4_2 is for Trigger status display
                  HF1cont(hname,FPGA_temp[3][1],FPGA_add[3][1]);
                  FPGA_temp[3][1]++;
                  if(FPGA_temp[3][1]>ED_xChan_maxNo)    FPGA_temp[3][1]=0;

                      FPGA_add[4][1]=FPGA4_HLAdd(onedata[11],onedata[12]);
                  sprintf(hname, "FPGA_5_2");	//FPGA_5_2 is for Trigger Number
                  HF1cont(hname,FPGA_temp[4][1],FPGA_add[4][1]);
                  FPGA_temp[4][1]++;
                  if(FPGA_temp[4][1]>ED_xChan_maxNo)    FPGA_temp[4][1]=0;

                  //for status_reg[5] and status_reg[6]
                      FPGA_add[5][1]=(onedata[13]&0xF0)>>4;
                  sprintf(hname, "FPGA_6_2");	//FPGA_6_2 is for times of current auto power-off
                  HF1cont(hname,FPGA_temp[5][1],FPGA_add[5][1]);
                  FPGA_temp[5][1]++;
                  if(FPGA_temp[5][1]>ED_xChan_maxNo)  FPGA_temp[5][1]=0;

                      FPGA_add[6][1]=FPGA4_HLAdd(onedata[13],onedata[14]);
                  sprintf(hname, "FPGA_7_2");	//FPGA_7_2 is for current detection value display
                  HF1cont(hname,FPGA_temp[6][1],FPGA_add[6][1]);
                  FPGA_temp[6][1]++;
                  if(FPGA_temp[6][1]>ED_xChan_maxNo)    FPGA_temp[6][1]=0;

                  //for status_reg[7]
                  FPGA_add[7][1]=(onedata[15]&0x0F);
                  sprintf(hname, "FPGA_8_2_1");	//FPGA_8_2_1 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[7][1],FPGA_add[7][1]);
                      FPGA_add[7][1]=(onedata[15]&0xF0)>>4;
                  sprintf(hname, "FPGA_8_2_2");	//FPGA_8_2_2 is for data management command injection
                  HF1cont(hname,FPGA_temp[7][1],FPGA_add[7][1]);
                  FPGA_temp[7][1]++;
                  if(FPGA_temp[7][1]>ED_xChan_maxNo)    FPGA_temp[7][1]=0;

                  //for status_reg[8] and status_reg[9]
                      FPGA_add[8][1]=FPGA8_HLAdd(onedata[16],onedata[17]);
                  sprintf(hname, "FPGA_9_2");	//FPGA_9_2 is for timer of data management command injection packages
                  HF1cont(hname,FPGA_temp[8][1],FPGA_add[8][1]);
                  FPGA_temp[8][1]++;
                  if(FPGA_temp[8][1]>ED_xChan_maxNo)   FPGA_temp[8][1]=0;

                  //for status_reg[10] and status_reg[11]
                      FPGA_add[10][1]=FPGA8_HLAdd(onedata[18],onedata[19]);
                  sprintf(hname, "FPGA_11_2");	//FPGA_11_2 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[10][1],FPGA_add[10][1]);
                  FPGA_temp[10][1]++;
                  if(FPGA_temp[10][1]>ED_xChan_maxNo)   FPGA_temp[10][1]=0;

                  //for status_reg[12] and status_reg[13]
                      FPGA_add[12][1]=FPGA4_HLAdd(onedata[20],onedata[21]);
                  sprintf(hname, "FPGA_13_2");	//FPGA_13_2 is for threshold of auto power off
                  HF1cont(hname,FPGA_temp[12][1],FPGA_add[12][1]);
                  FPGA_temp[12][1]++;
                  if(FPGA_temp[12][1]>ED_xChan_maxNo)   FPGA_temp[12][1]=0;

                  //for status_reg[14] and status_reg[15]
                      FPGA_add[14][1]=FPGA8_HLAdd(onedata[22],onedata[23]);
                  sprintf(hname, "FPGA_15_2");	//FPGA_15_2 is for timer of power-off
                  HF1cont(hname,FPGA_temp[14][1],FPGA_add[14][1]);
                  FPGA_temp[14][1]++;
                  if(FPGA_temp[14][1]>ED_xChan_maxNo)  FPGA_temp[14][1]=0;

                      //for status_reg[16]
                      FPGA_add[16][1]=(onedata[24]&0xF0)>>4;
                  sprintf(hname, "FPGA_17_2");	//FPGA_17_2 is for command and trigger channel display
                  HF1cont(hname,FPGA_temp[16][1],FPGA_add[16][1]);
                  FPGA_temp[16][1]++;
                  if(FPGA_temp[16][1]>ED_xChan_maxNo) FPGA_temp[16][1]=0;

                  //for status_reg[17] and status_reg[18]
                      FPGA_add[17][1]=FPGA8_HLAdd(onedata[25],onedata[26]);
                  sprintf(hname, "FPGA_18_2");	//FPGA_18_2 is for timer of odd erro
                  HF1cont(hname,FPGA_temp[17][1],FPGA_add[17][1]);
                  FPGA_temp[17][1]++;
                  if(FPGA_temp[17][1]>ED_xChan_maxNo)   FPGA_temp[17][1]=0;

                  //for status_reg[19] and status_reg[20]
                      FPGA_add[19][1]=FPGA8_HLAdd(onedata[27],onedata[28]);
                  sprintf(hname, "FPGA_20_2");	//FPGA_20_2 is for timer of accumative sum
                  HF1cont(hname,FPGA_temp[19][1],FPGA_add[19][1]);
                  FPGA_temp[19][1]++;
                  if(FPGA_temp[19][1]>ED_xChan_maxNo)   FPGA_temp[19][1]=0;

                  //for status_reg[21] and status_reg[22]
                      FPGA_add[21][1]=FPGA8_HLAdd(onedata[29],onedata[30]);
                  sprintf(hname, "FPGA_22_2");	//FPGA_22_2 is for timer of invalid data management injection command
                  HF1cont(hname,FPGA_temp[21][1],FPGA_add[21][1]);
                  FPGA_temp[21][1]++;
                  if(FPGA_temp[21][1]>ED_xChan_maxNo)   FPGA_temp[21][1]=0;

                  //for status_reg[23]
                      FPGA_add[23][1]=onedata[31]&0xFF;
                  sprintf(hname, "FPGA_24_2");	//FPGA_24_2 is for timer of high lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[23][1],FPGA_add[23][1]);
                  FPGA_temp[23][1]++;
                  if(FPGA_temp[23][1]>ED_xChan_maxNo)   FPGA_temp[23][1]=0;

                  //for status_reg[24]
                      FPGA_add[24][1]=onedata[32]&0xFF;
                  sprintf(hname, "FPGA_25_2");	//FPGA_24_2 is for timer of low lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[24][1],FPGA_add[24][1]);
                  FPGA_temp[24][1]++;
                  if(FPGA_temp[24][1]>ED_xChan_maxNo)   FPGA_temp[24][1]=0;

                  //for status_reg[25] and status_reg[26]
                      FPGA_add[25][1]=FPGA8_HLAdd(onedata[33],onedata[34]);
                  sprintf(hname, "FPGA_26_2");	//FPGA_26_2 is for scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[25][1],FPGA_add[25][1]);
                  FPGA_temp[25][1]++;
                  if(FPGA_temp[25][1]>ED_xChan_maxNo)   FPGA_temp[25][1]=0;

                  //for status_reg[27] and status_reg[28]
                      FPGA_add[27][1]=FPGA8_HLAdd(onedata[35],onedata[36]);
                  sprintf(hname, "FPGA_28_2");	//FPGA_28_2 is for accumulative sum of scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[27][1],FPGA_add[27][1]);
                  FPGA_temp[27][1]++;
                  if(FPGA_temp[27][1]>ED_xChan_maxNo)   FPGA_temp[27][1]=0;

                  //for status_reg[29] and status_reg[30]
                      FPGA_add[29][1]=FPGA8_HLAdd(onedata[37],onedata[38]);
                  sprintf(hname, "FPGA_30_2");	//FPGA_30_2 is for timer of trigger width erro
                  HF1cont(hname,FPGA_temp[29][1],FPGA_add[29][1]);
                  FPGA_temp[29][1]++;
                  if(FPGA_temp[29][1]>ED_xChan_maxNo)   FPGA_temp[29][1]=0;

                  FPGA_count_temp[1]++;
                  triggercounttotal2=FPGA_count_temp[1];

 //NYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNYNY


                      //for status_reg[0] NY
                      FPGA_add[0][3]=onedata[296]&0x3F;
                  sprintf(hname, "FPGA_1_4");	//FPGA_1_4 is for FEEID
                  HF1cont(hname,FPGA_temp[0][3],FPGA_add[0][3]);
                  FPGA_temp[0][3]++;
                  if(FPGA_temp[0][3]>ED_xChan_maxNo)    FPGA_temp[0][3]=0;

                  //for status_reg[1]
                  FPGA_add[1][3]=	onedata[297]&0x01;//
                  sprintf(hname, "FPGA_2_4_1");//FPGA_2_1 is for Auto-power off function, Bit0 if for trigger signal receive enable
                  HF1cont(hname,FPGA_temp[1][3],FPGA_add[1][3]);
                  FPGA_add[1][3]=	(onedata[297]&0x02)>>1;//Bit1 if for operating mode
                  sprintf(hname, "FPGA_2_4_2");//FPGA_2_2 is for Auto-power off function, Bit1 if for operation mode
                  HF1cont(hname,FPGA_temp[1][3],FPGA_add[1][3]);
                  FPGA_add[1][3]=	(onedata[297]&0x0C)>>2;//Bit2-3 if for data output mode
                  sprintf(hname, "FPGA_2_4_3");//FPGA_2_3 is for Auto-power off function, Bit2-3 if for data output mode
                  HF1cont(hname,FPGA_temp[1][3],FPGA_add[1][3]);
                  FPGA_add[1][3]=	(onedata[297]&0x10)>>4;//Bit4 if for trigger signal input selection
                  sprintf(hname, "FPGA_2_4_4");//FPGA_2_4 is for Auto-power off function, Bit4 if for trigger signal input selection
                  HF1cont(hname,FPGA_temp[1][3],FPGA_add[1][3]);
                  FPGA_add[1][3]=	(onedata[297]&0x40)>>6;//Bit6 if for power status display
                  sprintf(hname, "FPGA_2_4_5");//FPGA_2_5 is for Auto-power off function, Bit6 if for power status display
                  HF1cont(hname,FPGA_temp[1][3],FPGA_add[1][3]);
                  FPGA_add[1][3]=	(onedata[297]&0x80)>>7;//Bit7 if for auto power enable switch
                  sprintf(hname, "FPGA_2_4_6");//FPGA_2_6 is for Auto-power off function, Bit7 if for auto power enable switch
                  HF1cont(hname,FPGA_temp[1][3],FPGA_add[1][3]);

                  FPGA_temp[1][3]++;
                  if(FPGA_temp[1][3]>ED_xChan_maxNo)    FPGA_temp[1][3]=0;

                  //for status_reg[2]
                      FPGA_add[2][3]=onedata[298]&0xFF;
                  sprintf(hname, "FPGA_3_4");	//FPGA_3_4 is for peak delay
                  HF1cont(hname,FPGA_temp[2][3],FPGA_add[2][3]);
                  FPGA_temp[2][3]++;
                  if(FPGA_temp[2][3]>ED_xChan_maxNo)FPGA_temp[2][3]=0;

                  //for status_reg[3] and status_reg[4]
                      FPGA_add[3][3]=(onedata[299]&0xF0)>>4;
                  sprintf(hname, "FPGA_4_4");	//FPGA_4_4 is for Trigger status display
                  HF1cont(hname,FPGA_temp[3][3],FPGA_add[3][3]);
                  FPGA_temp[3][3]++;
                  if(FPGA_temp[3][3]>ED_xChan_maxNo)    FPGA_temp[3][3]=0;

                      FPGA_add[4][3]=FPGA4_HLAdd(onedata[299],onedata[300]);
                  sprintf(hname, "FPGA_5_4");	//FPGA_5_4 is for Trigger Number
                  HF1cont(hname,FPGA_temp[4][3],FPGA_add[4][3]);
                  FPGA_temp[4][3]++;
                  if(FPGA_temp[4][3]>ED_xChan_maxNo)    FPGA_temp[4][3]=0;

                  //for status_reg[5] and status_reg[6]
                      FPGA_add[5][3]=(onedata[301]&0xF0)>>4;
                  sprintf(hname, "FPGA_6_4");	//FPGA_6_4 is for times of current auto power-off
                  HF1cont(hname,FPGA_temp[5][3],FPGA_add[5][3]);
                  FPGA_temp[5][3]++;
                  if(FPGA_temp[5][3]>ED_xChan_maxNo)  FPGA_temp[5][3]=0;

                      FPGA_add[6][3]=FPGA4_HLAdd(onedata[301],onedata[302]);
                  sprintf(hname, "FPGA_7_4");	//FPGA_7_4 is for current detection value display
                  HF1cont(hname,FPGA_temp[6][3],FPGA_add[6][3]);
                  FPGA_temp[6][3]++;
                  if(FPGA_temp[6][3]>ED_xChan_maxNo)    FPGA_temp[6][3]=0;

                  //for status_reg[7]
                  FPGA_add[7][3]=(onedata[303]&0x0F);
                  sprintf(hname, "FPGA_8_4_1");	//FPGA_8_4_1 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[7][3],FPGA_add[7][3]);
                      FPGA_add[7][3]=(onedata[303]&0xF0)>>4;
                  sprintf(hname, "FPGA_8_4_2");	//FPGA_8_4_2 is for data management command injection
                  HF1cont(hname,FPGA_temp[7][3],FPGA_add[7][3]);
                  FPGA_temp[7][3]++;
                  if(FPGA_temp[7][3]>ED_xChan_maxNo)    FPGA_temp[7][3]=0;

                  //for status_reg[8] and status_reg[9]
                      FPGA_add[8][3]=FPGA8_HLAdd(onedata[304],onedata[305]);
                  sprintf(hname, "FPGA_9_4");	//FPGA_9_4 is for timer of data management command injection packages
                  HF1cont(hname,FPGA_temp[8][3],FPGA_add[8][3]);
                  FPGA_temp[8][3]++;
                  if(FPGA_temp[8][3]>ED_xChan_maxNo)   FPGA_temp[8][3]=0;

                  //for status_reg[10] and status_reg[11]
                      FPGA_add[10][3]=FPGA8_HLAdd(onedata[306],onedata[307]);
                  sprintf(hname, "FPGA_11_4");	//FPGA_11_4 is for timer of FEE command response
                  HF1cont(hname,FPGA_temp[10][3],FPGA_add[10][3]);
                  FPGA_temp[10][3]++;
                  if(FPGA_temp[10][3]>ED_xChan_maxNo)   FPGA_temp[10][3]=0;

                  //for status_reg[12] and status_reg[13]
                      FPGA_add[12][3]=FPGA4_HLAdd(onedata[308],onedata[309]);
                  sprintf(hname, "FPGA_13_4");	//FPGA_13_4 is for threshold of auto power off
                  HF1cont(hname,FPGA_temp[12][3],FPGA_add[12][3]);
                  FPGA_temp[12][3]++;
                  if(FPGA_temp[12][3]>ED_xChan_maxNo)   FPGA_temp[12][3]=0;

                  //for status_reg[14] and status_reg[15]
                      FPGA_add[14][3]=FPGA8_HLAdd(onedata[310],onedata[311]);
                  sprintf(hname, "FPGA_15_4");	//FPGA_15_4 is for timer of power-off
                  HF1cont(hname,FPGA_temp[14][3],FPGA_add[14][3]);
                  FPGA_temp[14][3]++;
                  if(FPGA_temp[14][3]>ED_xChan_maxNo)  FPGA_temp[14][3]=0;

                      //for status_reg[16]
                      FPGA_add[16][3]=(onedata[312]&0xF0)>>4;
                  sprintf(hname, "FPGA_17_4");	//FPGA_17_4 is for command and trigger channel display
                  HF1cont(hname,FPGA_temp[16][3],FPGA_add[16][3]);
                  FPGA_temp[16][3]++;
                  if(FPGA_temp[16][3]>ED_xChan_maxNo) FPGA_temp[16][3]=0;

                  //for status_reg[17] and status_reg[18]
                      FPGA_add[17][3]=FPGA8_HLAdd(onedata[313],onedata[314]);
                  sprintf(hname, "FPGA_18_4");	//FPGA_18_4 is for timer of odd erro
                  HF1cont(hname,FPGA_temp[17][3],FPGA_add[17][3]);
                  FPGA_temp[17][3]++;
                  if(FPGA_temp[17][3]>ED_xChan_maxNo)   FPGA_temp[17][3]=0;

                  //for status_reg[19] and status_reg[20]
                      FPGA_add[19][3]=FPGA8_HLAdd(onedata[315],onedata[316]);
                  sprintf(hname, "FPGA_20_4");	//FPGA_20_4 is for timer of accumative sum
                  HF1cont(hname,FPGA_temp[19][3],FPGA_add[19][3]);
                  FPGA_temp[19][3]++;
                  if(FPGA_temp[19][3]>ED_xChan_maxNo)   FPGA_temp[19][3]=0;

                  //for status_reg[21] and status_reg[22]
                      FPGA_add[21][3]=FPGA8_HLAdd(onedata[317],onedata[318]);
                  sprintf(hname, "FPGA_22_4");	//FPGA_22_4 is for timer of invalid data management injection command
                  HF1cont(hname,FPGA_temp[21][3],FPGA_add[21][3]);
                  FPGA_temp[21][3]++;
                  if(FPGA_temp[21][3]>ED_xChan_maxNo)   FPGA_temp[21][3]=0;

                  //for status_reg[23]
                      FPGA_add[23][3]=onedata[319]&0xFF;
                  sprintf(hname, "FPGA_24_4");	//FPGA_24_4 is for timer of high lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[23][3],FPGA_add[23][3]);
                  FPGA_temp[23][3]++;
                  if(FPGA_temp[23][3]>ED_xChan_maxNo)   FPGA_temp[23][3]=0;

                  //for status_reg[24]
                      FPGA_add[24][3]=onedata[320]&0xFF;
                  sprintf(hname, "FPGA_25_4");	//FPGA_24_4 is for timer of low lever threshold CRC erro
                  HF1cont(hname,FPGA_temp[24][3],FPGA_add[24][3]);
                  FPGA_temp[24][3]++;
                  if(FPGA_temp[24][3]>ED_xChan_maxNo)   FPGA_temp[24][3]=0;

                  //for status_reg[25] and status_reg[26]
                      FPGA_add[25][3]=FPGA8_HLAdd(onedata[321],onedata[322]);
                  sprintf(hname, "FPGA_26_4");	//FPGA_26_4 is for scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[25][3],FPGA_add[25][3]);
                  FPGA_temp[25][3]++;
                  if(FPGA_temp[25][3]>ED_xChan_maxNo)   FPGA_temp[25][3]=0;

                  //for status_reg[27] and status_reg[28]
                      FPGA_add[27][3]=FPGA8_HLAdd(onedata[323],onedata[324]);
                  sprintf(hname, "FPGA_28_4");	//FPGA_28_4 is for accumulative sum of scientific data transmission timeout
                  HF1cont(hname,FPGA_temp[27][3],FPGA_add[27][3]);
                  FPGA_temp[27][3]++;
                  if(FPGA_temp[27][3]>ED_xChan_maxNo)   FPGA_temp[27][3]=0;

                  //for status_reg[29] and status_reg[30]
                      FPGA_add[29][3]=FPGA8_HLAdd(onedata[325],onedata[326]);
                  sprintf(hname, "FPGA_30_4");	//FPGA_30_4 is for timer of trigger width erro
                  HF1cont(hname,FPGA_temp[29][3],FPGA_add[29][3]);
                  FPGA_temp[29][3]++;
                  if(FPGA_temp[29][3]>ED_xChan_maxNo)   FPGA_temp[29][3]=0;

                  FPGA_count_temp[3]++;
                  triggercounttotal3=FPGA_count_temp[3];

                                break;//7:end for telemetering of Status of FEE at Negative(-)
                            case 2570://8:for telemetering of HV and DV-DC power FEE at Negative(-)
                                break;//8:end for telemetering of HV and DV-DC power FEE at Negative(-)

                            case 2585://9:for telemetering of HV and DV-DC power FEE at Positive(+)
                                break;//9:end for telemetering of HV and DV-DC power FEE at Positive(+)

                            case 3083://10:for HV engineering data
                                break;//10:end for HV engineering data

                            case 3585://11:for data processor engineering data
                                break;//11:end for data processor engineering data

                            case 3602://12:for DM engineering data
                                break;//12:end for DM engineering data

                            case 3983://13:for stallite status
                                break;//13:end for stallite status

                            case 3996://14: for DAMPE temprature data
                                break; //14: end for DAMPE temprature data

                      }//end switch FILEID
        }
    }
        else break;//end synchronous code
        gSystem->ProcessEvents();
    }



   datafile->Write(0,TObject::kOverwrite);
   fClose = kTRUE;
   fclose(fp);

     switch(FILEID)
       {
            case 2067://1:for	scientific data
              break;//1:end for	scientific data

            case 2313://2:for telemetering of current of FEE at Positive(+)
          sprintf(outinfo,"   Decoded Successfully:\n");
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at PX(I).\n",triggercounttotal1);
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at PY(I).\n",triggercounttotal3);
          ShowText(outinfo);
                break;//2:end for telemetering of current of FEE at Positive(+)

            case 2330://3:for telemetering of temperature of FEE at Positive(+)
          sprintf(outinfo,"   Decoded Successfully:\n");
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at PX(T).\n",triggercounttotal1);
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at PY(T).\n",triggercounttotal3);
          ShowText(outinfo);
                break;//3:for telemetering of temperature of FEE at Positive(+)

            case 2348://4:for telemetering of Status of FEE at Positive(+)
          sprintf(outinfo,"   Decoded Successfully:\n");
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at PX(FPGA).\n",triggercounttotal1);
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at PY(FPGA).\n",triggercounttotal3);
          ShowText(outinfo);
                break;//4:end for telemetering of Status of FEE at Positive(+)

            case 2367://5:for telemetering of current of FEE at Negative(-)
         sprintf(outinfo,"   Decoded Successfully:\n");
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at NX(I).\n",triggercounttotal2);
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at NY(I).\n",triggercounttotal4);
          ShowText(outinfo);
                break;//5:end for telemetering of current of FEE at Negative(-)

            case 2383://6:for telemetering of temperature of FEE at Negative(-)
          sprintf(outinfo,"   Decoded Successfully:\n");
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at NX(T).\n",triggercounttotal2);
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at NY(T).\n",triggercounttotal4);
          ShowText(outinfo);
                break;//6:end for telemetering of temperature of FEE at Negative(-)

            case 2396://7:for telemetering of Status of FEE at Negative(-)
          sprintf(outinfo,"   Decoded Successfully:\n");
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at NX(FPGA).\n",triggercounttotal2);
          ShowText(outinfo);
          sprintf(outinfo,"   %i events are processed at NY(FPGA).\n",triggercounttotal4);
          ShowText(outinfo);
                break;//7:end for telemetering of Status of FEE at Negative(-)

            case 2570://8:for telemetering of HV and DV-DC power FEE at Negative(-)
                break;//8:end for telemetering of HV and DV-DC power FEE at Negative(-)

            case 2585://9:for telemetering of HV and DV-DC power FEE at Positive(+)
                break;//9:end for telemetering of HV and DV-DC power FEE at Positive(+)

            case 3083://10:for HV engineering data
                break;//10:end for HV engineering data

            case 3585://11:for data processor engineering data
                break;//11:end for data processor engineering data

            case 3602://12:for DM engineering data
                break;//12:end for DM engineering data

            case 3983://13:for stallite status
                break;//13:end for stallite status

            case 3996://14: for DAMPE temprature data
                break; //14: end for DAMPE temprature data
   }//end switch FILEID

   processing = false;
   fBonefile->SetEnabled(true);
   sprintf(outinfo,"Processing End!");
   ShowText(outinfo);

    delete[] onevent;
}

void GuiFrame::OnRawdataDecoding()
{
    if(fHBG_mode->IsEnabled()){
        if((fHBG_type->GetButton(kB_pedestal)->GetState() == kButtonDown))
        {
            fHBG_mode->SetState(kFALSE);
            rawdata_type=kNormal;
        }
        else if((fHBG_type->GetButton(kB_calibration)->GetState() == kButtonDown))
        {
            fHBG_mode->SetState(kFALSE);
            rawdata_type=kCalibration;
        }
    }
    else if((fHBG_type->GetButton(kB_mips)->GetState() == kButtonDown)){
        fHBG_mode->SetState(kTRUE);
        if(fHBG_mode->GetButton(kB_normal)->GetState() == kButtonDown){
            rawdata_type=kNormal;
        }
        else if(fHBG_mode->GetButton(kB_compressed)->GetState() == kButtonDown){
            rawdata_type=kCompressed;
        }
        else if(fHBG_mode->GetButton(kB_smaller)->GetState() == kButtonDown){
            rawdata_type=kSmaller;
        }
    }
    else if(fHBG_type->GetButton(kB_pedestal)->GetState() == kButtonDown){
        rawdata_type=kNormal;
    }
    else if(fHBG_type->GetButton(kB_calibration)->GetState() == kButtonDown){
        rawdata_type=kCalibration;
    }

    fFile_Sci_ref->SetState(kFALSE);
    fBopen_Sci_ref->SetState(kButtonDisabled);

    fGF_analyze->SetTitle("Decode");
    fTBN_analyzeSci->SetState(kButtonUp);
    fTBN_showResult->SetState(kButtonDisabled);


}

void GuiFrame::OnRootfileAnalyze()
{
    fHBG_mode->SetState(kFALSE);

    if((fHBG_type->GetButton(kB_pedestal)->GetState()==kButtonDown) ||
            (fHBG_type->GetButton(kB_mips)->GetState()==kButtonDown))
    {
        fFile_Sci_ref->SetState(kTRUE);
        fBopen_Sci_ref->SetState(kButtonUp);
    }
    else if((fHBG_type->GetButton(kB_calibration)->GetState()==kButtonDown))
    {
        fFile_Sci_ref->SetState(kFALSE);
        fBopen_Sci_ref->SetState(kButtonDisabled);
    }

    fGF_analyze->SetTitle("Analyze");
    fTBN_analyzeSci->SetState(kButtonDisabled);
    fTBN_showResult->SetState(kButtonUp);
}

void GuiFrame::OnProcessConfig(Int_t process_id)
{
    switch (process_id) {
    case kB_rawdatadecoding:
        OnRawdataDecoding();
        /*
        char outinfo[20];
        sprintf(outinfo,"%x",rawdata_type);
        ShowText(outinfo);
        */
        break;
    case kB_rootfileanalyze:
        OnRootfileAnalyze();
        break;
    default:
        break;
    }
}

void GuiFrame::OnTypeConfig(Int_t type_id)
{
    analyze_type=type_id;
    if(fHBG_process->GetButton(kB_rawdatadecoding)->GetState() == kButtonDown){
        switch (type_id) {
        case kB_mips:
            if(!fHBG_mode->IsEnabled()){
                fHBG_mode->SetState(kTRUE);
            }
            if(fHBG_mode->GetButton(kB_normal)->GetState() == kButtonDown){
                rawdata_type=kNormal;
            }
            else if(fHBG_mode->GetButton(kB_compressed)->GetState() == kButtonDown){
                rawdata_type=kCompressed;
            }
            else if(fHBG_mode->GetButton(kB_smaller)->GetState() == kButtonDown){
                rawdata_type=kSmaller;
            }
            break;
        case kB_pedestal:
            fHBG_mode->SetState(kFALSE);
            rawdata_type=kNormal;
            break;
        case kB_calibration:
            fHBG_mode->SetState(kFALSE);
            rawdata_type=kCalibration;
            break;
        }

    }
    else if(fHBG_process->GetButton(kB_rootfileanalyze)->GetState() == kButtonDown)
    {
        switch(type_id){
        case kB_calibration:
            fFile_Sci_ref->SetState(kFALSE);
            fBopen_Sci_ref->SetState(kButtonDisabled);
            //rawdata_type=kCalibration;
            break;
        case kB_pedestal:
            fFile_Sci_ref->SetState(kTRUE);
            fBopen_Sci_ref->SetState(kButtonUp);
            //rawdata_type=kNormal;
            break;
        case kB_mips:
            fFile_Sci_ref->SetState(kTRUE);
            fBopen_Sci_ref->SetState(kButtonUp);
        }
    }
    /*
    char outinfo[20];
    sprintf(outinfo,"%x",rawdata_type);
    ShowText(outinfo);
    */
}

void GuiFrame::OnModeConfig(Int_t mode_id)
{

    switch (mode_id) {
    case kB_normal:
        rawdata_type=kNormal;
        break;
    case kB_compressed:
        rawdata_type=kCompressed;
        break;
    case kB_smaller:
        rawdata_type=kSmaller;
        break;
    default:
        break;
    }
    /*
    char outinfo[20];
    sprintf(outinfo,"%x",rawdata_type);
    ShowText(outinfo);
    */
}

void GuiFrame::OnUseTimecode(Bool_t flag)
{
    if(flag){
        fNE_year->SetState(kTRUE);
        fNE_month->SetState(kTRUE);
        fNE_day->SetState(kTRUE);
        fNE_hour->SetState(kTRUE);
        fNE_minute->SetState(kTRUE);
        fNE_second->SetState(kTRUE);
        fTBN_starttime->SetState(kButtonUp);
        fTBN_stoptime->SetState(kButtonUp);
    }
    else{
        fNE_year->SetState(kFALSE);
        fNE_month->SetState(kFALSE);
        fNE_day->SetState(kFALSE);
        fNE_hour->SetState(kFALSE);
        fNE_minute->SetState(kFALSE);
        fNE_second->SetState(kFALSE);
        fTBN_starttime->SetState(kButtonDisabled);
        fTBN_stoptime->SetState(kButtonDisabled);
    }
}

void GuiFrame::OnGetStartTime()
{
    Int_t year,month,day,hour,minute,second;
    year=fNE_year->GetIntNumber();
    month=fNE_month->GetIntNumber();
    day=fNE_day->GetIntNumber();
    hour=fNE_hour->GetIntNumber();
    minute=fNE_minute->GetIntNumber();
    second=fNE_second->GetIntNumber();

    starttime.Form("%d",year);

    if(month<10){
        starttime.Append(Form("0%d",month));
    }
    else{
        starttime.Append(Form("%d",month));
    }

    if(day<10){
        starttime.Append(Form("0%d",day));
    }
    else{
        starttime.Append(Form("%d",day));
    }

    if(hour<10){
        starttime.Append(Form("0%d",hour));
    }
    else{
        starttime.Append(Form("%d",hour));
    }

    if(minute<10){
        starttime.Append(Form("0%d",minute));
    }
    else{
        starttime.Append(Form("%d",minute));
    }

    if(second<10){
        starttime.Append(Form("0%d",second));
    }
    else{
        starttime.Append(Form("%d",second));
    }

    ShowText("Start Time:");
    ShowText(starttime.Data());
}

void GuiFrame::OnGetStopTime()
{
    Int_t year,month,day,hour,minute,second;
    year=fNE_year->GetIntNumber();
    month=fNE_month->GetIntNumber();
    day=fNE_day->GetIntNumber();
    hour=fNE_hour->GetIntNumber();
    minute=fNE_minute->GetIntNumber();
    second=fNE_second->GetIntNumber();

    stoptime.Form("%d",year);

    if(month<10){
        stoptime.Append(Form("0%d",month));
    }
    else{
        stoptime.Append(Form("%d",month));
    }

    if(day<10){
        stoptime.Append(Form("0%d",day));
    }
    else{
        stoptime.Append(Form("%d",day));
    }

    if(hour<10){
        stoptime.Append(Form("0%d",hour));
    }
    else{
        stoptime.Append(Form("%d",hour));
    }

    if(minute<10){
        stoptime.Append(Form("0%d",minute));
    }
    else{
        stoptime.Append(Form("%d",minute));
    }

    if(second<10){
        stoptime.Append(Form("0%d",second));
    }
    else{
        stoptime.Append(Form("%d",second));
    }

    ShowText("Stop time:");
    ShowText(stoptime.Data());
}

void GuiFrame::OnInputSciData()
{
    new TGFileDialog(gClient->GetRoot(),this,kFDOpen,fi_Sci_input);
    input_sci_filename=fi_Sci_input->fFilename;
    fFile_Sci_input->Clear();
    fFile_Sci_input->SetText(input_sci_filename.Data());

    ShowText("Input:");
    ShowText(input_sci_filename.Data());
}

void GuiFrame::OnOutputSciData()
{
    new TGFileDialog(gClient->GetRoot(),this,kFDSave,fi_Sci_output);
    output_sci_filename=fi_Sci_output->fFilename;
    fFile_Sci_output->Clear();
    fFile_Sci_output->SetText(output_sci_filename.Data());

    ShowText("Output:");
    ShowText(output_sci_filename.Data());
}

void GuiFrame::OnRefSciData()
{
    TGFileDialog *dialog;
    dialog=new TGFileDialog(gClient->GetRoot(),this,kFDOpen,fi_Sci_ref);
    ref_sci_filename=fi_Sci_ref->fFilename;
    fFile_Sci_ref->Clear();
    fFile_Sci_ref->SetText(ref_sci_filename.Data());

    ShowText("Reference:");
    ShowText(ref_sci_filename.Data());

    //printf("%d\n",dialog);
    //delete dialog;
}

void GuiFrame::OnSciAnalyze()
{
    ShowText("Analyzing Begin.");
    TString calib_configfile;
    //---input dir---
    TString inputDirName,inputBaseName;
    inputBaseName=gSystem->BaseName(input_sci_filename.Data());
    inputDirName=input_sci_filename;
    inputDirName.ReplaceAll(inputBaseName,"");
    //---output dir--
    TString outputDirName,outputBaseName;
    outputBaseName=gSystem->BaseName(output_sci_filename.Data());
    outputDirName=output_sci_filename;
    outputDirName.ReplaceAll(outputBaseName,"");
    //---ref dir---
    TString refDirName,refBaseName;
    TString prefix;
    TString newped_filename;
    TString tempstr;

    //---analyze--
    FILE* fp;
    switch (analyze_type) {
    case kB_pedestal:
        //ShowText("pedestal");
        refBaseName=gSystem->BaseName(ref_sci_filename.Data());
        refDirName=ref_sci_filename;
        refDirName.ReplaceAll(refBaseName,"");
        fp=fopen(Form("%s/pedestal.log",outputDirName.Data()),"w");
        fprintf(fp,"Pedestal Analyze Log Info:\n");
        fprintf(fp,"\tRaw root file dir: %s\n",inputDirName.Data());
        fprintf(fp,"\tRaw root file name: %s\n",inputBaseName.Data());
        fprintf(fp,"\tReference pedestal root file dir: %s\n",refDirName.Data());
        fprintf(fp,"\tReference pedestal root file name: %s\n",refBaseName.Data());
        fprintf(fp,"Analyze Result:\n");
        fprintf(fp,"\tOutput file dir: %s\n",outputDirName.Data());
        fprintf(fp,"\tOutput Log file(this file): pedestal.log\n");
        fprintf(fp,"\tOutput root file: %s\n",outputBaseName.Data());
        prefix=outputBaseName;
        prefix.ReplaceAll(".root","");
        getpedseed_event(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),prefix.Data());
        prefix.Append("_analyze.root");
        draw_rel(output_sci_filename.Data(),ref_sci_filename.Data(),outputDirName.Data(),prefix.Data());
        draw_relp(output_sci_filename.Data(),ref_sci_filename.Data(),outputDirName.Data(),prefix.Data());
        draw_ped(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),prefix.Data());
        break;
    case kB_calibration:
        //ShowText("calibration");
        fp=fopen(Form("%s/calibration.log",outputDirName.Data()),"w");
        fprintf(fp,"Calibration Analyze Log Info:\n");
        fprintf(fp,"Working dir: %s\n",bin_dir.Data());
        fprintf(fp,"\tRaw root file dir: %s\n",inputDirName.Data());
        fprintf(fp,"\tRaw root file name: %s\n",inputBaseName.Data());
        fprintf(fp,"Analyze Result:\n");
        fprintf(fp,"\tOutput file dir: %s\n",outputDirName.Data());
        fprintf(fp,"\tOutput Log file(this file): calibration.log\n");
        fprintf(fp,"\tOutput root file: %s\n",outputBaseName.Data());
        outputBaseName.ReplaceAll(".root","");
        calib_configfile= bin_dir+"/calib.config";
        fit_calibration(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),outputBaseName.Data(),calib_configfile.Data());
        break;
    case kB_mips:
        //ShowText("mips");
        fp=fopen(Form("%s/mips.log",outputDirName.Data()),"w");
        refBaseName=gSystem->BaseName(ref_sci_filename.Data());
        refDirName=ref_sci_filename;
        refDirName.ReplaceAll(refBaseName,"");
        fprintf(fp,"MIPs Analyze Log Info:\n");
        fprintf(fp,"\tRaw root file dir: %s\n",inputDirName.Data());
        fprintf(fp,"\tRaw root file name: %s\n",inputBaseName.Data());
        fprintf(fp,"\tReference pedestal root file dir: %s\n",refDirName.Data());
        fprintf(fp,"\tReference pedestal root file name: %s\n",refBaseName.Data());
        fprintf(fp,"Analyze Result:\n");
        fprintf(fp,"\tOutput file dir: %s\n",outputDirName.Data());
        fprintf(fp,"\tOutput Log file(this file): mips.log\n");
        fprintf(fp,"\tOutput root file: %s\n",outputBaseName.Data());
        draw_channels(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),outputBaseName.Data());
        newped_filename= outputDirName;
        newped_filename.Append("/new_ped.root");
        draw_relp(newped_filename.Data(),ref_sci_filename.Data(),outputDirName.Data(),outputBaseName.Data());
        draw_mip(input_sci_filename.Data(),newped_filename.Data(),outputDirName.Data(),outputBaseName.Data());
        draw_mapping(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data());
        fit_dy58(input_sci_filename.Data(),newped_filename.Data(),outputDirName.Data(),"dy58");
        draw_pedVStime(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),outputBaseName.Data());
        extract_psd(input_sci_filename.Data(),newped_filename.Data(),outputDirName.Data(),"analysis.root");
        tempstr=outputDirName+"analysis.root";
        draw_hitnum(tempstr.Data(),"hit",5);
        bt_draw_meangeo(tempstr.Data());
        bt_efficiency(tempstr.Data(),25,25,5);
        break;
    default:
        break;
    }
    fclose(fp);
    ShowText("Analyzing End.");
}

void GuiFrame::OnSciDecode()
{
    //---input dir---
    TString inputDirName,inputBaseName;
    inputBaseName=gSystem->BaseName(input_sci_filename.Data());
    inputDirName=input_sci_filename;
    inputDirName.ReplaceAll(inputBaseName,"");
    //---output dir--
    TString outputDirName,outputBaseName;
    outputBaseName=gSystem->BaseName(output_sci_filename.Data());
    outputDirName=output_sci_filename;
    outputDirName.ReplaceAll(outputBaseName,"");
    //---convert----
    ShowText("Start Decoding:");
    
    FILE* fp=fopen(Form("%s/decode.log",outputDirName.Data()),"w");
    
    convert_psd_scidata(fp,rawdata_type,inputDirName.Data(),inputBaseName.Data(),
                        outputDirName.Data(),outputBaseName.Data());
    //convert_event(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),outputBaseName.Data());
    
    fclose(fp);
    /*
    switch (rawdata_type) {
    case kNormal:
        //convert_event(inputDirName.Data(),inputBaseName.Data(),outputDirName.Data(),outputBaseName.Data());
        convert_psd_scidata(fp,rawdata_type,inputDirName.Data(),inputBaseName.Data(),
                        outputDirName.Data(),outputBaseName.Data());
        break;
    case kCompressed:
        convert_psd_scidata(fp,rawdata_type,inputDirName.Data(),inputBaseName.Data(),
                        outputDirName.Data(),outputBaseName.Data());
        //ShowText("kCompressed");
        break;
    case kSmaller:
        convert_psd_scidata(fp,rawdata_type,inputDirName.Data(),inputBaseName.Data(),
                        outputDirName.Data(),outputBaseName.Data());
        break;
    default:
        break;
    }
*/
    ShowText("Decode End.");
    //ShowText(inputDirName.Data());
    //ShowText(inputBaseName.Data());
    //ShowText(outputDirName.Data());
    //ShowText(outputBaseName.Data());
    //ShowText(starttime.Data());
    //ShowText(Form("%d",convert_dateTotimecode(starttime.Data())));
}

void GuiFrame::OnInputLowthresh()
{
    //ShowText("what");
    //new TGFileDialog(gClient->GetRoot());
    new TGFileDialog(gClient->GetRoot(),this,kFDOpen,fi_lowthreshOpen);
    input_lowthresh_filename=fi_lowthreshOpen->fFilename;
    fFile_lowthreshOpen->Clear();
    fFile_lowthreshOpen->SetText(input_lowthresh_filename.Data());

    ShowText("Low Thresh reference file:");
    ShowText(input_lowthresh_filename.Data());
}

void GuiFrame::OnOutputLowthresh()
{
    //ShowText("what");
    new TGFileDialog(gClient->GetRoot(),this,kFDSave,fi_lowthreshSave);
    output_lowthresh_filename=fi_lowthreshSave->fFilename;
    fFile_lowthreshSave->Clear();
    fFile_lowthreshSave->SetText(output_lowthresh_filename.Data());

    ShowText("Low Thresh output:");
    ShowText(output_lowthresh_filename.Data());

}

void GuiFrame::OnOutputHighthresh()
{
    //ShowText("what");
    new TGFileDialog(gClient->GetRoot(),this,kFDSave,fi_highthreshSave);
    output_highthresh_filename=fi_highthreshSave->fFilename;
    fFile_highthreshSave->Clear();
    fFile_highthreshSave->SetText(output_highthresh_filename.Data());

    ShowText("High Thresh output:");
    ShowText(output_highthresh_filename.Data());
}

void GuiFrame::OnLowthreshEncoding()
{
    Int_t level=fNE_lowthresh->GetIntNumber();
    TString outputDirName,outputBaseName;
    outputBaseName=gSystem->BaseName(output_lowthresh_filename.Data());
    outputDirName=output_lowthresh_filename;
    outputDirName.ReplaceAll(outputBaseName,"");

    FILE* fp=fopen(output_lowthresh_filename.Data(),"w");
    fprintf(fp,"Low Threshold Configuraion:\n");
    fprintf(fp,"\t1) threshold: mean + %d*sigma\n",level);
    fprintf(fp,"\t2) reference file: %s\n",input_lowthresh_filename.Data());
    fprintf(fp,"\t3) output directory: %s\n",outputDirName.Data());
    fprintf(fp,"\t4) output binary threshold file:\n");
    fprintf(fp,"\t\txpos_lowthresh.bin\n");
    fprintf(fp,"\t\txneg_lowthresh.bin\n");
    fprintf(fp,"\t\typos_lowthresh.bin\n");
    fprintf(fp,"\t\tyneg_lowthresh.bin\n");
    fprintf(fp,"\t5) log file(this file): %s\n",outputBaseName.Data());

    getconfig_lowthresh(input_lowthresh_filename.Data(),outputDirName.Data(),level);
    fclose(fp);
    /*
    TString basename,dirname;
    basename=gSystem->BaseName(input_lowthresh_filename.Data());
    dirname=input_lowthresh_filename;
    dirname.ReplaceAll(basename,"");

    char buffer[200];
    sprintf(buffer,"fullname:%s",input_lowthresh_filename.Data());
    ShowText(buffer);
    sprintf(buffer,"dirname:%s",dirname.Data());
    ShowText(buffer);
    sprintf(buffer,"basename:%s",basename.Data());
    ShowText(buffer);
    */

}

void GuiFrame::OnHighthreshEncoding()
{
    Int_t value=fNE_highthresh->GetHexNumber();
    char highthresh=value&0xFF;
    TString outputDirName,outputBaseName;
    outputBaseName=gSystem->BaseName(output_highthresh_filename.Data());
    outputDirName=output_highthresh_filename;
    outputDirName.ReplaceAll(outputBaseName,"");

    FILE* fp=fopen(output_highthresh_filename.Data(),"w");
    fprintf(fp,"High Threshold Configuraion:\n");
    fprintf(fp,"\t1) threshold: %d",(unsigned char)highthresh);
    fprintf(fp,"\t2) output directory: %s\n",outputDirName.Data());
    fprintf(fp,"\t3) output binary threshold file:\n");
    fprintf(fp,"\t\txpos_highthresh.bin\n");
    fprintf(fp,"\t\txneg_highthresh.bin\n");
    fprintf(fp,"\t\typos_highthresh.bin\n");
    fprintf(fp,"\t\tyneg_highthresh.bin\n");
    fprintf(fp,"\t4) log file(this file): %s\n",outputBaseName.Data());

    getconfig_highthresh(outputDirName.Data(),highthresh);
    fclose(fp);
    /*
    char buffer[20];
    sprintf(buffer,"%d",highthresh);
    ShowText(buffer);
    */
}

void GuiFrame::OnDrawChannel()
{

}

void GuiFrame::OnScanChannels()
{

}

void GuiFrame::SetResultState(bool flag)
{
    if(flag){

        fTBN_DrawChannel->SetState(kButtonUp);
        fTBN_ScanChannels->SetState(kButtonUp);

        fCB_strip->SetEnabled(kTRUE);
        fCB_layer->SetEnabled(kTRUE);
        fCB_side->SetEnabled(kTRUE);
        fCB_dynode->SetEnabled(kTRUE);
    }
    else{
        fTBN_DrawChannel->SetState(kButtonDisabled);
        fTBN_ScanChannels->SetState(kButtonDisabled);

        fCB_strip->SetEnabled(kFALSE);
        fCB_layer->SetEnabled(kFALSE);
        fCB_side->SetEnabled(kFALSE);
        fCB_dynode->SetEnabled(kFALSE);
    }
}
