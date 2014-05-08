/////////////////////////////////////////////////////////
// File name: GuiFrame.h                            //
// Brief introduction:                                 //
//       This class create the main frame for          //
//       gui program of Process Files                  //
//                                                     //
// Version: V1.1                                       //
// Author: liu longxiang                               //
// Date: Nov. 2013                                     //
// Mail: liulongxiang@impcas.ac.cn                     //
/////////////////////////////////////////////////////////

#ifndef GuiFrame_H
#define GuiFrame_H

#include "TGFrame.h"
#include "TString.h"

class TGMainFrame;
class TGGroupFrame;
class TGVerticalFrame;
class TGWindow;
class TGFileContainer;
class TGFileInfo;
class TGLabel;
class TGListView;
class TGLVEntry;
class TGTextView;
class TGComboBox;
class TGTextButton;
class TGText;
class TCanvas;
class TString;
class TList;
class TIter;
class TObject;
class TGTextEntry;
class TGNumberEntry;
class TGHProgressBar;
class TGTextBuffer;
class TFile;
class TTree;
class TGCheckButton;
class TGHorizontalFrame;
class TRootEmbeddedCanvas;
class TGButtonGroup;
class TGHButtonGroup;


enum CMDIdentifiers
{
    kB_process,
    kB_stop,
    kB_exit,
    kB_resetcurr,
    kB_resetall,
    kB_privious,
    kB_next,
    kB_update,
    kB_integral,
    kB_openEngineer,
    kB_inputSciData,
    kB_outputSciData,
    kB_refSciData,
    kB_analyzeSci,
    kB_showResult,
    kB_drawChannel,
    kB_scanChannel,
    kB_lowthreshEncoding,
    kB_highthreshEncoding,
    kB_lowthreshOpen,
    kB_lowthreshSave,
    kB_highthreshSave,
    kB_getStarttime,
    kB_getStoptime
};

enum ProcessID{
    kB_rawdatadecoding,
    kB_rootfileanalyze
};

enum TypeID
{
    kB_pedestal,
    kB_calibration,
    kB_mips
};

enum ModeID
{
    kB_normal,
    kB_compressed,
    kB_smaller
};

enum NumEntryID
{
    kNE_year,
    kNE_month,
    kNE_day,
    kNE_hour,
    kNE_minute,
    kNE_second,
    kNE_lowthresh,
    kNE_highthresh
};

enum DataType
{
    kNormal=0x20,
    kCompressed=0x60,
    kCalibration=0xa0,
    kSmaller=0x30
};
class GuiFrame : public TGMainFrame
{
    ClassDef(GuiFrame,0)
protected:
    //TGVerticalFrame *fFFileList;
    TGGroupFrame *fFFileList;
    TGVerticalFrame *fFDrawPanel;
    TGGroupFrame *fFScidataPanel;
    TGHorizontalFrame *fFTextView;
    TGLabel  *fLFilename_Engineer;
    TGFileInfo *fi_engineer;
    TGTextEntry  *fFileEngineer;
    TGTextBuffer  *fTbEngineer;
    TGHProgressBar    *fHProgH;

    TGLabel *flabDrawOpt, *flabText, *flabDivPad;
    TGComboBox *fcomDrawOpt;
    TGListView  *flvFile;
    TGFileContainer  *fFileCont;
    TGTextButton *fBonefile,*fBstop, *fBresetcurr, *fBresetall, *fBprevious, *fBnext, *fBupdate, *fBintegral;
    TGTextButton *fBdisplayfile;
    TGComboBox *fcomDivPad;
    TGTextView *fviewText;
    TRootEmbeddedCanvas *fEmbeddedCanvas;

    TGHButtonGroup *fHBG_process;
    TGHButtonGroup *fHBG_type;
    TGHButtonGroup *fHBG_mode;
    TGLabel    *fLFilename_Sci_input;
    TGLabel    *fLFilename_Sci_output;
    TGLabel    *fLFilename_Sci_ref;
    TGTextBuffer  *fTb_Sci_input;
    TGTextBuffer  *fTb_Sci_output;
    TGTextBuffer  *fTb_Sci_ref;
    TGTextButton* fBopen_Sci_input;
    TGTextButton* fBopen_Sci_output;
    TGTextButton* fBopen_Sci_ref;
    TGTextEntry*  fFile_Sci_input;
    TGTextEntry*  fFile_Sci_output;
    TGTextEntry*  fFile_Sci_ref;
    TGFileInfo *fi_Sci_input;
    TGFileInfo *fi_Sci_output;
    TGFileInfo *fi_Sci_ref;

    TGNumberEntry *fNE_year;
    TGNumberEntry *fNE_month;
    TGNumberEntry *fNE_day;
    TGNumberEntry *fNE_hour;
    TGNumberEntry *fNE_minute;
    TGNumberEntry *fNE_second;
    TGTextButton *fTBN_starttime;
    TGTextButton *fTBN_stoptime;

    TGNumberEntry *fNE_lowthresh;
    TGNumberEntry *fNE_highthresh;
    TGTextButton *fTBN_lowthreshEncoding;
    TGTextButton *fTBN_highthreshEncoding;
    TGTextButton *fTBN_lowthreshOpen;
    TGTextEntry  *fFile_lowthreshOpen;
    TGTextBuffer *fTb_lowthreshOpen;
    TGFileInfo   *fi_lowthreshOpen;
    TGTextButton *fTBN_lowthreshSave;
    TGTextEntry  *fFile_lowthreshSave;
    TGTextBuffer *fTb_lowthreshSave;
    TGFileInfo   *fi_lowthreshSave;
    TGTextButton *fTBN_highthreshSave;
    TGTextEntry  *fFile_highthreshSave;
    TGTextBuffer *fTb_highthreshSave;
    TGFileInfo   *fi_highthreshSave;

    TGGroupFrame *fGF_analyze;
    TGTextButton *fTBN_analyzeSci;
    TGTextButton *fTBN_showResult;
    TGTextButton *fTBN_DrawChannel;
    TGTextButton *fTBN_ScanChannels;
    TGComboBox *fCB_layer;
    TGComboBox *fCB_strip;
    TGComboBox *fCB_side;
    TGComboBox *fCB_dynode;

    TCanvas *canpf;
    TObject *objcurr;
    TList *ObjList;
    TFile *datafile;
    TTree *tree;
    TList *TCmlist;   //TCavas Contextmenu list
    TList *THmlist;   //TH1F Contextmenu list
    TList *TH2mlist;  //TH2F Contextmenu list

    bool processing;
    unsigned short *onevent;
    Bool_t          fClose;

public:
    GuiFrame(const TGWindow *p, UInt_t w, UInt_t h);
    virtual ~GuiFrame();
    virtual void CloseWindow();
    //virtual bool ProcessMessage(Long_t msg, Long_t param1, Long_t);
    const char* GetDrawOpt();
    void ShowText(TGText *text);
    void ShowText(const char *text);
    virtual void CreateCanvas();
    void ClearTextView();
    virtual void OnDoubleClick(TGLVEntry* f, Int_t btn);
    virtual void DisplayFile(const TString &fname);
    virtual void DisplayObject(const TString& fname,const TString& name);
    virtual void DisplayDirectory(const TString &fname);
    virtual void ImB_update();
    virtual void ImB_next();
    virtual void ImB_previous();
    virtual void ImB_integral();
    virtual bool ObjListOK();
    virtual void DrawObj(TObject *obj);
    virtual void DrawObj(CMDIdentifiers id);
    virtual void ResetTH(TObject *obj);
    virtual void ImB_ResetAllTH();
    virtual void ImB_ResetCurrTH();
    virtual void GetDivPad(int &nx, int &ny);
    virtual void RemoveMenuEntry(const char *menuTitle, TList *mlist);
    virtual void MakeTcMenuList();
    virtual void MakeTH1MenuList();
    virtual void MakeTH2MenuList();
    friend void HBOOK1(const  char *hname, const char *title, int nxbin, float xlow, float xup, float vmx);
    friend void HBOOK2(const char *hname, const char *title, int nxbin, float xlow, float xup, int nybin, float ylow, float yup, float vmx);
    friend void HF1(const char *hname, float value, float weight);
    friend void HF1cont(const char *hname, int nbin, float value);
    friend void HF2(const char *hname, float x, float y, float weight);

 //   virtual bool ProcessMessage(Long_t msg, Long_t param1, Long_t);
    void Stop();
    void DecodeEngineerData(TString filename);
    void ProcessFiles();
    void OpenEngineerData();
    void OnInputSciData();
    void OnOutputSciData();
    void OnRefSciData();
    void OnRawdataDecoding();
    void OnRootfileAnalyze();
    void OnProcessConfig(Int_t process_id);
    void OnTypeConfig(Int_t type_id);
    void OnModeConfig(Int_t mode_id);
    void OnSciDecode();
    void OnSciAnalyze();
    void OnInputLowthresh();
    void OnOutputLowthresh();
    void OnOutputHighthresh();
    void OnLowthreshEncoding();
    void OnHighthreshEncoding();
    void OnDrawChannel();
    void OnScanChannels();
    void SetResultState(bool flag);
    void OnGetStartTime();
    void OnGetStopTime();
protected:
    GuiFrame(const GuiFrame &onf);
    GuiFrame& operator=(const GuiFrame &onf);

    TString input_sci_filename;
    TString output_sci_filename;
    TString ref_sci_filename;
    Int_t year,month,day,hour,minute,second;
    TString input_lowthresh_filename;
    TString output_lowthresh_filename;
    TString output_highthresh_filename;
    Int_t rawdata_type;
    Int_t analyze_type;
    TString starttime;
    TString stoptime;
};

#endif //#ifndef GuiFrame_H
