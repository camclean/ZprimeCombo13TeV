#ifndef STEERPARSER_H
#define STEERPARSER_H

#include <TObject.h>
#include <TObjArray.h>
#include <TString.h>

class BaseSteer;
class TClass;
class TObjLink;

class SteerParser: public TObject {

public:
    SteerParser();
    virtual ~SteerParser();

    Bool_t   CreateSteering(TString& theclassname, TString& objectname, TString& block, 
			    const char* nameSpace = NULL);
    Int_t    GetNumSteerings(TClass* theclass);
    Int_t    GetNumSteerings() { return fSteerings.GetEntriesFast(); }
    BaseSteer* GetSteer(TClass* cl);
    BaseSteer* GetSteer(TClass* cl, const char* name);
    Bool_t   ParseFile(const char* filename);
    virtual void     Print(Option_t* option="") const;
    virtual void     Print(Option_t* option, Int_t first, Int_t last) const;
    TString  ReadString(const char* keyword);
    Float_t  ReadFloat(const char* keyword);
    Int_t    ReadInt(const char* keyword);
    Int_t    ReadArray(const char* keyword,Float_t array[]);
    Int_t    ReadArray(const char* keyword,Int_t array[]);
    Int_t    ReadArray(const char* keyword, TObjArray& array, Bool_t checkNamespace = kFALSE);
    void     ReadFastArray(const char* keyword,Float_t array[]);
    void     ReadFastArray(const char* keyword,Int_t array[]);
    void     SetFilename(const char* filename="stdin");
    void     SetClassname(const char* classname, const char* objectname="") ;
    Int_t     NextNameSpace();
    const char* GetCurrentNameSpace() const;
    Bool_t    SetCurrentNameSpace(const char* name);
    
private:
    void ParseLine(TString& line);
    Bool_t FindKeywordInLine(TString& line, const TString& keyword);
    Bool_t FindClassnameInLine(TString& line) const;
    Bool_t LineEnd(const TString& line) const;
    TString ValueOfLine(const TString& line) const;
    void RemoveComments(TString& line) const;
    void RemoveWhiteSpace(TString& line) const;
    void CountBraces(TString& line);
    Bool_t   FindArrayInFile(const char* keyword, Bool_t checkNamespace = kFALSE);
    void    SelectNamespace(TList& steerLines, const char* nameSpace) const;
    TObjLink* RemoveLink(TList& list, TObjLink* link) const;
    TString   fFilename;   // current filename
    TString   fNamespace;  // current namespace
    TString   fClassname;  // current classname
    TString   fObjectname; // current name of the steering object
    TString   fBlock;      // current block
    TString   fLine;       // current line
    Bool_t    fInNamespace;// flag
    Bool_t    fInBlock;    // flag

    TObjArray fSteerings;  // list of steering objects
    Int_t     fBraces;     // current  total of '{' minus '}'   
    Int_t     fArrayLength;// current length of array defined in steering

    TObjArray fNameSpaces; // storing 'namespaces' if any required
    Int_t     fCurrentNameIndex;// no of current 'namespace', -1 if none required

};

#endif
