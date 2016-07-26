#ifndef BASESTEER_H
#define BASESTEER_H

#include <TNamed.h>

class TString;
class TObjArray;
class TArrayI;
class TArrayD;
class TArrayF;

class BaseSteer: public TNamed {
 public:
    virtual void  SetValues(TString& block);
    virtual Int_t SplitExpression(const TString& exp, TString& left,
                                  TString& right) const;
    virtual void  SplitString(const char* str, const char* splitStr,
			      TObjArray* result) const;
    virtual Int_t StringToArray(const char* str, TArrayI& result) const;
    virtual Int_t StringToArray(const char* str, TArrayD& result) const;
    virtual Int_t StringToArray(const char* str, TArrayF& result) const;

    virtual void PrintArray(const char* description, const TArrayI& result) const;
    virtual void PrintArray(const char* description, const TArrayD& result) const;
    virtual void PrintArray(const char* description, const TArrayF& result) const;

 private:
    ClassDef(BaseSteer,0)
};

#endif
