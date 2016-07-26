/********************************************************************
 * 
 * @Package: Plotter
 *
 *
 * Class for setting values from a steering file.
 * The parsing of the steering file is done via SteerParser
 * that delivers blocks of a file to the concrete subclasses
 *
 * Derived classes have to provide
 *
 * - datamembers
 * - getter and setter functions for those datamembers (setters must NOT be inline!)
 * - the ROOT dictionary (ClassDef/ClassImp)
 * - Print(Option_t* option="") const
 *
 * Everything else is handled by this class via the ROOT RTTI features.
 *
 * Steerings can have a name (GetName()) and can belong to a namespace (GetTitle()).
 * This enables quite flexible ways of 'multiple' steering.
 *
 * Example of a derived steering class and its part in the steering file:
 * 
 *
 * in MySteer.h
 *==============================================
 * class MySteer: public BaseSteer {
 * public:
 *   MySteer();
 *   virtual ~MySteer();
 *
 *   void Print(Option_t* option="") const;
 *
 *   Float_t GetMyFloat();
 *   TArrayI GetArrayI();
 *   const char* GetOneString();
 *   TObjArray* GetStrings();
 * 
 *   void SetMyFloat(Float_t value);
 *   void SetArrayI(const char* arrayString);
 *   void SetOneString(const char* string);
 *   void SetStrings(const char* strings);
 *
 * private:
 *   Float_t   fMyFloat;
 *   TArrayI   fArrayI;
 *   TString   fOneString;
 *   TObjArray fStrings;
 *
 *   ClassDef(MySteer,0)  * steering class for the analysis
 * };
 *
 * and in MySteer.C:
 * ===========================================================
 * 
 * ClassImp(MySteer)
 *
 * MySteer::MySteer()
 * {
 *  // Here you have to set the default values which will be used,
 *  // if this steering class is required but not given in the steering file
 * }
 *
 * MySteer::~MySteer()
 * {
 *   fStrings.Delete();
 * }
 * void MySteer::Print(Option_t* option) const
 * {
 * // Prints all settings of the steering
 *   cout << "--MySteer values:\n";
 *       // Whatever you want to print to screen after reading in 
 *       // the values from file.
 *   cout << endl;
 * }
 *
 * void MySteer::SetMyFloat(Float_t value)
 * {
 *   fMyFloat = value;
 * }
 *
 * Float_t MySteer::GetMyFloat()  {return fMyFloat;}
 *
 * void MySteer::SetArrayI(const char* in)
 * {
 *   this->StringToArray(in, fArrayI);
 * }
 *
 * TArrayI MySteer::GetArrayI()  {return fArrayI;}
 *
 * void MySteer::SetOneString(const char* in) {fOneString = in;}
 *
 * const char* MySteer::GetOneString()  {return fOneString.Data();}
 *
 * void MySteer::SetStrings(const char* in)
 * {
 *   this->SplitString(in,",",&fStrings);
 * }
 *
 * // This method will return a TObjArray containing TObjString's.
 * // TObjString::GetName() returns the string
 * TObjArray* MySteer::GetStrings()  {return &fStrings;}
 *
 *
 *
 * in mysteer.steer (to be read in by SteerParser):
 *=====================================================
 * MySteer(){
 *
 *  fMyFloat   = 7.;
 *
 *  fArrayI    = "6, 8, 9, 23";
 *
 *  fOneString = "theOne";
 *  fStrings   = "first,second, third";
 *
 * }
 *
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <stdlib.h>

#include <TClass.h>
#include <TDataMember.h>
#include <TMethodCall.h>
#include <TObjString.h>
#include <TString.h> 
#include <TObjArray.h>
#include <TArrayI.h>
#include <TArrayD.h>
#include <TArrayF.h>

#include "BaseSteer.h"

ClassImp(BaseSteer)

using namespace std;

//__________________________________________________________
Int_t BaseSteer::SplitExpression(const TString& expr, TString& left,
                               TString& right) const
{
// Splits the expression exp into its left hand side and right hand side
// (assuming that there is just one operator!)
//
// Returns the type ID (0 if not successful)
// 1: =
//
// In future updates one could use the negative numbers to indicate that
// there are additional operators in the expression
    Int_t index,exptype=0;
    
    if ( ( index=expr.Index("=") ) >= 0 ) { // assignment

        left  = expr(0,index);
        right = expr(index+1,expr.Length()-index-2);
        
        exptype=1;
    }
    return exptype;
}

//__________________________________________________________
void BaseSteer::SetValues(TString& block)
{
// Default implementation that supports simple assignments
//
// Loops over all expressions in the block and for each datamember
// that appears in this block, SetValues calls the corresponding setter
// method.
//
// All this is done via the ROOT RTTI system, so be sure to implement
// setters for all datamembers in derived classes as non-inline functions!
//
// In case a derived class needs more sophisticated operations than
//
//   fDatamember = value;
//
// one has to override SetValues
    Int_t index,type;
    TString expression,variable,value;
    TClass *myclass = this->IsA();
    TDataMember *thevar;
    TMethodCall *setter;

    while ( ( index=block.Index(";") ) >= 0 ) { // expressions left
        expression = block(0,index+1);
        block = block(index+1,block.Length()-index-1); // rest after ';'

            // split the expression into left hand side and right hand side
        type = SplitExpression(expression,variable,value);

        if ( type != 1 ) continue; // only assignments are treated so far

	Bool_t setVariable = kFALSE;
        thevar = myclass->GetDataMember(variable);

        if (thevar) { // variable in this class
	    setter = thevar->SetterMethod(myclass);
	    //setter = thevar->SetterMethod();

	    if (setter) {
	        if(value.BeginsWith("\"") && value.EndsWith("\"")){
		    // work around for strings longer than CINT-limit: transfer pointer only
  		    value[value.Length()-1] = '\0'; // ...but cut last and first '"'
		    setter->Execute(this,Form("(const char*) %p", value.Data()+1));
		} else {
		    setter->Execute(this,value);
		}
		setVariable = kTRUE;
	    }
        } else {      // variable might be in base class
            TClass *base = myclass->GetBaseDataMember(variable);
            if (base) { // variable in a base class
                thevar = base->GetDataMember(variable);
		setter = thevar->SetterMethod(base);
		//setter = thevar->SetterMethod();

		if (setter) {
		    if(value.BeginsWith("\"") && value.EndsWith("\"")){
		        // work around for strings longer than CINT-limit: transfer pointer only
		        value[value.Length()-1] = '\0'; // ... but cut last and first '"'
			setter->Execute(this,Form("(const char*) %p", value.Data()+1));
		    } else {
		        setter->Execute(this,value);
		    }
		    setVariable = kTRUE;
		}
	    }
        }
	if (!setVariable) {
	    this->Error("SetValues", "cannot set member %s", variable.Data());
	}
    } // loop over expressions in block 
}

//__________________________________________________________
void BaseSteer::SplitString(const char* str, const char* splitStr,
			  TObjArray* result) const
{
// Splits up 'str' into an TObjArray ('result') of substrings (TObjString).
// First, all entries in 'result' are deleted, then the new strings
// are created on the heap. They have to be deleted by the user when not
// used anymore.
//
// The character 'splitStr' is used as substring delimiter (note:
// spaces are not allowed!)
//
// Example:
//
// str      == "value1,value2,value3"
// splitStr == ","
// -> result contains three TObjStrings: "value1", "value2", "value3"

  result->Delete();
  char* buffer = new char[strlen(str)+1];
  strcpy(buffer, str);
  char* p = buffer;

  while (strlen(p) > 0) {
//    while (p[0] == ' ') p++;
    char* name = p;
    UInt_t pos = strcspn(p, splitStr);
    if (pos < strlen(p)) {
      p[pos] = 0;
      p += pos + 1;
    } else
      p += strlen(p);
    if (strlen(name)) result->Add(new TObjString(name));
  }
  delete[] buffer;
}

//__________________________________________________________
Int_t BaseSteer::StringToArray(const char* str, TArrayI& result) const
{
// Fill integer values in 'str' into 'result', returns new length of result.
//
// Every non-digit may serve as delimiter between two numbers, but take care that
// white spaces will be removed while parsing by SteerParser!
// 
// The conversion is done by strtol, for details on format cf. 'man strtol'.
//
// Example:
//
// str      == "1,78t-6"
// -> result will have the length 3 with 
//    result[0] == 1, result[1] == 78, result[2] == -6

  char* buffer = new char[strlen(str)+1];
  strcpy(buffer, str);
  char* p = buffer;
  Int_t entries = 0;

  while (strlen(p) > 0){
    char* nextValue = p;
    if(result.GetSize() <= entries) result.Set(2 * result.GetSize() + 20);
    result[entries++] = strtol(nextValue, &p, 10);   // post '++' !
    if(*p == '\0' && *nextValue != '\0') break; // already end of whole string!
    if(nextValue == p++){ // post '++' !
      --entries;
      this->Error("StringToArray", "Ignore incorrect char '%s'!", 
		  TString(*(p-1)).Data());
    }
  }

  delete[] buffer;
  result.Set(entries); // shorten to desired length
  return result.GetSize();
}


//__________________________________________________________
Int_t BaseSteer::StringToArray(const char* str,  TArrayD& result) const
{
// Fill double/float values in 'str' into 'result', returns new length of result.
//
// Every non-digit may serve as delimiter between two numbers, but take care that
// white spaces will be removed while parsing by SteerParser!
// 
// Conversion is done by strtod and 'man strtod' shows which format is expected:
//
//        The expected form of the string is optional leading  white
//        space  as  checked by isspace(3), an optional plus (``+'')
//        or minus sign (``-'') followed by  a  sequence  of  digits
//        optionally  containing  a decimal-point character, option-
//        ally followed by an exponent.  An exponent consists of  an
//        ``E''  or  ``e'',  followed  by  an optional plus or minus
//        sign, followed by a non-empty sequence of digits.  If  the
//        locale  is  not  "C"  or "POSIX", different formats may be
//        used.
//
//
// Example: (works at least at DESY site...)
//
// str      == "1.1,-78e6,1.34E-2"
// -> result will have the length 3 with 
//    result[0] == 1.1, result[1] == -7.8 * 10^7, result[2] == 0.0134

  char* buffer = new char[strlen(str)+1];
  strcpy(buffer, str);
  char* p = buffer;
  Int_t entries = 0;

  while (strlen(p) > 0){
    char* nextValue = p;
    if(result.GetSize() <= entries) result.Set(2 * result.GetSize() + 20);
    result[entries++] = strtod(nextValue, &p);   // post '++' !
    if(*p == '\0' && *nextValue != '\0') break; // already end of whole string!
    if(nextValue == p++){ // post '++' !
      --entries;
      this->Error("StringToArray", "Ignore incorrect char '%s'!", 
		  TString(*(p-1)).Data());
    }
  }

  delete[] buffer;
  result.Set(entries); // shorten to desired length
  return result.GetSize();
}

//__________________________________________________________
Int_t BaseSteer::StringToArray(const char* str,  TArrayF& result) const
{
  // cf. StringToArray(const char* str,  TArrayD& result), 
  // but numbers will be stored in float format, not double.

  TArrayD temp;
  result.Set(this->StringToArray(str, temp));
  for(Int_t i = 0; i < temp.GetSize(); ++i){
    result[i] = static_cast<Float_t>(temp[i]);
  }
  return result.GetSize();
}

//__________________________________________________________
void BaseSteer::PrintArray(const char* description, const TArrayI& array) const
{
  // Prints "'decription' = {'array[0]', 'array[1]',... }" on screen,
  // usefull for the print method.
  cout << description << " = {";
  for(Int_t i = 0; i < array.GetSize(); ++i) {
    cout << array.At(i);
    if (i != array.GetSize()-1) cout << ", ";
  }
  cout << "}" << endl;
}

//__________________________________________________________
void BaseSteer::PrintArray(const char* description, const TArrayD& array) const
{
  // Prints "'decription' = {'array[0]', 'array[1]',... }" on screen,
  // usefull for the print method.
  cout << description << " = {";
  for(Int_t i = 0; i < array.GetSize(); ++i) {
    cout << array.At(i);
    if (i != array.GetSize()-1) cout << ", ";
  }
  cout << "}" << endl;
}

//__________________________________________________________
void BaseSteer::PrintArray(const char* description, const TArrayF& array) const
{
  // Prints "'decription' = {'array[0]', 'array[1]',... }" on screen,
  // usefull for the print method.
  cout << description << " = {";
  for(Int_t i = 0; i < array.GetSize(); ++i) {
    cout << array.At(i);
    if (i != array.GetSize()-1) cout << ", ";
  }
  cout << "}" << endl;
}
