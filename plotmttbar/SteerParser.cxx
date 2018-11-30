
/********************************************************************
 * 
 * @Package: Plotter
 *
 * Class to parse steering files in C++. Classes have to be enclosed
 * in curly brackets and start with <classname>. Example:
 * SteerExample(){ .... };
 * Steering blocks can also be named, e.g. SteerExample("default") and
 * SteerExample("modify").
 * The different versions of the class can then be obtained with 
 * GetSteer(TClass* cl) or GetSteer(TClass* cl, const char* name).
 *
 *******************************************************************/

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <TClass.h>
#include <TROOT.h>
#include <RVersion.h>
#include <TObjString.h>
#include <TList.h>
#include <TSystem.h>

#include "SteerParser.h"
#include "BaseSteer.h"

using namespace std;

//__________________________________________________________
SteerParser::SteerParser()
  : fFilename((Ssiz_t)1024),fNamespace((Ssiz_t)256),fClassname((Ssiz_t)1024), fObjectname((Ssiz_t)1024), fBlock((Ssiz_t)1024), fLine((Ssiz_t)1024),fInNamespace(kFALSE),fInBlock(kFALSE),fBraces(0),fArrayLength(0),fNameSpaces(0),fCurrentNameIndex(-1)
{
  // constructor
}

//__________________________________________________________
SteerParser::~SteerParser()
{
// Deletes the list of steerings

    fNameSpaces.Delete();
    fSteerings.Delete();
}

//__________________________________________________________
void SteerParser::SetFilename(const char* filename)
{
// Searches for the file 'filename' in $PWD
  
  char *tmpFileName = gSystem->ExpandPathName("$(PWD)/");
  fFilename = tmpFileName;
  delete tmpFileName;
  if(fFilename !="/") {
    fFilename.Append(filename);
    cout <<"SteerParser: Searching for file "<< fFilename<<endl;
    ifstream steer1(fFilename,ios::in);
    if(!steer1.fail()) {
      cout <<"SteerParser: File "<< fFilename<<" opened for reading"<<endl;
      return;
    }
  }
  cout <<"SteerParser: Could not find "<< filename
       <<" . "<< endl;  
  return;
 
}

//__________________________________________________________
void SteerParser::SetClassname(const char* classname, const char* objectname) 
{
  //Set class name with optional object name argument
  fClassname=classname;
  fObjectname=objectname;
  if(fObjectname=="") 
    cout <<"SteerParser: Classname set to "<< fClassname<<endl;
  else
    cout <<"SteerParser: Classname set to " << fClassname 
   	 << "(\"" << fObjectname <<"\")" << endl;
}

//__________________________________________________________
Bool_t SteerParser::ParseFile(const char* filename)
{
// This method divides the steering file into blocks
// and creates BaseSteer objects for each block found.
// If no prior namespace sets the current namespace to the first found.
    SetFilename( filename);
    Bool_t status = kTRUE;

    cout << "\n======== SteerParser: Reading from file '"<<fFilename<<"' ============="<<endl;
    
    ifstream steer(fFilename,ios::in); // open steering file

    if ( steer.is_open() ) {
        TString line;
	const Int_t numSteerBefore = fSteerings.GetEntriesFast();

        while ( ! steer.eof() ) {
            line.ReadLine( steer );
            if ( (line.Length() == 0) || (line.BeginsWith("//")) ) continue;

            RemoveComments(line);
            RemoveWhiteSpace(line);

            cout << "----- "<< line << endl; // we want to have that in log files

            ParseLine(line);
        }
	if(fInNamespace){
	    Warning("ParseFile","missing closing brace of namespace %s!", fNamespace.Data());
	}
	Print(fFilename, numSteerBefore, fSteerings.GetLast()); // print summary at the end
    } else { // steering file not open
        Warning("ParseFile","Could not open steering file %s!\nWill use defaults!",filename);
        status = kFALSE;
    }
    //cout << "======== SteerParser: Done reading from file '"<<fFilename<<"' ========\n"<<endl;

    if(fCurrentNameIndex == -1) fCurrentNameIndex = (fNameSpaces.GetEntriesFast() ? 0 : -1);
    return status;
}

//__________________________________________________________
Float_t SteerParser::ReadFloat(const char* keyword)
{
  // Get value as a float
  Int_t index;
  Float_t val=0;
  if(FindArrayInFile(keyword)) {
    // Statement found
    if ( ( index=fLine.Index(",") ) >= 0 ) {
      TString value=fLine(0,index);// 0=beginning, index=length (not end!)
      //convert value to a float
      sscanf(value,"%f",&val);
      //      cout <<"keyword found " <<keyword<<" = "<<val<<endl;
      return val;
    }
  }
  return val;
}

//__________________________________________________________
Int_t SteerParser::ReadInt(const char* keyword)
{
  // Get value as an int
  Int_t index;
  Int_t val=0;
  if(FindArrayInFile(keyword)) {
    // Statement found
    if ( ( index=fLine.Index(",") ) >= 0 ) {
      TString value=fLine(0,index);// 0=beginning, index=length (not end!)
      //convert value to an int
      sscanf(value,"%i",&val);
      //      cout <<"keyword found " <<keyword<<" = "<<val<<endl;
      return val;
    }
  }
  return val;
}

//__________________________________________________________
TString SteerParser::ReadString(const char* keyword)
{
  // Get value as a string
  Int_t index;
  TString value("");
  if(FindArrayInFile(keyword)) {
    // Statement found
    if ( ( index=fLine.Index(",") ) >= 0 ) {
      value=fLine(0,index);// 0=beginning, index=length (not end!)
      //      cout <<"keyword found " <<keyword<<" = "<<value<<endl;
      return value;
    }
  }
  return value;
}

//__________________________________________________________
Int_t SteerParser::ReadArray(const char* keywordin,Float_t array[])
{
  // Split line up into a list of floats and store results in array
  TString value("");
  Float_t val;
  Int_t indarray=0;
  Int_t index=0;
  if(FindArrayInFile(keywordin)) {
  // Statement found
    while( fLine.Length()>0){
    
      if ( ( index=fLine.Index(",") ) >= 0 ) {
	value=fLine(0,index);// 0=beginning, index=length (not end!)
	fLine = fLine(index+1,fLine.Length()-index-1); //reduce line
	//convert value to a float
	sscanf(value,"%f",&val);
	array[indarray]=val;
	if(indarray<fArrayLength)array[indarray]=val;
	else if(fArrayLength>0) {
	  cout<<"SteerParser: Mismatch in declared and actual number of elements and for Array" <<keywordin<< endl;

	  cout<<"SteerParser: Array index=" <<indarray<< 
	      ", Array size=" <<fArrayLength<< endl;
	}

	//	cout <<"found value" <<val<<endl;
	indarray++;
      }
    }
  }
  if(indarray>0&&indarray<fArrayLength&&fArrayLength>0){
    cout<<"SteerParser: Mismatch in declared and actual number of elements and for Array" <<keywordin<< endl;

    cout<<"SteerParser: Array index=" <<indarray<< 
      ", Array size" <<fArrayLength<< endl;
  }
  return fArrayLength;

}

//__________________________________________________________
Int_t SteerParser::ReadArray(const char* keywordin,Int_t array[])
{
  // Split line up into a list of ints and store results in array
  TString value("");
  Int_t val;
  Int_t indarray=0;
  Int_t index=0;
  if(FindArrayInFile(keywordin)) {
  // Statement found
    while( fLine.Length()>0){
    
      if ( ( index=fLine.Index(",") ) >= 0 ) {
	value=fLine(0,index);// 0=beginning, index=length (not end!)
	fLine = fLine(index+1,fLine.Length()-index-1); //reduce line
	//convert value to a float
	sscanf(value,"%i",&val);
	array[indarray]=val;
	if(indarray<fArrayLength)array[indarray]=val;
	else if(fArrayLength>0) {
	  cout<<"SteerParser: Mismatch in declared and actual number of elements and for Array" <<keywordin<< endl;
	  cout<<"SteerParser: Array index=" <<indarray<< 
	      ", Array size" <<fArrayLength<< endl;
	}

	//	cout <<"found value" <<val<<endl;
	indarray++;
      }
    }
  }
  if(indarray>0&&indarray<fArrayLength&&fArrayLength>0) {
    cout<<"SteerParser: Mismatch in declared and actual number of elements and for Array" <<keywordin<< endl;
    cout<<"SteerParser: Array index=" <<indarray<< 
      ", Array size" <<fArrayLength<< endl;
  }

  return fArrayLength;

}

//__________________________________________________________
Int_t SteerParser::ReadArray(const char* keywordin, TObjArray& array, Bool_t checkNamespace)
{
  // Split line up into a list of TObjStrings and store results in array

  TString value("");
  TObjArray val;
  Int_t indarray=0;
  Int_t index=0;

  if(FindArrayInFile(keywordin, checkNamespace)) {
    // Statement found
    while( fLine.Length()>0){
    
      if ( ( index=fLine.Index(",") ) >= 0 ) {
	value=fLine(0,index);// 0=beginning, index=length (not end!)
	fLine = fLine(index+1,fLine.Length()-index-1); //reduce line
	value.ReplaceAll("};", "");
	value.ReplaceAll("\"", "");
	value.ReplaceAll(";", "");
	val.AddLast(new TObjString(value));
	indarray++;
      }
    }
  }
  array.Clear();
  array.AddAll(&val);

  if(indarray>0&&indarray<fArrayLength&&fArrayLength>0){
    cout<<"SteerParser: Mismatch in declared and actual number of elements and for Array " <<keywordin<< endl;

    cout<<"SteerParser: Array index=" <<indarray<< 
      ", Array size" <<fArrayLength<< endl;
  }
  return fArrayLength;
}

//__________________________________________________________
void SteerParser::ReadFastArray(const char* keyword,Float_t array[])
{
  // Put the array in the file into the array array.

    TString Keyword(keyword); 
    fBraces=0;
    Bool_t classfound=kFALSE;
    if(fFilename=="") {
      cout << "SteerParser: You  must specify a file first with SetFilename()"<<endl;
      return;
    }
      
    ifstream steer(fFilename,ios::in); // open steering file

    if ( steer.is_open() ) {

      while ( ! steer.eof() ) {
	fLine.ReadLine( steer );
	if ( (fLine.Length() == 0) || (fLine.BeginsWith("//")) ) continue;

	RemoveComments(fLine);
	RemoveWhiteSpace(fLine);
	//Look for class heading
	if(FindClassnameInLine(fLine)) {
	  //	  cout <<"SteerParser: Found class "<< fClassname <<endl;
	  classfound=kTRUE;
	}
	if(classfound) {
	  //Count braces to check we are still in class description
	  CountBraces(fLine);
	  if(FindKeywordInLine(fLine,Keyword)) {

	    // now know the array size, read in that many floats from the file.
	    Int_t counter = 0;
	    Float_t tmp = 0;
	    while(( counter != fArrayLength ) && ! steer.eof()){
	      steer >> tmp;
	      array[counter] = tmp;
	      counter++;
	    }
	    
	    // We have read all the information so return
	    return;
	  }
	  // Keep these error messages
	  if(fBraces==0&&classfound) {
	    cout<<"SteerParser: Could not find variable "<< 
	      keyword <<" in class "<<fClassname<<endl;
	    return;
	  }
	}
      }
      if(!classfound) {
	cout<<"SteerParser: Could not find class "<<fClassname 
	    <<" with object "<<fObjectname<<endl;
	return;
      }
    } else { // steering file not open
      cout << "SteerParser: Could not open steering file "<<fFilename <<endl;
      return;
    }
    cout <<"SteerParser:Could not find keyword "<< 
      keyword<< " in file " <<fFilename<<endl;
    if(fBraces !=0) cout<<"SteerParser: Braces in file "<< 
		      fFilename <<" don't match "<<fBraces<<endl;
    
    return;

  return;
}

//__________________________________________________________
void SteerParser::ReadFastArray(const char* keyword,Int_t array[])
{
    TString Keyword(keyword); 
    fBraces=0;
    Bool_t classfound=kFALSE;
    if(fFilename=="") {
      cout << "SteerParser: You  must specify a file first with SetFilename()"<<endl;
      return;
    }
      
    ifstream steer(fFilename,ios::in); // open steering file

    if ( steer.is_open() ) {

      while ( ! steer.eof() ) {
	fLine.ReadLine( steer );
	if ( (fLine.Length() == 0) || (fLine.BeginsWith("//")) ) continue;

	RemoveComments(fLine);
	RemoveWhiteSpace(fLine);
	//Look for class heading
	if(FindClassnameInLine(fLine)) {
	  //	  cout <<"SteerParser: Found class "<< fClassname <<endl;
	  classfound=kTRUE;
	}
	if(classfound) {
	  //Count braces to check we are still in class description
	  CountBraces(fLine);
	  if(FindKeywordInLine(fLine,Keyword)) {

	    // now know the array size, read in that many ints from the file.
	    Int_t counter = 0;
	    Int_t tmp = 0;
	    while(( counter != fArrayLength ) && ! steer.eof()){
	      steer >> tmp;
	      array[counter] = tmp;
	      counter++;
	    }
	    // We have read all the information so return
	    return;
	  }
	  // Keep these error messages
	  if(fBraces==0&&classfound) {
	    cout<<"SteerParser: Could not find variable "<< 
	      keyword <<" in class "<<fClassname<<endl;
	    return;
	  }
	}
      }
      if(!classfound) {
	cout<<"SteerParser: Could not find class "<<fClassname 
	    <<" with object "<<fObjectname<<endl;
	return;
      }
    } else { // steering file not open
      cout << "SteerParser: Could not open steering file "<<fFilename <<endl;
      return;
    }
    cout <<"SteerParser:Could not find keyword "<< 
      keyword<< " in file " <<fFilename<<endl;
    if(fBraces !=0) cout<<"SteerParser: Braces in file "<< 
		      fFilename <<" don't match "<<fBraces<<endl;
    
    return;

  return;

}

//__________________________________________________________
Bool_t SteerParser::FindArrayInFile(const char* keywordin, Bool_t checkNamespace)
{
//
// Finds an array or single value in the steering file and stores numbers/strings
// separated by commas in fLine
//
// If(checkNamespace) only takes into account lines from the namespace 'fNamespace'
    TString keyword(keywordin); 
    fBraces=0;
    if(fFilename=="") {
        Warning("FindArrayInFile","You  must specify a file first with SetFilename()");
        return kFALSE;
    }
      
        //    cout << "\n======== SteerParser: Reading from file '"<<fFilename<<"' ============="<<endl;
    
    ifstream steer(fFilename,ios::in); // open steering file

    if ( steer.is_open() ) {
        Bool_t classfound=kFALSE;
	TList lines;
        while ( ! steer.eof() ) {
	    TObjString *newLine = new TObjString;
	    newLine->String().ReadLine( steer );
	    lines.Add(newLine);
	}

	if(checkNamespace && !fNamespace.IsNull()) this->SelectNamespace(lines, fNamespace);

	TIter nextLine(&lines);
	while(TObject* line = nextLine()){
	    fLine = line->GetName(); // content of string
	    

            if ( (fLine.Length() == 0) || (fLine.BeginsWith("//")) ) continue;
            RemoveComments(fLine);
            RemoveWhiteSpace(fLine);
                //Look for class heading
            if(FindClassnameInLine(fLine)) {
                classfound=kTRUE;
            }
            if(classfound) { // true for all lines after class has been found
                    //Count braces to check we are still in class description
                CountBraces(fLine);
                if(FindKeywordInLine(fLine,keyword)) {
                        //If line does not end in ';' Keep reading till we find it
                    while( ! LineEnd(fLine) ){
			TObject* newNextLine = nextLine();
			if(!newNextLine) break;
			TString linetmp = newNextLine->GetName();
                        if ( (linetmp.Length() == 0) 
                             || (linetmp.BeginsWith("//")) ) continue;
                        RemoveComments(linetmp);
                        RemoveWhiteSpace(linetmp);
                        CountBraces(fLine);
                        fLine +=linetmp;
                    }
                        //line should now contain a list of strings
                    fLine +=","; //add comma to end of line to simplify string separations
		  
                        // We have read all the information so return
                    return kTRUE;
                }
                    // if block is closed and keyword not found print out an error
                    // second condition is needed when "SteerClass()" is on a line by itself
                    // if the steering has no name, line is just ')', else it's empty
                if(fBraces==0 &&
                   ! ( fLine.CompareTo(")")==0 || fLine.CompareTo("")==0 ) ) {
                    Warning("FindArrayInFile","Could not find variable %s in class %s",
                            keywordin,fClassname.Data());
                    return kFALSE;
                }
            }
        }
	lines.Delete();
        if(!classfound) {
            Warning("FindArrayInFile","Could not find class %s with object %s",
                    fClassname.Data(),fObjectname.Data());
	    return kFALSE;
        }
    } else { // steering file not open
        Warning("FindArrayInFile","Could not open steering file %s",fFilename.Data());
        return kFALSE;
    }
    Warning("FindArrayInFile","Could not find keyword %s in file %s",keywordin,fFilename.Data());
    if(fBraces !=0)
        Warning("FindArrayInFile","Braces in file %s don't match: %d",fFilename.Data(),fBraces);

    return kFALSE;
}

//__________________________________________________________
Bool_t SteerParser::FindKeywordInLine(TString& line, const TString& keyword)
{
  //
  // Searches for a keyword in the line and returns the remainder of the line
  // Returns false if line does not contain keyword
  //
    Int_t   index;
    TString keywordtest;
    fArrayLength=0;
    //    if ( ! fInBlock ) { // we're not yet inside a block (enclosed in {})
    if ( ( index=line.Index("=") ) >= 0 ) { // found a keyword = statement
      // check charecters  before "=" sign for compatibility with keyword
      keywordtest = line(0,index);
      if(keywordtest==keyword) {
	line = line(index+1,line.Length()-index-1);
	return kTRUE;
      }
      Int_t   index1=line.Index("[");
      Int_t   index2=line.Index("]");
      // perhaps we have an array with declared # elements (eg arr[99]= ...)
      if ( index1>= 0 && index2>index1 ) {
	// we have found the opening and closing braces
	keywordtest = line(0,index1);
	if(keywordtest==keyword) {
	  //Now read the string between the braces and convert into
	  // an integer which denotes the size of the array
	  TString arraylength=line(index1+1,index2-index1-1);
	  line = line(index+1,line.Length()-index-1);
	  //COnvert string to integer
	  if(arraylength.Length()>0) sscanf(arraylength,"%i",&fArrayLength);
	  else fArrayLength=0;

// 	  cout <<"SteerParser: Array "<<keyword<<" found ";
// 	  if(fArrayLength>0)  cout << " with length " <<fArrayLength <<endl;
// 	  else cout << " with unspecified length "<<endl;
  
	  return kTRUE;
	}
      }
    }
    return kFALSE;
}

//__________________________________________________________
Bool_t SteerParser::FindClassnameInLine(TString& line) const
{
  //
  // Searches line for a classname and object name (if specified) to match
  // those given in fClassname and fObjectname
  // class and object names are stripped so that line afterwards contains the
  // remainder from ')' to the end
  // Returns false if line does not contain classname and object name
  //
    Int_t   index;
    TString keywordtest;
    if ( ( index=line.Index("(") ) >= 0 ) { // found a keyword ( statement
      // check characters  before "(" sign for compatibility with classname
      keywordtest = line(0,index);
      if(keywordtest==fClassname) {
	line = line(index+1,line.Length()-index-1);
	//no Object name specified so block is already found
	if(fObjectname=="") return kTRUE;
	//Now check Object name
	else if ( line.Length()>2&&( index=line.Index(")") ) >= 0 ) { 
	  // found a keyword ) statement
	  // Do not compare 1st +last elements as these contain " characters
	  keywordtest = line(1,index-2); 
	  if(keywordtest==fObjectname) {
	    line = line(index+1,line.Length()-index-1);
	    return kTRUE;
	  }
	}
      }
    }
    return kFALSE;
}
//__________________________________________________________
Bool_t SteerParser::LineEnd(const TString& line) const
{
  // true if line ends in a ';' symbol
  return line.Index(";")  >= 0 ?  kTRUE :  kFALSE;
}

//__________________________________________________________
TString SteerParser::ValueOfLine(const TString& line) const
{
// 
//  Copies the contents of a string until it reaches the ';' character
//  This is usually the value of the parameter or for arrays a set
//  of values separated by commas
//
  Int_t   index;
  TString value("");
  if ( ( index=line.Index(";") ) >= 0 ) { // Look for ';' 
    value = line(0,index);// 0=beginning, index=length (not end!)
  }
  else {
    cout <<"SteerParser::ValueOfLine: line does not end with ;" << endl;
    cout <<"                             Check Steering file " <<endl;
  }
  return value;
}


//__________________________________________________________
void SteerParser::ParseLine(TString& line)
{
// Splits lines into components, calls different methods
// according to the internal state
  
    Int_t   index;

    // First try to find a namespace
    if( ! fInNamespace ) { 
        if(line.BeginsWith("namespace")){
	    index = line.Index("{");
	    if(index == kNPOS) index = line.Length(); // '{' in next line, take care below
	    fNamespace = line(9, index - 9); // 9=beginning, index-9=length
	    if(fNamespace.Length()) {
	        if(!fNameSpaces.FindObject(fNamespace)) 
		    fNameSpaces.Add(new TObjString(fNamespace));
	    } else {
	        Warning("ParseLine", "Expect a name between 'namespace' and '{', use global!");
	    }
	    fInNamespace = kTRUE;
	    line = line(index+1,line.Length()-index-1); // rest after '{'
	    if ( line.Length() == 0 ) return; // line doesn't contain other parts else go on
	}
    } else if( ! fInBlock && line[0] == '}'){ // ending namespace
        fNamespace = "";
	fInNamespace = kFALSE;
	if(line.Length() != 1) ParseLine(line.Remove(0, 1)); // parse the rest after '}'
	return;
    }

    TString expression;
    // now we can deal with the steerings
    if ( ! fInBlock ) { // we're not yet inside a block (enclosed in {})
        if(fInNamespace && line[0] == '{') { // assuming starting of namespace
	    // might be wrong, but than it is a wrong '{' anyway!
	    if(line.Length() == 1) return; // line doesn't contain other parts
	    else line.Remove(0, 1);        // remove '{' and go on
	}

        if ( ( index=line.Index("(") ) >= 0 ) { // classname
            fClassname = line(0,index); // 0=beginning, index=length (not end!)
            line = line(index+1,line.Length()-index-1); // rest after '('

            if ( line.Length() > 0 ) ParseLine(line); // line might contain other parts
        } else if ( ( index=line.Index(")") ) >= 0 ) { // beginning of block
            fInBlock = kTRUE;
            if (index == 0){ // unnamed steering
                line = line(index+1,line.Length()-index-1); // rest after ')'

                if ( line.Length() > 0 ) ParseLine(line); // line might contain other parts
            } else {
                if (line[0]=='"' && line[index-1]=='"'){ // name in double ticks
                    fObjectname = line(1,index-2); // extract name

                } else {
                    Warning("ParseLine",
			    "Expected a name like '%s' in quotation marks! Don't use a name...",
			    line(1,index-1).Data());
                }
                line = line(index+1,line.Length()-index-1); // rest after ')'
                
                if ( line.Length() > 0 ) ParseLine(line); // line might contain other parts
            }
        } else {
            Warning("ParseLine","Line '%s' makes no sense here! Skipping...",line.Data());
        }
    } else { // we're inside a block
	Int_t indOpen  = line.Index("{"); 
	Int_t indClose = line.Index("}"); 
	Int_t indSemi  = line.Index(";"); 
        if ( indOpen >= 0 && (indClose == kNPOS || indOpen < indClose)
	     && (indSemi == kNPOS || indOpen < indSemi)) { // starting '{'
	    index = indOpen;
            line = line(index+1,line.Length()-index-1); // rest after '{'
            
            if ( line.Length() > 0 ) ParseLine(line); // line might contain other parts
        } else if ( indSemi >= 0 && (indClose == kNPOS || indSemi < indClose)) { // expression
	    index = indSemi;
            expression = line(0,index+1);
            fBlock.Append(expression); // add expression to block
            line = line(index+1,line.Length()-index-1); // rest after ';'

            if ( line.Length() > 0 ) ParseLine(line); // line might contain other parts
        } else if ( indClose >= 0 ) { // end of block
	    index = indClose;
	    CreateSteering(fClassname,fObjectname,fBlock,fNamespace);
            fInBlock = kFALSE;

            fBlock="";
            fClassname="";
            fObjectname="";
	    if(index < line.Length()-1){
	        line = line(index+1,line.Length()-index-1); // rest after '}'
		ParseLine(line);
	    }

        } else { // assume that expression is continued on next line
            fBlock.Append(line); // add entire line to block
        }
    }
}

//__________________________________________________________
Bool_t SteerParser::CreateSteering(TString& theclassname, TString& objectname, TString& block, const char* nameSpace)
{
// Gets called by ParseLine when the closing brace of a block is found
// or might be directly called by a user
//
// First parameter is the classname of a concrete steering class
// e.g. BaseSteerCreateScatElec
//
// Second parameter is the name of the steering object
//
// Third parameter is a block of expressions to be evaluated by the
// concrete steering class
// e.g. 'value1=0.8; value2=4711; value3=42.0;' and so on...
//
// Optional fourth parameter is a namespace this steering should belong to.
// (Looking for namespace might cause trouble for very long strings IF this method 
//  was called 'standalone', because in the long string work around in BaseSteer::SetValues
//  it is assumed that the correct classname etc. are in the corresponding data members
//  of SteerParser.)

    TClass *myclass;

    RemoveWhiteSpace(theclassname);
    RemoveWhiteSpace(objectname);
    RemoveWhiteSpace(block);
    
    myclass = gROOT->GetClass(theclassname);

    if ( myclass ) {
        BaseSteer *mysteer = static_cast<BaseSteer*>( myclass->New() );
        mysteer->SetValues(block);
        mysteer->SetName(objectname);
	if(nameSpace) mysteer->SetTitle(nameSpace);
        fSteerings.AddLast( mysteer );
        return kTRUE;
    } else {
        Warning("CreateSteering","Cannot find steering class '%s'! Skipping...",
                theclassname.Data() );
        return kFALSE;
    }    
}

//__________________________________________________________
BaseSteer* SteerParser::GetSteer(TClass* cl)
{
// If a steering corresponding to class 'cl' and without a name 
// is in the list of steerings, it is simply returned.
// But if namespaces are defined, a steering from that namespace is preferred
// to a steering without namespace, which is the 'fall back' solution.
//
// If no corresponding steering object exists, a default one is created
// and appended to the list of steerings. Default values are printed to
// screen.
    BaseSteer* noNameSpaceSteer = NULL;

    TIter next(&fSteerings);

    while ( BaseSteer* steer = (BaseSteer*)next() ) {
        if ( steer->IsA() == cl && strlen(steer->GetName()) == 0 ){// GetName never is NULL!
	    if(fCurrentNameIndex == -1){ // no namespace required
                return steer;
	    } else {
	        const char* nameSpace = steer->GetTitle();
		if(nameSpace && !strcmp(nameSpace, fNameSpaces[fCurrentNameIndex]->GetName())){
		    return steer;
		} else if(!nameSpace || !strcmp(nameSpace, "")){// no namescpace of steering
		    noNameSpaceSteer = steer; //if no matching is found, take this without name
		}
	    }
	}
    }
    if(noNameSpaceSteer) return noNameSpaceSteer;

    BaseSteer* thesteer = (BaseSteer*)cl->New();

    cout << endl
         << "====== SteerParser: Using default values for class '"
         << cl->GetName()
         <<"'  ======\n";
    thesteer->Print();
    cout << "====== SteerParser: End of defaults for class '"
         << cl->GetName()
         <<"'  ===========\n"
         << endl;
    
    fSteerings.AddLast(thesteer);
    
    return thesteer; // return default object if not explicitly given in file
}

//__________________________________________________________
BaseSteer* SteerParser::GetSteer(TClass* cl, const char* name)
{
// If a steering corresponding to class 'cl' is in the list of steerings
// it is returned if the name matches.
// But if namespaces are defined, a steering from that namespace is preferred
// to a steering without namespace, which is the 'fall back' solution.
//
// If no corresponding steering object exists, NULL is returned.

    BaseSteer* noNameSpaceSteer = NULL;

    TIter next(&fSteerings);

    while ( BaseSteer* steer = (BaseSteer*)next() ) {
        if ( steer->IsA() == cl && !strcmp(steer->GetName(), name) ){
	    if(fCurrentNameIndex == -1){ // no namespace required
                return steer;
	    } else {
	        const char* nameSpace = steer->GetTitle();
		if(nameSpace && !strcmp(nameSpace, fNameSpaces[fCurrentNameIndex]->GetName())){
		    return steer;
		} else if(!nameSpace || !strcmp(nameSpace, "")){// no namescpace of steering
		    noNameSpaceSteer = steer; //if no matching is found, take this without name
		}
	    }
	}
    }

    if(!noNameSpaceSteer) 
        Error("GetSteer", "No steering type %s with name %s!", cl->GetName(), name);
    return noNameSpaceSteer; // maybe NULL
}


//__________________________________________________________
void SteerParser::RemoveComments(TString& line) const
{
// Removes comments from line (everything after '//')
    Int_t index;
    
    if ( ( index=line.Index("//") ) >= 0 ) { // line has comment
        line = line(0,index);
    }
}

//__________________________________________________________
void SteerParser::RemoveWhiteSpace(TString& line) const
{
// Removes all spaces and tabs from line
    Int_t index;
    
    while ( ( index=line.Index(" ") ) >= 0 ) { // space at position index
        line = line.Remove(index,1);
    }
    while ( ( index=line.Index("\t") ) >= 0 ) { // tab at position index
        line = line.Remove(index,1);
    }
}
//__________________________________________________________
void SteerParser::CountBraces(TString & line)
{
//
// Counts and removes curly braces from line
// fBraces is the running total of open - close  braces 
// so fBraces=0 is a complete block
    Int_t index;
    while ( ( index=line.Index("{") ) >= 0 ) { // space at position index
        line = line.Remove(index,1);
	fBraces++;
    }
    while ( ( index=line.Index("}") ) >= 0 ) { // space at position index
        line = line.Remove(index,1);

	fBraces--;
    }
}

//__________________________________________________________
Int_t SteerParser::GetNumSteerings(TClass* theclass)
{
// Returns the number of accepted steerings of type theclass
// (the method without argument returns the number of all accepted steerings)
    Int_t i=0, num=0, max=fSteerings.GetEntriesFast();

    for (i=0; i<max; ++i ){
        if (fSteerings.UncheckedAt(i)->IsA() == theclass ) ++num;
    }
    return num;
}

//__________________________________________________________
void SteerParser::Print(Option_t* option) const
{
// Prints the given steering values to stdout, taking 'option' as filename if given
    this->Print(option, 0, fSteerings.GetLast());
}

//__________________________________________________________
void SteerParser::Print(Option_t* option, Int_t first, Int_t last) const
{
// Prints the values of steerings 'first' to 'last' to stdout, 
// taking 'option' as filename if given
    last = TMath::Min(last, fSteerings.GetLast());
    cout << "====== SteerParser: Accepted values";
    if(option && strlen(option)) cout << " from file '" << option << "'";
       cout << " ======" << endl;
    if(fNameSpaces.GetEntriesFast() == 0){      
	for(Int_t i = first; i <= last; ++i) fSteerings[i]->Print();
    } else {
        for(Int_t i = 0; i < fNameSpaces.GetEntriesFast(); ++i){
	    cout << "\n  ===== namespace '" << fNameSpaces[i]->GetName() << "':" << endl;
	    for(Int_t j = first; j <= last; ++j){
	        if(!strcmp(fSteerings[j]->GetTitle(), fNameSpaces[i]->GetName())) 
		    fSteerings[j]->Print();
	    }
	}
	cout << "\n  ===== without namespace:"<< endl;
	for(Int_t j = first; j <= last; ++j){
	    if(!strcmp(fSteerings[j]->GetTitle(), "")) fSteerings[j]->Print();
	}
    }
}

//__________________________________________________________
Int_t SteerParser::NextNameSpace()
{
  // switches to next 'namespace' of steerings
  // 1 if success, 0 if no further 'namespace' (or no 'namespaces' given at all)
  if(fCurrentNameIndex == -1) return 0;

  if(++fCurrentNameIndex < fNameSpaces.GetEntriesFast()) return 1;  // prefix increment!

  fCurrentNameIndex = -1;
  return 0;
}

//__________________________________________________________
const char* SteerParser::GetCurrentNameSpace() const
{
  // giving name of namespace that steering objects currently need to belong to,
  // if non required return NULL
  if(fCurrentNameIndex == -1) return NULL;
  else return fNameSpaces.At(fCurrentNameIndex)->GetName();
}

//__________________________________________________________
Bool_t SteerParser::SetCurrentNameSpace(const char* name)
{
  // makes namespace 'name' the current one, kFALSE if no such namespace found

  TObject* newNameSpace = fNameSpaces.FindObject(name); // it's a TObjString
  if(newNameSpace){
    fCurrentNameIndex = fNameSpaces.IndexOf(newNameSpace);
    return kTRUE;
  }else{
    return kFALSE;
  }
}

//__________________________________________________________
void SteerParser::SelectNamespace(TList& steerLines, const char* nameSpace) const
{
  // SteerLines is assumed to contain TObjStrings containing lines read in from steerig file,
  // removes all lines and part of lines that do not belong to the requested namespace.

  Bool_t inNamespace = kFALSE;
  Int_t nBraces = 0;
  
  TObjLink *lnk = steerLines.FirstLink();
  while (lnk) {                                                
    TString &line = static_cast<TObjString*>(lnk->GetObject())->String(); // reference!
    if (line.Length() == 0 || line.BeginsWith("//")) {
      lnk = this->RemoveLink(steerLines, lnk);
      continue;
    }
    this->RemoveComments(line);
    this->RemoveWhiteSpace(line);

    if(!inNamespace){
      Int_t index = line.Index("namespace");
      if(index == kNPOS){
	lnk = this->RemoveLink(steerLines, lnk);
	continue;
      }

      Int_t indexEnd = line.Index("{");
      if(indexEnd == kNPOS) indexEnd = line.Length();

      if(TString(line(index + 9, indexEnd - index - 9)) == nameSpace){
	inNamespace = kTRUE;
	nBraces = (indexEnd == line.Length() ? 0 : 1); // if no opening '{', look later...
      }

      line = line(indexEnd+1,line.Length()-indexEnd-1); // rest after '{'
      if(line.IsNull()){
	lnk = this->RemoveLink(steerLines, lnk);
      } // else just go on with rest of line
    } else { // in namespace 
      if(nBraces == 0 && line[0] != '{'){
	this->Error("SelectNamespace", "Assuming '{' after 'namespace %s', found %s!", 
		    nameSpace, line.Data());
	nBraces = 1;
      } 

      const Int_t indOpen  = line.Index("{");
      const Int_t indClose = line.Index("}");
      Int_t indEnd = kNPOS;
      if(indOpen == kNPOS && indClose == kNPOS){
	lnk = lnk->Next();
	continue; // go on
      } else if(indOpen == kNPOS) {
	indEnd = indClose;
	--nBraces;
      } else if(indClose == kNPOS || indOpen < indClose){
	indEnd = indOpen;
	++nBraces;
      } else {
	indEnd = indClose;
	--nBraces;
      }
      
      TString rest = line(indEnd+1,line.Length()-indEnd-1); // rest of line
      if(rest.IsNull()){
      } else {
	line = line(0,indEnd+1); // start: 0, length: indEnd+1
	steerLines.AddAfter(lnk, new TObjString(rest));
      }

      if(nBraces == 0){ // end of namespace!
	inNamespace = kFALSE;
	line.Remove(indEnd,1); //stays as empty string if it was the last line with '}' only...
      }

      lnk = lnk->Next();
    } // end of checking namespace or not
  } // end while 
}


//__________________________________________________________
TObjLink* SteerParser::RemoveLink(TList& list, TObjLink* link) const
{
  // Removes the link from the list, deletes its object and returns the next,
  // if link was not part of list, no delete but next is returned anyhow.
  if(!link) return NULL;

  TObjLink* newLink = link->Next();
  delete list.Remove(link); // if link is not part of list ==> delete NULL
  return newLink;
}
