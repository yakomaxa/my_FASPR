/*******************************************************************************************************************************
This file is a part of the protein side-chain packing software FASPR

Copyright (c) 2020 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/
#include "Search.h"
#include <sstream>
#include <cctype>
#include <cstdlib>

using namespace std;
string PROGRAM_PATH=(string)".";
string ROTLIB2010=(string)"dun2010bbdep.bin";
int ROTAMER_TOPN=-1;//-1 means no restriction (default FASPR behavior)
vector<SiteOverride> SITE_OVERRIDES;//per-site -n overrides, from inline "-n CHAIN:POS:N" or a "-m" site-map file
bool AROMACHI2_ENABLED=false;//LimitAromaChi2: off by default (preserves default FASPR behavior)
float AROMACHI2_MIN=70.0;
float AROMACHI2_MAX=110.0;
bool AROMACHI2_INCLUDE_TRP=false;

static bool IsAllDigits(const string&s)
{
  if(s.empty())return false;
  for(size_t k=0;k<s.size();k++){
    if(!isdigit((unsigned char)s[k]))return false;
  }
  return true;
}

static bool ParsePositiveInt(const string&s,int&out)
{
  if(!IsAllDigits(s))return false;
  long v=strtol(s.c_str(),NULL,10);
  if(v<=0)return false;
  out=(int)v;
  return true;
}

//chain token: a single character, or "_" meaning the blank/space chain id
static bool ParseChainToken(const string&s,char&chID)
{
  if(s=="_"){chID=' ';return true;}
  if(s.size()==1){chID=s[0];return true;}
  return false;
}

//residue position token: digits (optional leading '-') with an optional single trailing insertion-code letter, e.g. "15" or "15A"
static bool ParsePosIcodeToken(const string&s,int&pos,char&ins)
{
  if(s.empty())return false;
  string body=s;
  ins=' ';
  char last=s[s.size()-1];
  if(!isdigit((unsigned char)last)){
    if(!isalpha((unsigned char)last))return false;
    ins=last;
    body=s.substr(0,s.size()-1);
  }
  if(body.empty())return false;
  size_t start=0;
  if(body[0]=='-')start=1;
  if(start>=body.size())return false;
  for(size_t k=start;k<body.size();k++){
    if(!isdigit((unsigned char)body[k]))return false;
  }
  pos=atoi(body.c_str());
  return true;
}

static bool ParseSiteTriple(const vector<string>&fields,SiteOverride&ov,string&err)
{
  if(fields.size()!=3){
    err="expected 3 fields: CHAIN POS[ICODE] N";
    return false;
  }
  if(!ParseChainToken(fields[0],ov.chID)){
    err="invalid chain id '"+fields[0]+"' (use a single character, or '_' for a blank chain id)";
    return false;
  }
  if(!ParsePosIcodeToken(fields[1],ov.pos,ov.ins)){
    err="invalid residue position '"+fields[1]+"' (expected e.g. '15' or '15A')";
    return false;
  }
  if(!ParsePositiveInt(fields[2],ov.topn)){
    err="invalid N '"+fields[2]+"' (must be a positive integer)";
    return false;
  }
  return true;
}

static bool ParseReal(const string&s,double&out)
{
  if(s.empty())return false;
  char* endptr;
  double v=strtod(s.c_str(),&endptr);
  if(endptr!=s.c_str()+s.size())return false;//trailing garbage
  if(endptr==s.c_str())return false;//nothing parsed
  out=v;
  return true;
}

static vector<string> SplitOnChar(const string&s,char sep)
{
  vector<string> fields;
  size_t start=0,pos;
  while((pos=s.find(sep,start))!=string::npos){
    fields.push_back(s.substr(start,pos-start));
    start=pos+1;
  }
  fields.push_back(s.substr(start));
  return fields;
}

int main(int argc,char** argv)
{
  cout<<"###########################################################################"<<endl;
  cout<<"                    FASPR (Version 20200309)                 "<<endl;
  cout<<"  A method for fast and accurate protein side-chain packing, "<<endl;
  cout<<"which is an important problem in protein structure prediction"<<endl;
  cout<<"and protein design."<<endl;
  cout<<endl;
  cout<<"Copyright (c) 2020 Xiaoqiang Huang"<<endl;
  cout<<"Yang Zhang Lab"<<endl;
  cout<<"Dept. of Computational Medicine and Bioinformatics"<<endl;
  cout<<"Medical School"<<endl;
  cout<<"University of Michigan"<<endl;
  cout<<"Email:tommyhuangthu@foxmail.com, xiaoqiah@umich.edu"<<endl;
  cout<<"###########################################################################"<<endl;
  clock_t start,finish;
  float duration;
  start = clock();

  if(argc<2){
    cout<<"Usage: ./FASPR -i input.pdb -o output.pdb\n";
    cout<<"[-s sequence.txt] to load a sequence file\n";
    cout<<"[-n N] to restrict packing to the N highest-probability rotamers per residue\n";
    cout<<"[-n CHAIN:POS[ICODE]:N] (repeatable) to set a per-site rotamer limit, e.g. -n A:15:3\n";
    cout<<"[-m sitemap.txt] a site-map file of 'CHAIN POS[ICODE] N' lines, as an alternative to repeated -n CHAIN:POS:N (mutually exclusive with it)\n";
    cout<<"[-a CHI2MIN:CHI2MAX[:trp]] LimitAromaChi2: only keep PHE/TYR/HIS (and TRP if ':trp' given) rotamers whose\n";
    cout<<"    chi2 (folded into [0,180) by +180 if negative) falls in [CHI2MIN,CHI2MAX], e.g. -a 70:110 or -a 70:110:trp\n";
    return 0;
  }

  /**********************************************************/
  /* get program path and check if rotamer library exists   */
  /**********************************************************/
  char fullpath[2048];
  strcpy(fullpath,argv[0]);
  for(int i = strlen(fullpath); i >= 0; --i){
    if(fullpath[i]=='/' || fullpath[i]=='\\'){
      fullpath[i + 1]='\0';
      break;
    }
  }
  PROGRAM_PATH=(string)fullpath;
  string rotfile=PROGRAM_PATH+"/"+ROTLIB2010;
  fstream infile(rotfile.c_str(),ios::in|ios::binary);
  if(!infile){
    cerr<<"error! cannot find rotamer library "<<ROTLIB2010<<endl;
    exit(0);
  }
  else{
    infile.close();
  }

  string pdbin=(string)"example/1mol.pdb";
  string pdbout=(string)"example/1mol_FASPR.pdb";
  string seqfile=(string)"void";

  bool sflag=false;
  bool haveInlineSiteOverride=false;
  bool haveSiteMapFile=false;
  int i;
  for(i=1;i<argc-1;i++){
   if(argv[i][0]=='-'){
     if(argv[i][1]=='i'){
       i++;
       pdbin=argv[i];
     }
     else if(argv[i][1]=='o'){
       i++;
       pdbout=argv[i];
     }
     else if(argv[i][1]=='s'){
       i++;
       seqfile=argv[i];
       sflag=true;
     }
     else if(argv[i][1]=='n'){
       i++;
       string tok=argv[i];
       if(tok.find(':')!=string::npos){
         SiteOverride ov;
         string err;
         if(!ParseSiteTriple(SplitOnChar(tok,':'),ov,err)){
           cerr<<"error! invalid -n site override '"<<tok<<"': "<<err<<endl;
           return 1;
         }
         SITE_OVERRIDES.push_back(ov);
         haveInlineSiteOverride=true;
       }
       else{
         int nval;
         if(!ParsePositiveInt(tok,nval)){
           cerr<<"error! invalid value for -n: '"<<tok<<"'. -n requires a positive integer.\n";
           return 1;
         }
         ROTAMER_TOPN=nval;
       }
     }
     else if(argv[i][1]=='a'){
       i++;
       string tok=argv[i];
       vector<string> fields=SplitOnChar(tok,':');
       if(fields.size()<2||fields.size()>3){
         cerr<<"error! invalid -a value '"<<tok<<"': expected CHI2MIN:CHI2MAX or CHI2MIN:CHI2MAX:trp\n";
         return 1;
       }
       double minv,maxv;
       if(!ParseReal(fields[0],minv) || !ParseReal(fields[1],maxv)){
         cerr<<"error! invalid -a value '"<<tok<<"': CHI2MIN and CHI2MAX must be real numbers\n";
         return 1;
       }
       if(minv<0. || maxv>180. || minv>=maxv){
         cerr<<"error! invalid -a range '"<<tok<<"': require 0 <= CHI2MIN < CHI2MAX <= 180\n";
         return 1;
       }
       bool trp=false;
       if(fields.size()==3){
         if(fields[2]!="trp"){
           cerr<<"error! invalid -a value '"<<tok<<"': third field must be exactly 'trp'\n";
           return 1;
         }
         trp=true;
       }
       AROMACHI2_ENABLED=true;
       AROMACHI2_MIN=(float)minv;
       AROMACHI2_MAX=(float)maxv;
       AROMACHI2_INCLUDE_TRP=trp;
     }
     else if(argv[i][1]=='m'){
       i++;
       string mapfile=argv[i];
       ifstream mf(mapfile.c_str());
       if(!mf){
         cerr<<"error! cannot open site-map file: "<<mapfile<<endl;
         return 1;
       }
       string line;
       int lineno=0;
       while(getline(mf,line)){
         lineno++;
         size_t h=line.find('#');
         if(h!=string::npos)line=line.substr(0,h);
         size_t b=line.find_first_not_of(" \t\r\n");
         if(b==string::npos)continue;
         istringstream iss(line);
         vector<string> fields;
         string f;
         while(iss>>f)fields.push_back(f);
         SiteOverride ov;
         string err;
         if(!ParseSiteTriple(fields,ov,err)){
           cerr<<"error! invalid line "<<lineno<<" in site-map file '"<<mapfile<<"': "<<err<<endl;
           return 1;
         }
         SITE_OVERRIDES.push_back(ov);
       }
       mf.close();
       haveSiteMapFile=true;
     }
   }
  }
  if(haveInlineSiteOverride && haveSiteMapFile){
    cerr<<"error! cannot combine inline -n CHAIN:POS:N site overrides with a -m site-map file; choose one mechanism.\n";
    return 1;
  }

  Solution faspr;
  faspr.ReadPDB(pdbin);
  if(sflag) faspr.LoadSeq(seqfile);
  else faspr.LoadSeq();
  faspr.BuildSidechain();
  faspr.CalcSelfEnergy();
  faspr.CalcPairEnergy();
  faspr.Search();
  faspr.WritePDB(pdbout);
  finish = clock();
  duration = (float)(finish-start)/CLOCKS_PER_SEC;
  cout<<"#computational time: "<<duration<<" seconds"<<endl;

  return 0;
}
