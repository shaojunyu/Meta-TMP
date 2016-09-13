// Updated at Aug 02, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <unistd.h>

#include "utility.h"
#include "version.h"

#define MAX_BOOT 1000
#define DEF_BOOT 200

using namespace std;

string Listfile;
string Listprefix;

vector <string> Infilename;
vector <string> Outfilename;
vector <string> Sam_name;

string database_path;
string database = "gg_13";

string Outpath = "Rare_Out";

int Seq_depth = 0;
int Bootstrap = DEF_BOOT;

int Mode = -1;

map <string, string> OTU_Taxa;

int printhelp(){
    
    cout << "Random_rare : " << Version << endl;
    cout << "Usage : " << endl;
    cout << "rand-rare [-options] value" << endl;
    cout << "\toption : " << endl;

    cout << "\t-i Input single file name [Conflict with -l]" << endl;
    cout << "\t-l Input filename lis [Conflict with -i]" << endl;
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    cout << "\t-o Output path, default is \"Rare_Out\"" << endl;
    cout << "\t-s Rarefaction depth [Required]" << endl;
    cout << "\t-b Bootstrap for sequence number normalization, default is " << DEF_BOOT << ", maximum is " << MAX_BOOT << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    database_path =  Check_Env() + "/databases/" + database + "/";
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
  
    Mode = -1;
    
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
             case 'i':
                 
                 if (Mode != -1) {
                     
                     cerr << "Error: -i conflicts with -l" << endl;
                     exit(0);
                     
                 }
                 Infilename.push_back(argv[i+1]);
                 Sam_name.push_back("Sample");
                 Mode = 0;
                 break;
                 
             case 'l':
                 
                 if (Mode != -1){
                     
                     cerr << "Error: -l conflicts with -i" << endl;
                     exit(0);
                     
                 }
                 Listfile = argv[i+1];
                 Mode = 1;
                 break;
         
             case 'p': Listprefix = argv[i+1]; break;    
             case 'o': Outpath = argv[i+1]; break;
             case 's': Seq_depth = atoi(argv[i+1]); break;
             case 'b': Bootstrap = atoi(argv[i+1]); break;
                 
             case 'h': printhelp(); break;
             default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
             }
          
         i+=2;
         }
    if (Seq_depth <=0 ){
          cerr << "Please assign the rarefaction sequence depth by -s" << endl;
          exit(0);
          }
    if (Bootstrap <= 0){
          cerr << "Warning: The minimum bootstrap is 1, change to default " << DEF_BOOT << endl;
          Bootstrap = DEF_BOOT;
          }
    
    if (Bootstrap > MAX_BOOT){
        cerr << "Warning: The maximum bootstrap is " << MAX_BOOT << ", change to default " << DEF_BOOT << endl;
        Bootstrap = DEF_BOOT;
    }
    
    Check_Path(Outpath.c_str(), 1);
    
    if (Mode == 1) //list
       Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix);
    
    return 0;
    }

int Load_Taxa(const char * infilename, map <string, string> & table){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Open Taxonomy Annotation file error : " <<infilename << endl;
                 exit(0);
                 }
                              
    unsigned int count = 0;   
                  
    string buffer;
    getline(infile, buffer); //title
    
    while (getline(infile, buffer)){
          
          if (buffer.size() == 0) continue;
          
          string id;
          string taxa;
          
          int first_table =  buffer.find('\t', 0); //location of first '\t'
          int last_table = buffer.rfind('\t', buffer.size()-1);//location of last '\t'
          
          id = buffer.substr(0, first_table);
          taxa = buffer.substr(last_table + 1, buffer.size() - last_table -1 );
          
          if (table.count(id) > 0) {
                               
                               cerr << "Error: Loading Taxonomy Annotation error : Duplicate ID : " << id << endl;
                               exit(0);
                               
                               }
          else (table)[id] = taxa;
          
          count ++;
          
          }
    
    infile.close();
    infile.clear();
    
    return count;
    
    }
/*
void Get_Random(int n, int * order){
    
    srand((int)time(NULL));
    
    int * order_table = new int [n];
    
    for(int i = 0; i< n; i++){
            order_table[i] = i;
            }
    for (int i = 0; i < n; i++ ){
         
         int r =(int)((float) (n-1-i)* rand()/(RAND_MAX+1.0)); 
         int temp = order_table[n-1-i]; //last 
         order_table[n-1-i] = order_table[r]; 
         order_table[r] = temp;
         }

    for (int i = 0; i < n; i++)
        order[i] = order_table[i];
    
    delete [] order_table;
    
    return;
    }
*/
void Get_Random(int n, int * order, int s, int loop){
    
    srand((int)time(NULL) + loop * 3000);
    
    int * order_table = new int [n];
    
    for(int i = 0; i< n; i++){
            order_table[i] = i;
            }
    for (int i = 0; i < s; i++ ){
         
         int r =(int)((float) (n-1-i)* rand()/(RAND_MAX+1.0)); 
         int temp = order_table[n-1-i]; //last 
         order_table[n-1-i] = order_table[r]; 
         order_table[r] = temp;
         }

    for (int i = 0; i < s; i++)
        order[i] = order_table[n-1-i];
    
    delete [] order_table;
    
    return;
    }

int Get_Int(float f){
    int n = (int) f;
    n += (f - n >= 0.5) ? 1 : 0;
    return n;
    }

unsigned int Load(const char * infilename, vector <string> & otus){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open input file: " << infilename << endl;
                 return 0;
                 }
    string buffer;
    getline(infile, buffer); //title
    
    while(getline(infile, buffer)){
                          stringstream strin(buffer);
                          string otu;
                          strin >> otu >> otu;
                          otus.push_back(otu);
                          }
    infile.close();
    infile.clear();
    return otus.size();
    }

int Output(const char * outfilename, map <string, float> abd){
    
     ofstream outfile(outfilename, ofstream::out);
     if (!outfile){
                   cerr << "Error: Cannot open output file: " << outfilename << endl;
                   return 0;                   
                   }
     
     outfile << "#Sequence_Id\tDatabase_Otu\tFLAG\tMAPQ\tClassification" << endl; //title
     
     int count = 0;
     
     for (map <string, float> ::iterator miter = abd.begin(); miter != abd.end(); miter ++){
         float abd = miter->second;
         int seq_num = Get_Int(abd);
         
         if (seq_num <=0) continue;
         
         string taxa;
         if (OTU_Taxa.count(miter->first) != 0)
             taxa = OTU_Taxa[miter->first];
         else{
             cerr << "Warning: Cannot find OTU " << miter->first << endl;
             taxa = "Unclassified";
            }
         char seq_output[1000];
         
         for (int i = 0; i < seq_num; i ++){
             sprintf(seq_output, "Rand_seq_%d\t%s\t0\t1\t%s", count, (miter->first).c_str(), taxa.c_str());
             outfile << seq_output << endl;
             count ++;
             }
         }
         
     outfile.close();
     outfile.clear();
     return count;
     }

int Output(string outpath, string samname, int mode, map <string, float> abd){
    
    string outfilename = outpath + "/";
    
    if (mode == 0)
            outfilename += "classification.rare.txt"; //Single mode
    else{
        Check_Path((outpath + "/" + samname).c_str(), 1);
        outfilename = outfilename + "/" + samname + "/" + "classification.rare.txt";
        }
    return Output(outfilename.c_str(), abd);
    }

void Rare_Single_Boot(vector <string> otus, map <string, float> & abd, int s, int b){
    
    if (s >= otus.size()){
          cerr << "Warning: Rarefaction sequence depth must be larger than input sequence number" << endl;
          return;
          } 
    
    for (int i = 0; i < b; i ++){
    
        int * order  = new int [s];
        Get_Random(otus.size(), order, s, i);
        for (int j = 0; j < s; j ++){
            if (abd.count(otus[order[j]]) == 0)
               abd[otus[order[j]]] == 0;
             abd[otus[order[j]]] ++;
             }
        }
    
    for (map <string, float> :: iterator miter = abd.begin(); miter != abd.end(); miter ++){
        miter->second /= (float) b;
        miter->second = (float) Get_Int(miter->second);
        }
    
    //correct abd
    float sum = 0;
    for (map <string, float> :: iterator miter = abd.begin(); miter != abd.end(); miter ++)
        sum += miter->second;
    
    float rate = (float) s / sum;
    //cout << rate << endl; //debug
    
    for (map <string, float> :: iterator miter = abd.begin(); miter != abd.end(); miter ++)
        miter->second *= rate;
    }

int main(int argc, char * argv[]){        
    
    Parse_Para(argc, argv);
    
    Load_Taxa((database_path + "/taxonomy_annotation.txt").c_str(), OTU_Taxa);
    
    int count_out = 0;
    
    for (int i = 0; i < Infilename.size(); i ++){
        vector <string> otus;
        map <string, float> abd;
        int n = Load(Infilename[i].c_str(), otus);
        cout << n << " input seqs loaded" << endl;
        if (Seq_depth >= n){
            cerr << "Warning: Sample " << i + 1 << ": Rarefaction sequence depth must be larger than input sequence number, removed from the output" << endl;
            continue;
            }
        Rare_Single_Boot(otus, abd, Seq_depth, Bootstrap);
        cout << Output(Outpath, Sam_name[i], Mode, abd) << " random seqs output" << endl;
        count_out ++;
        }
    cout << Infilename.size() << " Sample(s) input" << endl;
    cout << count_out << " Sample(s) output" << endl;
    return 0;
    }
