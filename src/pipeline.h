// Updated at July 7, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>

#include <sys/stat.h>
#include <omp.h>

#include "utility.h"
#include "version.h"

#define BUFFER_SIZE 5000
#define TLevN 7
#define FLevN 4

#define MAX_BOOT 1000
#define DEF_BOOT 200

#define MIN_SAM_NUM 3

#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Min(a,b) (((a) < (b)) ? (a) : (b))

#ifndef PIPELINE_H    
#define PIPELINE_H

string Bin_path;
string R_path;

string TLevel[TLevN] = {"phylum", "class", "order", "family", "genus", "species", "OTU"};
string FLevel[FLevN] = {"l1", "l2", "l3", "KO"};

bool TLevel_Set[TLevN] = {false};
bool FLevel_Set[FLevN] = {false};

//list
//string Id_file;
string Meta_file;
string Seq_list_file;
//string Group_file;
string Taxa_list_file;
string Func_list_file;
string Report_file;
string Error_file;

string List_prefix;

//Folder
string Abd_dir = "Abundance_Tables";
string Dist_dir = "Distance_Matrix";
string Clust_dir = "Clustering";
string Marker_dir = "Markers";
string Network_dir = "Network";
string Alpha_dir = "A_Diversity";
string Beta_dir = "B_Diversity";
string Sampleview_dir = "Sample_Views";
string Singlesample_dir = "Single_Sample";
string Singlesamplerare_dir = "Single_Sample.Rare";
string Singlesamplelist_dir = "Single_Sample.List";
string Temp_dir = "Temp";

//para
string Out_path;
char Seq_type = 'r';
string Align_mode = "very-sensitive-local";
string Paired_mode = "fr";
int Length_t = 0;
int Cluster = 2;
float Network_t = 0.5;
char Is_pair = 'F';
char Is_format_check = 'F';

bool Is_taxa = true;
bool Is_func = true;
bool Is_rare = false;
bool Is_rare_curve = false;
bool Is_paired_seq = false;

int Rare_depth = 0;
int Bootstrap = DEF_BOOT;

//step
int Step = -1;

int Coren = 0;

vector <string> Seq_files;
vector <string> Ids;

//Sample info
int Input_sam_num = 0;
int Filter_sam_num = 0;

int printhelp(){
    
    cout << "Welcome to Parallel-META Pipeline" << endl;
    cout << "Version : " << Version << endl;
    cout << "Usage : " << endl;
    cout << "pipeline [Options] Value" << endl;
    cout << "Option : " << endl << endl;
    
    cout << "\tInput options:" << endl;
    cout << "\t  -i Sequence list file, pair-ended sequences are supported [Conflict with '-l']" << endl;
    cout << "\t  -m Meta data file [Required]" << endl;    
    cout << "\t  -l Taxonomic analysis results list [Conflict with '-i']" << endl;
    cout << "\t  -p List file path prefix for '-i' and '-l' [Optional]" << endl;
    cout << endl;
    
    cout << "\tOutput options:" << endl;
    cout << "\t  -o Output path, default is \"default_out\"" << endl;
    cout << endl;
    
    cout << "\tProfiling parameters:" << endl;
    cout << "\t  -M Sequence type, T (Shotgun) or F (16S rRNA), default is F" << endl;
   	cout << "\t  -e Alignment mode, 0: very fast, 1: fast, 2: sensitive, 3: very-sensitive, default is 3" << endl;
   	cout << "\t  -P Pair-end sequence orientation, 0: Fwd & Rev, 1: Fwd & Fwd, 2: Rev & Fwd, default is 0" << endl;
    cout << "\t  -r rRNA length threshold of rRNA extraction. 0 is disabled, default is 0" << endl;
    cout << "\t  -k Sequence format check, T(rue) or F(alse), default is F" << endl;
    cout << "\t  -f Functional analysis, T(rue) or F(alse), default is T" << endl;
    cout << endl;
    
    cout << "\tStatistic parameters:" << endl;
    cout << "\t  -L Taxonomical levels (1-6: Phylum - Species). Multiple levels are accepted" << endl;
    cout << "\t  -F Functional levels (Level 1, 2, 3 or 4 (KO number)). Multiple levels are accepted" << endl;
    cout << "\t  -s Sequence number normalization depth, 0 is disabled, default is disable" << endl;
    cout << "\t  -b Bootstrap for sequence number normalization, default is " << DEF_BOOT << ", maximum is " << MAX_BOOT << endl;
    cout << "\t  -R Rarefaction curve, T(rue) or F(alse), default is F" << endl;
    cout << "\t  -E If the samples are paired, T(rue) or F(alse), default is F" << endl;
    cout << "\t  -c Cluter number, default is 2" << endl;
    cout << "\t  -T Network analysis edge threshold, default is 0.5" << endl;
    cout << endl;
    
    cout << "\tOther options:" << endl;
    cout << "\t  -t cpu core number, default is auto" << endl;
    cout << "\t  -h help" << endl;
    
    exit(0);
    
    return 0;
    
    };

void Print_Report(const char * outfilename){
     
     ofstream outfile(outfilename, ios::out);
     if(!outfile){
                  cerr << "Error: Cannot open report file: " << outfilename << endl;
                  return;
                  }
     
     outfile << "Parallel-META Pipeline Analysis Report" << endl;
     outfile << "Version " << Version << endl;
     
     if (Step == 0){ //Profiling info
              outfile << "Input Type : Sequence list" << endl;
              outfile << "Sequence List : " << Seq_list_file << endl; 
              outfile << "Sequence Type : ";
              if (Seq_type == 'm') outfile << "Metagenomic Shotgun Sequences" << endl;
              else outfile << "rRNA Targeted Sequences" << endl; 
              
              outfile << "Pair-end Sequences : ";
              if (Is_paired_seq) outfile << "Yes" << endl;
              else outfile << "No" << endl;                                                        
              }
     else {
          outfile << "Input Type : Taxonomic analysis results list" << endl;
          outfile << "Results List : " << Taxa_list_file << endl;
          }
     
     if (Is_rare){                  
        outfile << "Sequence Normalization Depth: " << Rare_depth << endl;
        outfile << "Sequence Normalization Bootstrap : " << Bootstrap << endl;         
        }
     
     outfile << "Functional Anlaysis : ";
     if (Is_func) outfile << "Yes" << endl;
     else outfile << "No" << endl;
               
     outfile << "Number of Input Sample : " << Input_sam_num << endl;
     outfile << "Number of Output Sample : " << Filter_sam_num << endl;
     
     outfile.close();
     outfile.clear();
     }

int Parse_TLevel(string levels){
    
    int count = 0;
    
    for (int i = 0; i < TLevN; i++)
        TLevel_Set[i] = false;
        
    for(int i = 0; i < levels.size(); i ++){
            
            switch(levels[i]){
                              case '1':
                              case 'P': 
                              case 'p': TLevel_Set[0] = true; count ++; break;
                              case '2':
                              case 'C': 
                              case 'c': TLevel_Set[1] = true; count ++; break;
                              case '3':
                              case 'O': 
                              case 'o': TLevel_Set[2] = true; count ++; break;
                              case '4':
                              case 'F':
                              case 'f': TLevel_Set[3] = true; count ++; break;
                              case '5':
                              case 'G':
                              case 'g': TLevel_Set[4] = true; count ++; break;
                              case '6':
                              case 'S':
                              case 's': TLevel_Set[5] = true; count ++; break;
                              //case '7':
                              //case 'U':
                              //case 'u': TLevel_Set[6] = true; count ++; break;
                              
                              default:
                                      cerr << "Warning: Unrec taxa level: " << levels[i] << endl; break; 
                              }
            
            }
    
    if (count <= 0){
              cerr << "Error: At least select 1 taxa level by '-L'" << endl;
              exit(0);
              }
    return count;
    }

int Parse_FLevel(string levels){
    
    int count = 0;
    
    for (int i = 0; i < FLevN; i++)
        FLevel_Set[i] = false;
    
    for(int i = 0; i < levels.size(); i ++){
            
            switch(levels[i]){
                              case '1': FLevel_Set[0] = true; count ++; break;
                              case '2': FLevel_Set[1] = true; count ++; break;
                              case '3': FLevel_Set[2] = true; count ++; break;
                              case '4': FLevel_Set[3] = true; count ++; break;
                              
                              default:
                                      cerr << "Warning: Unrec func level: " << levels[i] << endl; break; 
                              }
            
            }
    
    if (count <= 0){
              cerr << "Error: At least select 1 func level by '-F'" << endl;
              exit(0);
              }
    return count;
    }

int Parse_Para(int argc, char * argv[]){

    Bin_path = Check_Env() + "/bin/";
    R_path = Check_Env() + "/Rscript";
    
    List_prefix = "";
    
    Step = -1;
    Cluster = 2;
    Seq_type = 'r';
    Coren = 0;
    Out_path = "default_out/";
    
    //init TLevel & FLevel
    TLevel_Set[0] = true; //phylum
    TLevel_Set[4] = true; //genus
    //TLevel_Set[6] = true; //OTU
    
    FLevel_Set[1] = true; //l2
    FLevel_Set[2] = true; //l3
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            //basic args
                            case 'm': Meta_file = argv[i+1]; break;
                            case 'i': Seq_list_file = argv[i+1]; 
                                      if (Step == 1) {
                                                     cerr << "Error: -i conflicts with -l and -f" << endl;
                                                     exit(0);
                                                     }
                                      else Step = 0;
                                      break;
                            //case 'g': Group_file = argv[i+1]; break;
                            //case 'a': Id_file = argv[i+1]; break;        
                            case 'l': Taxa_list_file = argv[i+1]; 
                                      if (Step == 0) {
                                                     cerr << "Error: -l conflicts with -i" << endl;
                                                     exit(0);
                                                     }
                                      else Step = 1;
                                      break;
                            
                            case 'p': List_prefix = argv[i+1]; break;
                                      
                            case 'o': Out_path = argv[i+1]; break;
                            
                            case 'e' : switch (argv[i+1][0]){
                                                             case '0': Align_mode = "very-fast-local"; break;
                                                             case '1': Align_mode = "fast-local"; break;
                                                             case '2': Align_mode = "sensitive-local"; break;
                                                             case '3': Align_mode = "very-sensitive-local"; break;
                                                             default: Align_mode = "very-sensitive-local"; break;                   
                                              }
                                        break; //Default is "very-sensitive-local"
                            case 'P': switch (argv[i+1][0]){
                                                             case '0': Paired_mode = "fr"; break;
                                                             case '1': Paired_mode = "ff"; break;
                                                             case '2': Paired_mode = "rf"; break;
                                                             default: Paired_mode = "fr"; break;
                                                             }
                                       break; // Default is "fr"
                                       
                            case 'M': if ((argv[i+1][0] == 'T') || (argv[i+1][0]) == 't' ) Seq_type = 'm'; break;
                            case 'r': Length_t = atoi(argv[i+1]); break; 
                            case 'k': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T' )) Is_format_check = 'T'; 
                                      else if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F' )) Is_format_check = 'F';
                                      break; 
                            
                            //adv args                                                                                                                  
                            case 'f': if ((argv[i+1][0] == 'f') || (argv[i+1][0]) == 'F' ) Is_func = false; break;
                            
                            case 'L': Parse_TLevel(argv[i+1]); break;
                            case 'F': Parse_FLevel(argv[i+1]); break;                                                                
                            
                            case 'R': if ((argv[i+1][0] == 't') || (argv[i+1][0]) == 'T' ) Is_rare_curve = true; break;
                            case 's': Rare_depth = atoi(argv[i+1]); Is_rare = true; break;
                            case 'b': Bootstrap = atoi(argv[i+1]); break;
                            case 'E': if ((argv[i+1][0] == 'T') || (argv[i+1][0]) == 't' ) Is_pair = 'T'; break;
                            case 'c': Cluster = atoi(argv[i+1]); break;
                            case 'T': Network_t = atof(argv[i+1]); break;
                            
                            //other args
                            case 't': Coren = atoi(argv[i+1]); break;         
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    
    if ((Step > 2) || (Step < 0)){
              cerr << "Warning: Step (-s) must be 0-1, change to default (0)" << endl;
              Step = 0; 
              }
    
    Abd_dir = Out_path + "/Abundance_Tables";
    Dist_dir = Out_path + "/Distance_Matrix";
    Clust_dir = Out_path + "/Clustering";
    Marker_dir = Out_path + "/Markers";
    Network_dir = Out_path + "/Network";
    Alpha_dir = Out_path + "/Alpha_Diversity";
    Beta_dir = Out_path + "/Beta_Diversity";
    Sampleview_dir = Out_path + "/Sample_Views";
    Singlesample_dir = Out_path + "/Single_Sample";
    Singlesamplerare_dir = Out_path + "/Single_Sample.Rare";
    Singlesamplelist_dir = Out_path + "/Single_Sample.List";
    Temp_dir = Out_path + "/Temp";
    
    Report_file = Out_path + "/Analysis_Report.txt";
    Error_file = Out_path + "/error.log";
    remove(Error_file.c_str());
    
    int Max_Core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > Max_Core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren = Max_Core_number;
                    }
    if (Rare_depth <= 0) Is_rare = false;
    
    if (Bootstrap <= 0){
        cerr << "Warning: The minimum bootstrap is 1, change to default " << DEF_BOOT << endl;
        Bootstrap = DEF_BOOT;
    }
    if (Bootstrap > MAX_BOOT) {
        cerr << "Warning: The maximum bootstrap is " << MAX_BOOT << ", change to default " << DEF_BOOT << endl;
        Bootstrap = DEF_BOOT;
        }
    }  

void Run_With_Error(char * command, const char * error){
     
     string command_with_error;
     command_with_error = command;
     command_with_error +=" 2>>";
     command_with_error += error;
     system(command_with_error.c_str());
     }

void Run_With_Error(const char * command, const char * program, const char * error){
     
     string command_with_error;
     //make error title
     command_with_error = "echo \"";
     command_with_error += program;
     command_with_error += ":\" >>";
     command_with_error += error;
     system(command_with_error.c_str());
     
     //make command
     command_with_error = command;
     command_with_error +=" 2>>";
     command_with_error += error;
     system(command_with_error.c_str());
     }

bool Check_Ids(vector <string> ids){
     
     if (ids.size() ==0) return false;
     for (int i = 0; i < ids.size(); i ++)
         if ((ids[i][0] >= '0') && (ids[i][0] <= '9')) //started by number
            return false;
     return true;
     }

void Check_Ids_By_Path(const char * path, vector <string> & ids){
    //check ids
    vector <string> newids;
    newids.push_back(ids[0]); //title
    
    for (int i = 1; i < ids.size(); i ++){
        string a_path = path;
        a_path += "/";
        a_path += ids[i];
        if (Check_Path(a_path.c_str())){
            newids.push_back(ids[i]);
        }
    }
    
    if (ids.size() != newids.size())
        ids = newids;

    }

bool Check_Metadata_By_Ids(vector <string> ids, string & metafile, string outmetafile){
    
    //update metadata
    ifstream infile(metafile.c_str(), ifstream::in);
    if (!infile){
        cerr << "Error: Cannot open input metadata file: " << metafile << endl;
        return false;
    }
    
    ofstream outfile(outmetafile.c_str(), ofstream::out);
    if (!outfile){
        cerr << "Error: Cannot open output metadata file: " << outmetafile << endl;
        return false;
    }

    map <string, string> metadata;
    string buffer;
    getline(infile, buffer); //title
    outfile << buffer << endl; //title
    while(getline(infile, buffer)){
        stringstream strin(buffer);
        string id;
        strin >> id;
        metadata[id] = buffer;
        }
    
    for (int i = 1; i < ids.size(); i ++){
        
        if (metadata.count(ids[i]))
            outfile << metadata[ids[i]] << endl;
        else{
            cerr << "Error: Cannot find metadata for sample: " << ids[i] << endl;
            return false;
            }
        }
    
    infile.close();
    infile.clear();
    outfile.close();
    outfile.clear();
    
    metafile = outmetafile;
    
    return true;
    }

#endif
