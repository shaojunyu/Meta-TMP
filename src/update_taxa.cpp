//Update Taxonomy of GG
//For PM 3 or above
//For GG_13-8
//Updatd at Aug 01, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>
#include "utility.h"
#include "version.h"
#include "hash.h"

using namespace std;

string database = "gg_13";
string database_path;

string Listfilename;
string Listprefix;
vector <string> Infilelist;

int Mode = -1; //0: single; 1: list

int printhelp(){

    cout << "Update-taxa: " << Version << endl;
    cout << "Usage : update-taxa [Options] Value" << endl;
    cout << "\toption : " << endl;
    cout << "\t-i Input filename [Conflict with -l]" << endl;
    cout << "\t-l Input filename list [Conflict with -i]" << endl;
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    cout << "\t-h Help" << endl;

    exit(0);

    return 0;

    };

void Parse_Para(int argc, char * argv[]){
    
    database_path = Check_Env() + "/databases/" + database + "/"; 
     
    if(argc==1){
                printhelp();
                }

    int i=1;
    while(i<argc){
                 if(argv[i][0]!='-'){
                                     cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                                     exit(0);
                                     }
                 switch(argv[i][1]){
                                    case 'i' : if (Mode == 1){
                                                        cerr << "Error: -i conflists with -l" << endl;
                                                        exit(0);
                                                        } 
                                               else {                                                    
                                                    Infilelist.push_back(argv[i+1]);
                                                    //Infilename = argv[i+1]; 
                                                    Mode = 0;
                                                    }                                               
                                               break;
                                               
                                    case 'l' : if (Mode == 0){
                                                        cerr << "Error: -l conflists with -i" << endl;
                                                        exit(0);
                                                        } 
                                               else {
                                                    Listfilename = argv[i+1];
                                                    Mode = 1;
                                                    }
                                               break;
                                    case 'p' : Listprefix = argv[i+1]; break;
                                    case 'h' : printhelp(); break;           
                                    default : cerr << "Error: Unrec argument " << argv[i] << endl; break;
					      printhelp();
					      break;
                                    }
		i=i+2;
                  
                  }
    
    if (Mode == 1)
       Load_List(Listfilename.c_str(), Infilelist, Listprefix); 
    } 



unsigned int Read_Taxonomy(const char * infilename, hash_map<string, string, std_string_hash> & table){
    
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
    
unsigned int Parse_Taxonomy(const char * infilename, const char * outfilename, hash_map<string, string, std_string_hash> * table){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Open Mapping File error : " << infilename << endl;
                 exit(0);
                 }

    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                  cerr <<"Error: Open output file error : " << outfilename << endl;
                  exit(0);
                  }
    
    outfile << "#Sequence_Id\tDatabase_OTU\tFLAG\tMAPQ\tClassification" << endl;
    
    unsigned int match_count = 0;
    unsigned int total_count = 0;
    
    string buffer;
        
    while(getline(infile, buffer)){
                          if (buffer.size()==0) continue;
                          stringstream strin(buffer);
                          string seq_id, database_id, flag;
                          int mapq;
                          strin >> seq_id >> database_id >> flag >> mapq;
                          
                          if (table->count(database_id) != 0){        
                          //Output the taxonomy annotation                                             
                             outfile << seq_id << "\t" << database_id << "\t" << flag << "\t" << mapq << "\t";
                             string taxa_buffer = (*table)[database_id];
                             outfile << taxa_buffer << endl;       
                             match_count ++;
                             }                            
                                                                                                                                                                                                                    
                          total_count ++;                                                    
                          }
    
    infile.close();
    infile.clear();
    outfile.close();
    outfile.clear();
    
    //cout << total_count << " sequences are loaded " << endl;
    
    return match_count;
    
    }

int main(int argc, char * argv[]){
            
    Parse_Para(argc, argv);
    string taxonomyfilename = database_path + "taxonomy_annotation.txt"; 
    
    hash_map<string, string, std_string_hash> table;
    
    Read_Taxonomy(taxonomyfilename.c_str(), table);
    
    for (int i = 0; i < Infilelist.size(); i ++){
         //copy the inputfile
         string command = "cp "+Infilelist[i] + " " + Infilelist[i] + ".bk";
	     system(command.c_str());
	     cout << Parse_Taxonomy((Infilelist[i] + ".bk").c_str(), Infilelist[i].c_str(), &table) << " sequence output" << endl;
        }
    
    return 0;
    }
