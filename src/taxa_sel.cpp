// Updated at July 28, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
//version 3.1 or above with _Table_Format

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "taxa_sel.h"
#include "utility.h"
#include "version.h"

#define TLevN 8

using namespace std;

//_Taxa_Table Taxa_table;
_Table_Format_Seq Taxa_table;

string Listfilename;
string Listprefix;

vector <string> Infilename;
vector <string> Sam_name;

//string Listfilename;
string Outfilename = "taxaonomy_selection";

int Level = 5; ////(0: Kingdom, 5: Genus, 6: Species, 7: OTU)
float Max_abd = 0.001;
float Min_abd = 0;
float No_zero_rate = 0.1;
float Ave_t = 0.001;

int Min_seq = 2;

bool Is_cp_correct = true;
bool Is_print = false;

string database = "gg_13";
string database_path;

string Taxa_level[TLevN] = {"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"};

void Print_Help(){
    
    cout << "Taxa-sel version" << Version << endl; 
    cout << "Usage : " << endl;
    cout << "taxa-sel [-option] value" << endl;
    cout << "\toption : " << endl;    
    cout << "\t-l Input list [Required]" << endl;  
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    cout << "\t-o Output file name, default is \"taxaonomy_selection\"" << endl;
    cout << "\t-L Taxonomical level (1-6: Phylum - Species, 7: OTU). default is 5" << endl;
    cout << "\t-r 16s rRNA copy number correction, T(rue) or F(alse), default is T" << endl;
    cout << "\t-P Print distribution barchart, T(rue) or F(alse), default is F" << endl;
    cout << "\t-q Minimum sequence count threshold, default is 2" << endl;
    cout << "\t-m Maximum abundance threshold, default is 0.001 (0.1%)" << endl;
    cout << "\t-n Minimum abundance threshold, default is 0.0 (0%)" << endl;
    cout << "\t-z Minimum No-Zero abundance threshold, default is 0.1 (10%)" << endl;
    cout << "\t-v Minimum average abundance threshold, default is 0.001 (0.1%)" << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    }

int Get_TLevel(string levels){
    switch(levels[0]){
                              case '1':
                              case 'P': 
                              case 'p': return 1; break;
                              case '2':
                              case 'C': 
                              case 'c': return 2; break;
                              case '3':
                              case 'O': 
                              case 'o': return 3; break;
                              case '4':
                              case 'F':
                              case 'f': return 4; break;
                              case '5':
                              case 'G':
                              case 'g': return 5; break;
                              case '6':
                              case 'S':
                              case 's': return 6; break;
                              case '7':
                              case 'U':
                              case 'u': return 7; break;
                              
                              default:
                                      cerr << "Warning: Unrec taxa level: " << levels[0] << endl; return 0; break; 
                              }
            
    return 0;
    }

void Parse_Para(int argc, char * argv[]){
      
      database_path = Check_Env() + "/databases/" + database + "/";
      
      int i = 1;
      
      if (argc ==1)
		Print_Help();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){                            
                            case 'l': Listfilename = argv[i+1]; break; 
                            case 'p': Listprefix = argv[i+1]; break;                             
                            case 'o': Outfilename = argv[i+1]; break;
                            case 'L': Level = Get_TLevel(argv[i+1]); break; 
                            case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_cp_correct = false; break;   
                            case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_print = true; break;   
                            case 'q': Min_seq = atoi(argv[i+1]); break;   
                            case 'm': Max_abd = atof(argv[i+1]); break;
                            case 'n': Min_abd = atof(argv[i+1]); break;
                            case 'z': No_zero_rate = atof(argv[i+1]); break;
                            case 'v': Ave_t = atof(argv[i+1]); break;
                            
                            default : printf("Error: Unrec argument %s\n", argv[i]); Print_Help(); break; 
                            }
         i+=2;
         }         
     
     //check
     if ((Level <= 0) || (Level > 7)){
               cerr << "Warning: Taxonomical level must be 0-7, change to default (5)" << endl;
               Level = 5;
               }
     
     if (Max_abd < Min_abd){
           cerr << "Warning: m must be larger than n, change to default (0, 0.001)" << endl;
           Max_abd = 0.001;
           Min_abd = 0.0;
           }
     
      Level ++; //0: Kingdom
     
     //load
     Load_List(Listfilename.c_str(), Infilename, Sam_name, Listprefix);      
     }

int Load_Taxa(vector <string> infilename, int level){     
    
    hash_map <string, float, std_string_hash> cp_number;
    Load_Copy_Number(database_path, cp_number);   
    
    for (int i = 0; i < infilename.size(); i ++){    
        ifstream infile(infilename[i].c_str(), ifstream::in);
        if (!infile){
                     cerr << "Error: Cannot open input file : " << infilename[i] << endl;
                     continue;
                     }
                      
        int seq_count = 0;
        string buffer;
        getline(infile, buffer);
        while(getline(infile, buffer)){
                              
                              if (buffer.size() == 0) continue;
                                                            
                              stringstream strin(buffer);
                              string id;
                              string otu;
                              string taxa;
                              vector <string> taxa_buffer;
                              string temp;
                              
                              strin >> id >> otu >> temp >> temp;
                              
                              if (level == 8){ //otu
                                        taxa = "OTU_";
                                        taxa += otu;
                                        }
                                         
                              else for (int j = 0; j < level; j++){
                                            strin >> taxa;
                                            
                                            if (taxa.find("Unclassified") != string::npos){
                                                                          taxa = "Unclassified";
                                                                          break;
                                                                          }
                                            if (taxa.find("otu_") != string::npos){
                                                                          taxa = "Unclassified";
                                                                          break;
                                                                          }
                                            
                                            if(taxa[taxa.size()-1] != ';'){
                                                             strin >> temp;
                                                             taxa += "_";
                                                             taxa += temp;                                                                          
                                                             }
                                            taxa_buffer.push_back(taxa);
                                            
                                            //add genus info for species
                                            if (j == 6){ 
                                               taxa = taxa_buffer[5];
                                               taxa = taxa.substr(0, taxa.size()-1);
                                               taxa += "_";
                                               taxa += taxa_buffer[6];
                                               }
                                            
                                            //add prefix for genus
                                            if (j == 5){
                                                  string taxa_prefix = taxa_buffer[j-1];
                                                  taxa_prefix = taxa_prefix.substr(0, taxa_prefix.size() - 1);
                                                  if (taxa_prefix.size() > 3) taxa_prefix = taxa_prefix.substr(0, 3);
                                                  taxa = taxa_prefix + "_" + taxa;
                                                  }
                                            }

                              
                              float cp_num = 1.0;
                              if (Is_cp_correct)
                                 if (cp_number.count(otu) != 0)
                                                     cp_num = cp_number[otu];
                                                                                                                                 
                              Taxa_table.Add_Feature(taxa, i, cp_num);
                              seq_count ++;                             
                                                            
                              }
                        
        infile.close();
        infile.clear();
        
        cout << seq_count << " sequences are loaded" << endl; 
        }
                 
    return 0;
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);     
    
    Taxa_table = _Table_Format_Seq(Sam_name);  
    
    cout << "Taxonomical Selection Starts" << endl << endl;
    cout << "Level is " << Taxa_level[Level-1] << endl;
    cout << "Total Sample Number is " << Infilename.size() << endl;
        
    Load_Taxa(Infilename, Level);
    
    Taxa_table.Filter_Seq_Count(Min_seq);
    
    Taxa_table.Normalization();
    
    Taxa_table.Filter_Abd(Max_abd, Min_abd, No_zero_rate, Ave_t);
    
    cout << "Total Taxa Number is " << Taxa_table.Get_Taxa_Size() << endl;
    //output abd
    cout << "Total Output Taxa Numer is " << Taxa_table.Output_Abd((Outfilename + ".Abd").c_str()) << endl;
    Taxa_table.Output_Count((Outfilename + ".Count").c_str());
    
    if (Is_print){
                  char command[BUFFER_SIZE];
                  sprintf(command, "Rscript %s/Rscript/PM_Distribution.R -i %s -o %s", Check_Env().c_str(), (Outfilename + ".Abd").c_str(), (Outfilename + ".Abd.distribution.pdf").c_str());
                  system(command);
                  }
    
    return 0;
    }
