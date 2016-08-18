// Updated at May 13, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
//version 3.1 or above with _Table_Format

#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <unistd.h>

#include "utility.h"
#include "hash.h"
#include "table_format.h"
#include "version.h"

using namespace std;

string database = "gg_13";
string database_path;

hash_map <string, double, std_string_hash> Nsti_table;
hash_map <string, float, std_string_hash> Cp_number;

string Listfilename;
string Listprefix;
vector <string> Infilename;
vector <string> Sam_name; 

string Tablefilename;

string Outfilename = "NSTI.out";
bool Is_abd = true;

string Otu_nsti_file = "otu_nsti.tab";

int Coren = 0;

int Mode = -1; //0: single 1: multi 2: multi_table

int printhelp(){
    
    cout << "Class-func-NSTI version " << Version << endl;
    cout << "Usage : class-func-nsti [Options] Value" << endl;
    cout << "Options : " << endl;
    cout << "\t-i Input single filename [Conflict with -l and -T]" << endl;
    //cout << "\t-a Optional ID file for multi-samples [Optional]" << endl;
    cout << "\t-l Input filename list [Conflict with -i and -T]" << endl;
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    cout << "\t-T Input OTU table format (*.Abd) [Conflict with -i and -T]" << endl;
    cout << "\t-o Output file, default is \"NSTI.out\"" << endl;
    cout << "\t-t Cpu core number, default is auto" << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){
    
    database_path =  Check_Env() + "/databases/" + database + "/";
    Otu_nsti_file = database_path + Otu_nsti_file;
    
    Listprefix = "";
    
    if (argc ==1) 
		printhelp();
    
    int i = 1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'i':                                  
                                 if (Mode != -1) {
                                             
                                             cerr << "Error: -i conflicts with -l and -T" << endl;
                                             exit(0);
                                             
                                             }                                  
                                  Infilename.push_back(argv[i+1]);
                                  Sam_name.push_back("Sample");
                                  Mode = 0;                                  
                                  break;
                                                                                  
                            case 'l':                                       
                                      if (Mode != -1){                                                  
                                                  cerr << "Error: -l conflicts with -i and -T" << endl;
                                                  exit(0);                                                  
                                                  }
                                     
                                      Listfilename = argv[i+1];
                                      Mode = 1; 
                                      break;
                            
                            case 'p': Listprefix = argv[i+1]; break;
                                         
                            case 'T': 
                                 if (Mode != -1){                                                  
                                                  cerr << "Error: -T conflicts with -i and -l" << endl;
                                                  exit(0);                                                  
                                                  }
                                 Tablefilename = argv[i+1]; 
                                 Mode = 2; 
                                 break;
                                 
                            case 'r': if ((argv[i+1][0] == 'f') || (argv[i+1][0] == 'F')) Is_abd = false; break;                                              
                            case 'o': Outfilename = argv[i+1]; break;    
                            case 't': Coren = atoi(argv[i+1]); break;                                                                                                           
                            case 'h': printhelp(); break;

                            default : printf("Error: Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }

    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren = max_core_number;
                    }     
    return 0;
    }

int Load_Nsti(const char * nsti_table_name, hash_map <string, double, std_string_hash> & nsti_table){
    
    ifstream infile(nsti_table_name, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open NSTI table file: " << nsti_table_name << endl;
                 return 0;
                 }
    int count;
    string buffer;
    while(getline(infile, buffer)){
                          if (buffer[0] == '#') continue;
                          stringstream strin(buffer);
                          string otu;
                          double otu_nsti_value;
                          strin >> otu >> otu_nsti_value;
                          nsti_table[otu] = otu_nsti_value;
                          count ++;
                          }
    infile.close();
    infile.clear();
    return nsti_table.size();
    }

double Get_Nsti(const char * infilename){
      
      hash_set <string, std_string_hash> Otu_set;
      
      ifstream infile(infilename, ifstream::in);
      if (!infile){
                   cerr << "Error: Cannot open input file: " << infilename << endl;
                   return -1;
                   }
      
      double count = 0;
      string buffer;
      
      double nsti_value = 0;
      
      getline(infile, buffer); //label
      while(getline(infile, buffer)){
                            
                            if (buffer[0] == '#') continue;
                            stringstream strin(buffer);
                            string otu;
                            strin >> otu >> otu;
                            
                            double cp_number = 1.0;
                            
                            if ((Nsti_table).count(otu) == 0) continue;                            
                            
                            if (!Is_abd) 
                               if (Otu_set.count(otu) != 0) continue; //ignore abd
                            
                            if (Cp_number.count(otu) != 0)
                                                     cp_number = Cp_number[otu];
                            
                            //debug
                            //cout << cp_number << endl;
                                                        
                            count += (1.0 / cp_number);
                            nsti_value += (Nsti_table[otu] / cp_number);
                            Otu_set.insert(otu);                                                                                    
                            }
                 
      infile.close();
      infile.clear();
               
      return nsti_value /  count;
      }

double Get_Nsti(vector <string> OTU, vector <float> Abd){
      
      double count = 0;
      double nsti_value = 0;
      
      for (int i = 0; i < OTU.size(); i ++){
                            
                            if (Abd[i] == 0.0) continue;
                            
                            float abd = Abd[i];
                            string otu = OTU[i].substr(4, OTU[i].size() - 4);                            
                            
                            if ((Nsti_table).count(otu) == 0) continue;                            
                            
                            if (!Is_abd) abd = 1.0;
                                                        
                            count += abd;
                            nsti_value += Nsti_table[otu] * abd;
                            }
             
      return nsti_value /  count;
      }

void Multi_Nsti(){
     
    Load_List(Listfilename.c_str(), Infilename, Sam_name, Listprefix); 
    
    ofstream outfile(Outfilename.c_str(), ofstream::out);
    if (!outfile){
                  cerr << "Error: Cannot open output file: " << Outfilename << endl;
                  return;
                  }
    
    double * nsti_values = new double [Infilename.size()];
    
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < Infilename.size(); i ++)
        nsti_values[i] = Get_Nsti(Infilename[i].c_str());
    
    outfile << "Sample\tNSTI_Value"<< endl; //title
    
    for (int i = 0; i < Infilename.size(); i ++)
        outfile << Sam_name[i] << "\t" << nsti_values[i] << endl;
        
    outfile.close();
    outfile.clear();
    }

void Single_Nsti(){
    
    cout << "The NSTI value is " << Get_Nsti(Infilename[0].c_str()) << endl;
    
    }

void Multi_Nsti_Table(){
     _Table_Format table(Tablefilename.c_str());
     
     ofstream outfile(Outfilename.c_str(), ofstream::out);
    if (!outfile){
                  cerr << "Error: Cannot open output file: " << Outfilename << endl;
                  return;
                  }
    
    double * nsti_values = new double [table.Get_Sample_Size()];
    vector <string> OTUs = table.Get_Feature_Names();
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < table.Get_Sample_Size(); i ++){
        vector <float> abd = table.Get_Abd(i);
        nsti_values[i] = Get_Nsti(OTUs, abd);
        }
        
    outfile << "Sample\tNSTI_Value"<< endl; //title
    
    for (int i = 0; i < table.Get_Sample_Size(); i ++)
        outfile << (table.Get_Sample_Names())[i] << "\t" << nsti_values[i] << endl;
        
    outfile.close();
    outfile.clear();
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    Load_Nsti(Otu_nsti_file.c_str(), Nsti_table);
    Load_Copy_Number(database_path, Cp_number); 
    
    switch (Mode){
           case 0: Single_Nsti(); break;
           case 1: Multi_Nsti(); break;
           case 2: Multi_Nsti_Table(); break;
           }
    return 0;
    }

