// Updated at Nov 25, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <string.h>

#include "version.h"

using namespace std;

string Infilename;
string Metafilename;
string Outfilename;

map <string, vector <string> > Data;
map <string, int> Id_to_order;
string Key;

int printhelp(){
    
    cout << "Split-Matrix Version : " << Version << endl;
    cout << "Usage : " << endl;
    cout << "split_matrix [-options] value" << endl;
    cout << "\toption : " << endl;

    cout << "\t-i Input matrix file [Required]" << endl;
    cout << "\t-m Meta data file [Required]" << endl;
    cout << "\t-c Group key [Required]" << endl;
    cout << "\t-o Result output file, default is \"Out\"" << endl;

    cout << "\t-h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    Outfilename = "Out";
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error: Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'i': Infilename = argv[i+1]; break;
                            case 'm': Metafilename = argv[i+1]; break;
                            case 'o': Outfilename = argv[i+1]; break;
                            case 'c': Key = argv[i+1]; break;
                                                                                       
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    
    if (Key.size() == 0){
                   
                   cerr << "Error: Please input the group key by -c" << endl;
                   exit(0);
                   }
    }

int Load_Meta(const char * infilename, string key){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << infilename << endl;
                 return 0;
                 }
    
    string buffer;
    getline(infile, buffer);
    stringstream strin(buffer);
    string key_temp;
    int key_index = 0;
    bool is_key = false;
    while(strin >> key_temp){
                if (key_temp == key){
                             is_key = true;
                             break;
                             }
                key_index ++;
                }
    if (!is_key){
                 cerr << "Error: Cannot find meta-data : " << key << endl;
                 return 0;
                 }
    
    int count = 0;
    while(getline(infile, buffer)){
                          if (buffer.size() == 0) continue;
                          strin.clear();
                          strin.str(buffer);
                          string id;
                          string group;
                          strin >> id;
                          for (int i = 0; i < key_index; i ++)
                              strin >> group;
                          Data[group].push_back(id);
                          count ++;
                          } 
    
    infile.close();
    infile.clear();
    
    return count;
    }

int Load_Matrix(const char * infilename, float * matrix, int n){
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << infilename << endl;
                 return 0;
                 }
    
    string buffer;
    getline(infile, buffer);
    int count  = 0;
    
    for (int i = 0; i < n; i ++){
                          getline(infile, buffer);                         
                          stringstream strin(buffer);
                          string id;
                          strin >> id;
                          Id_to_order[id] = i;
                          float temp = 0;
                          for (int j = 0; j < n; j ++){
                              strin >> temp;
                              matrix[i * n + j] = temp;                                                            
                              }
                                                                              
                         
                          } 
                          
    infile.close();
    infile.clear();
    return 0;
    }

int Output_Data(const char * outfilename, float * matrix, int n){
    
    for (map<string, vector <string> > :: iterator miter = Data.begin(); miter != Data.end(); miter ++){
        
        string outname = outfilename;
        outname += ".";
        outname += miter->first;
        
        ofstream outfile(outname.c_str(), ofstream::out);
        if (!outfile){
                      cerr << "Error: Cannot open file : " << outname << endl;
                      return 0;
                      }
        
        
        for(int i = 0; i < (miter->second).size(); i ++)
                outfile << "\t" << (miter->second)[i];
        outfile << endl;
        for(int i = 0; i < (miter->second).size(); i ++){
                outfile << (miter->second)[i];
                for(int j = 0; j < (miter->second).size(); j ++)
                outfile << "\t" << matrix[Id_to_order[(miter->second)[i]] * n + Id_to_order[(miter->second)[j]]];
                outfile << endl;
                }
        
        outfile.close();
        outfile.clear();
        }
    return Data.size();
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    int file_count = Load_Meta(Metafilename.c_str(), Key);
    float * sim_matrix = new float [file_count * file_count];
    memset(sim_matrix, 0, file_count * file_count * sizeof(float));
    cout << file_count << " Samples in total" << endl;
    
    Load_Matrix(Infilename.c_str(), sim_matrix, file_count);
    cout << Output_Data(Outfilename.c_str(), sim_matrix, file_count) << " Groups output" << endl;
    return 0;
    }
