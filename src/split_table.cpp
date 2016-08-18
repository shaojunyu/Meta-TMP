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

map <string, string> Id_to_group;
map <string, vector <string> > Data;
string Title;
string Key;

int printhelp(){
    
    cout << "Split-Table Version : " << Version << endl;
    cout << "Usage : " << endl;
    cout << "split-table [-options] value" << endl;
    cout << "\toption : " << endl;

    cout << "\t-T Feature table (*.Abd or *.Count) file [Required]" << endl;
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
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
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
    
    while(getline(infile, buffer)){
                          if (buffer.size() == 0) continue;
                          strin.clear();
                          strin.str(buffer);
                          string id;
                          string group;
                          strin >> id;
                          for (int i = 0; i < key_index; i ++)
                              strin >> group;
                          Id_to_group[id] = group;
                          } 
    
    infile.close();
    infile.clear();
    
    return Id_to_group.size();
    }

int Load_Data(const char * infilename){
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << infilename << endl;
                 return 0;
                 }
    
    string buffer;
    getline(infile, Title);
    int count  = 0;
    
    while(getline(infile, buffer)){
                          
                          stringstream strin(buffer);
                          string id;
                          strin >> id;
                          string group;
                          if (Id_to_group.count(id) != 0){
                                                    group = Id_to_group[id];
                                                    Data[group].push_back(buffer);
                                                    count ++;
                                                    }                          
                         
                          } 
                          
    infile.close();
    infile.clear();
    return count;
    }

int Output_Data(const char * outfilename){
    
    for (map<string, vector <string> > :: iterator miter = Data.begin(); miter != Data.end(); miter ++){
        
        string outname = outfilename;
        outname += ".";
        outname += miter->first;
        
        ofstream outfile(outname.c_str(), ofstream::out);
        if (!outfile){
                      cerr << "Error: Cannot open file : " << outname << endl;
                      return 0;
                      }
        
        outfile << Title << endl; 
        for (int i = 0; i < (miter->second).size(); i ++)
            outfile << (miter->second)[i] << endl;
        
        outfile.close();
        outfile.clear();
        }
    return Data.size();
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    Load_Meta(Metafilename.c_str(), Key);
    cout << Load_Data(Infilename.c_str()) << " Samples loaded" << endl;
    cout << Output_Data(Outfilename.c_str()) << " Groups output" << endl;
    return 0;
    }
