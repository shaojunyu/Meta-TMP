// Updated at Mar 13, 2016
// Updated by Xiaoquan Su
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

using namespace std;

class _Barcode{
      
      public:
             _Barcode(){
                        }
             _Barcode(const char * infilename, const char * outpath, int format){                    
                     Load(infilename);                     
                     Output_Init(outpath, format);                     
                     Seq_count = vector <unsigned int> (Sample.size(), 0);             
                     }                   
             
            bool Get_Match(vector <string> seq){//trim seq                   
                      
                      if (seq.size() < 2) return false;
                      
                      for (set <int>::iterator siter = Barcode_length.begin(); siter != Barcode_length.end(); siter++){
                          
                          if (seq.size() < *siter) continue;
                          string barcode = seq[1].substr(0, *siter);                           
                          if (Hash.count(barcode) != 0) {
                                                  int index = Hash[barcode];
                                                  string outfilename = Outfile[index];
                                                  ofstream outfile(outfilename.c_str(), ios::app);
                                                  if (!outfile){
                                                                cerr << "Error: Cannot open output file: " << outfilename << endl;
                                                                return false;
                                                                }
                                                  seq[1] = seq[1].substr(*siter, seq[1].size() - *siter);
                                                  for (int i = 0; i < seq.size(); i ++)
                                                      outfile << seq[i] << endl;
                                                      
                                                  outfile.close();
                                                  outfile.clear();
                                                  
                                                  Seq_count[index] ++;
                                                  return true;
                                                  }                                     
                                     }
                      return NULL;
                      }                                   
           
           int Get_Size(){
               return Barcode.size();
               }
                                
          void Output_Map(const char * outfilename){
               
               ofstream outfile(outfilename, ios::app);
               if (!outfile){
                             cerr << "Error: Cannot open map file: " << outfilename << endl;
                             return;
                             }
               
               outfile << "SampleID\tBarcode\tSequenceCount" << endl;
               
               for (int i = 0; i < Outfile.size(); i ++)
                   outfile << Sample[i] << "\t" << Barcode[i] << "\t" << Seq_count[i] << endl;
               
               outfile.close();
               outfile.clear();
               }
                      
      private:              
              vector <string> Barcode;
              set <int> Barcode_length;
              vector <string> Sample;
              vector <string> Outfile;
              vector <unsigned int> Seq_count;
              map<string, int> Hash;
                                          
              void Load(const char * infilename){
                   ifstream infile(infilename, ifstream::in);
                   if (!infile){
                                cerr << "Error: Cannot open input barcode file : " << infilename << endl;
                                exit(0);
                                return;
                                }
                   string buffer;
                   unsigned int count = 0;
                   while(getline(infile, buffer)){
                                         stringstream strin(buffer);
                                         string barcode, sample;
                                         strin >> barcode >> sample;
                                         Barcode.push_back(barcode);
                                         Barcode_length.insert(barcode.size());
                                         Sample.push_back(sample);
                                         Hash[barcode] = count;
                                         count ++;
                                         }                   
                   infile.close();
                   infile.clear();
                   }
              
              void Output_Init(const char * outpath, int format){
                   for (int i = 0; i < Sample.size(); i ++){                       
                       string outfilename = outpath;
                       outfilename += "/";
                       outfilename += Sample[i];   
                       if (format == 0)
                          outfilename += ".fa";
                       else 
                            outfilename += ".fq";                
                       
                       Outfile.push_back(outfilename);
                                                                                        
                       ofstream  outfile_temp(outfilename.c_str(), ofstream::out);
                       if (!outfile_temp){
                             cerr << "Error: Cannot open file : " << outfilename << endl;
                             return;
                             }
                       outfile_temp.close();
                       outfile_temp.clear();
                       }
                   }                                        
      };
      
class _Group{
      
      public:
             _Group(){}
             
             _Group(const char * infilename, const char * outpath, int format){                     
                                          
                     Load(infilename);                     
                     Output_Init(outpath, format);                                                       
                     }
             /*
             ofstream * Get_Match(string id){//trim seq                                                                            
                        if (Hash_sample.count(id) != 0){ 
                           Hash_count[Hash_sample[id]] ++;                                                         
                           return Hash_outfile[Hash_sample[id]];
                           }                                                             
                        else return NULL;                                     
                      }                  
           */
           
           bool Get_Match(vector <string> seq){
                if (seq.size() < 2) return false;
                string id = seq[0].substr(1, seq[0].size() - 1);
                if (Hash_sample.count(id) == 0) return false;
                string outfilename = Hash_outfile[Hash_sample[id]];
                ofstream outfile(outfilename.c_str(), ios::app);
                if (!outfile){
                              cerr << "Error: Cannot open output file: " << outfilename << endl;
                              return false;
                              } 
                for (int i = 0; i < seq.size(); i ++)
                    outfile << seq[i] << endl;
                Hash_count[Hash_sample[id]] ++;
                return true;
                }
           
           int Get_Size(){
               return Hash_outfile.size();
               }
                                         
           void Output_Map(const char * outfilename){
               
               ofstream outfile(outfilename, ios::app);
               if (!outfile){
                             cerr << "Error: Cannot open map file: " << outfilename << endl;
                             return;
                             }
               outfile << "SampleID\tSequenceCount" << endl;
               
               for (hash_map<string, unsigned int, std_string_hash> ::iterator miter = Hash_count.begin(); miter != Hash_count.end(); miter ++)
                   outfile << miter->first << "\t" << miter->second << endl;
               
               outfile.close();
               outfile.clear();
               }                     
      private:              
              
              hash_map<string, string, std_string_hash> Hash_sample; // seq id to sample
              hash_map <string, string, std_string_hash> Hash_outfile; // sample to file
              hash_map <string, unsigned int, std_string_hash> Hash_count; // sample to count
              
              void Load(const char * infilename){
                   ifstream infile(infilename, ifstream::in);
                   if (!infile){
                                cerr << "Error: Cannot open input group file : " << infilename << endl;
                                exit(0);
                                return;
                                }
                   string buffer;
                   unsigned int count = 0;
                   while(getline(infile, buffer)){
                                         stringstream strin(buffer);
                                         string id, sample;
                                         strin >> id >> sample;
                                         
                                         if (Hash_sample.count(id) != 0){
                                                                   cerr << "Warning: Duplicate sequence id : " << id << endl;
                                                                   continue;
                                                                   }
                                         Hash_sample[id] = sample;
                                         Hash_outfile[sample] = sample;
                                         Hash_count[sample] = 0;
                                         count ++;
                                         }                   
                   infile.close();
                   infile.clear();
                   }
              
              void Output_Init(const char * outpath, int format){
                   for (hash_map<string, string, std_string_hash> ::iterator miter = Hash_outfile.begin(); miter != Hash_outfile.end(); miter ++){                       
                       string outfilename = outpath;
                       outfilename += "/";
                       outfilename += miter->first;   
                       if (format == 0)
                          outfilename += ".fa";
                       else 
                            outfilename += ".fq";                
                                                                                        
                       miter->second = outfilename;
                       ofstream outfile_temp(outfilename.c_str(), ofstream::out);
                       if (!outfile_temp){
                             cerr << "Error: Cannot open file : " << outfilename << endl;
                             return;
                             }
                       outfile_temp.close();
                       outfile_temp.clear();
                       }
                   }
      };

string Infilename;
string Barcodefilename;
string Groupfilename;
string Outpath;
int Format = 0;
int Split_Mode = -1;  //0:barcode; 1: group

int printhelp(){
    
    cout << "Split-Seq: " << Version << endl;
    cout << "Usage : " << endl;
    cout << "split-seq [-options] value" << endl;
    cout << "\toption : " << endl;

    cout << "\t-i Input sequence file in FASTA or FASTQ format [Required]" << endl;
    cout << "\t-b Input barcode file [Conflict with -g]" << endl;
    cout << "\t-g Input group file [Conflict with -b]" << endl;
    cout << "\t-o Result output path, default is \"Out\"" << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    Outpath = "Out";
    Format = 0;
    Split_Mode = -1;
    
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
                            case 'b': Barcodefilename = argv[i+1]; 
                                      if (Split_Mode == 1){
                                         cerr << "Error: -b conflicts with -g" << endl;
                                         exit(0);
                                         }
                                      Split_Mode = 0;
                                      break;
                            case 'g': Groupfilename = argv[i+1]; 
                                      if (Split_Mode == 0){
                                         cerr << "Error: -g conflicts with -b" << endl;
                                         exit(0);
                                         }
                                      Split_Mode = 1;
                                      break;
                            case 'o': Outpath = argv[i+1]; break;                                                  
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }        
    
    Check_Path(Outpath.c_str(), 1);
    Format = Check_Format(Infilename.c_str());
    if (Format < 0) exit(0);
    }

void Split_By_Barcode(){
     
     _Barcode barcode(Barcodefilename.c_str(), Outpath.c_str(), Format);
    
    int sample_count = barcode.Get_Size();
    cout << sample_count << " Samples loaded" << endl;
    
    ifstream infile(Infilename.c_str(), ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open infile : " << Infilename << endl;
                 return;
                 }
    string buffer;
    unsigned int count = 0;
    unsigned int total_count = 0;
    unsigned int map_count = 0;
    while(getline(infile, buffer)){
                          
                          if ((buffer[0] != '>') && (buffer[0] != '@')){
                                         cerr << "Warning: Sequence format warning at line : " << count + 1 << endl;
                                         continue;
                                         } 
                          count ++;
                          vector <string> seq;
                          seq.push_back(buffer);                                          
                          getline(infile, buffer);
                          seq.push_back(buffer);
                          count ++;
                          
                          if (Format == 1){ //fastq
                                     getline(infile, buffer);
                                     seq.push_back(buffer);
                                     count ++;
                                     getline(infile, buffer);
                                     seq.push_back(buffer);
                                     count ++;
                                     }
                          
                          if (barcode.Get_Match(seq))
                             map_count ++;
                                                    
                          total_count ++;
                          }
    
    infile.close();
    infile.clear();
    
    cout << total_count << " Sequences are loaded" << endl;
    cout << map_count << " Sequences are output" << endl;
    
    //print report
    ofstream outfile((Outpath + "/Analysis_Report.txt").c_str(), ofstream::out);
    if (!outfile){
                  cerr << "Error: Cannot open file : " << Outpath + "/Analysis_Report.txt" << endl;
                  return;
                  }
    
    outfile << "Input sequence: " << Infilename << endl;
    outfile << "Barcode file: " << Barcodefilename << endl;
    outfile << "Number of samples: " << sample_count << endl;
    outfile << "Number of input sequences: " << total_count << endl;
    outfile << "Number of output sequences: " << map_count << endl;
    
    outfile.close();
    outfile.clear(); 
    
    barcode.Output_Map((Outpath + "/Analysis_Report.txt").c_str());        
    }

void Split_By_Group(){
     
     _Group group(Groupfilename.c_str(), Outpath.c_str(), Format);
    
    int sample_count = group.Get_Size();
    cout << sample_count << " Samples loaded" << endl;
    
    ifstream infile(Infilename.c_str(), ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open infile : " << Infilename << endl;
                 return;
                 }
    string buffer;
    unsigned int count = 0;
    unsigned int total_count = 0;
    unsigned int map_count = 0;
    while(getline(infile, buffer)){
                          
                          if ((buffer[0] != '>') && (buffer[0] != '@')){
                                         cerr << "Warning: Sequence format warning at line : " << count + 1 << endl;
                                         continue;
                                         } 
                          count ++;
                          vector <string> seq;
                          seq.push_back(buffer);                                          
                          getline(infile, buffer);
                          seq.push_back(buffer);
                          count ++;
                          
                          if (Format == 1){ //fastq
                                     getline(infile, buffer);
                                     seq.push_back(buffer);
                                     count ++;
                                     getline(infile, buffer);
                                     seq.push_back(buffer);
                                     count ++;
                                     }
                          
                          if (group.Get_Match(seq))
                             map_count ++;
                                                    
                          total_count ++;
                          }
    
    infile.close();
    infile.clear();
    
    cout << total_count << " Sequences are loaded" << endl;
    cout << map_count << " Sequences are output" << endl;
    
    //print report
    ofstream outfile((Outpath + "/Analysis_Report.txt").c_str(), ofstream::out);
    if (!outfile){
                  cerr << "Error: Cannot open file : " << Outpath + "/Analysis_Report.txt" << endl;
                  return;
                  }
    
    outfile << "Input sequence: " << Infilename << endl;
    outfile << "Group file: " << Groupfilename << endl;
    outfile << "Number of samples: " << sample_count << endl;
    outfile << "Number of input sequences: " << total_count << endl;
    outfile << "Number of output sequences: " << map_count << endl;
    
    outfile.close();
    outfile.clear(); 
    
    group.Output_Map((Outpath + "/Analysis_Report.txt").c_str());
    }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
                        
    if (Split_Mode == 0) Split_By_Barcode();
    else if (Split_Mode == 1) Split_By_Group();
    
    return 0;
    }
