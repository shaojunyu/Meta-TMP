// Updated at May 26, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include <stdlib.h>
#include <string.h>
#include <sys/dir.h>
#include <sys/stat.h>

#define BUFFER_SIZE 5000

#include "hash.h"
using namespace std;

string Check_Env(){
    
    if (getenv("ParallelMETA") == NULL){
                               
                               cerr << "Error: Please set the environment variable \"ParallelMETA\" to the directory" << endl;
                               exit(0);
                               
                                   }
    
    string path =  getenv("ParallelMETA");
    return path;
    }

int Check_Path(const char * path, int type){
    
    if (strlen(path) < 1) return 0;
    
    DIR *pDir = opendir(path);
            
    if(pDir!=NULL){
                  closedir(pDir);                   
                  if (type == 0){
                     string command = "rm -rf ";
                     command += path;
                     system(command.c_str());
                     mkdir(path, 0755);                          
                     }
                  }                  
    else 
         mkdir(path, 0755);           
    return 0;
    
    }

bool Check_Path(const char * path){
    
    if (strlen(path) < 1) return false;
    
    DIR *pDir = opendir(path);
    
    if (pDir != NULL){
             closedir(pDir);
             return true;
             }
    
    return false;
    }

bool Check_File(const char * file){
    
    fstream infile(file, ifstream::in);
    
    if (!infile){
                 
                 cerr << "Error: Cannot open file : " << file << endl;
                 return false;
                 }
    
    infile.close();
    infile.clear();
    
    return true;
    
    }

unsigned int Get_Count(const char * infilename){
         
         ifstream infile(infilename, ifstream::in);
         
         if (!infile){
                      
                      cerr << "Error: Cannot open file : " << infilename << endl;
                      return 0;
                      
                      }
         
         string buffer;
         unsigned int count = 0;
         
         while (getline(infile, buffer)){
               
               if (buffer[0] == '>') count ++;
               
               }
         
         infile.close();
         infile.clear();
         
         return count;
         }

int Check_Format(const char * infilename){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " <<infilename << endl;
                 return 0;
                 }
    
    string label;
    string seq;
    
    getline(infile,label);
    getline(infile, seq);
    
    infile.close();
    infile.clear();
    
    
    
    if ((label[0] == '>') && ((seq[0] == 'a') || (seq[0] == 't') || (seq[0] == 'c') || (seq[0] == 'g') || (seq[0] == 'A') || (seq[0] == 'T') || (seq[0] == 'C') || (seq[0] == 'G') || (seq[0] == 'N') || (seq[0] == 'n')))
    
       return 0;
    
    else if ((label[0] == '@') && ((seq[0] == 'a') || (seq[0] == 't') || (seq[0] == 'c') || (seq[0] == 'g') || (seq[0] == 'A') || (seq[0] == 'T') || (seq[0] == 'C') || (seq[0] == 'G')|| (seq[0] == 'N') || (seq[0] == 'n')))
         
         return 1;
         
    else {
         
         cerr << "File format error : " << infilename << " is not in FASTA or FASTQ format" << endl;
         return -1;
         
         }          
    }    

int Load_Copy_Number(string database_path, hash_map <string, float, std_string_hash> & cp_number){
    
    string copy_number_file = database_path + "16s_copy_number.txt";

    ifstream in_cp_number(copy_number_file.c_str(), ifstream::in);
    if (!in_cp_number){
                       cerr << "Error: Cannot open copy number table : " << copy_number_file << endl;
                       return 0;
                       } 
    
    int count = 0;
    string buffer;
    getline(in_cp_number, buffer);
    while(getline(in_cp_number, buffer)){
                                
                                stringstream strin(buffer);
                                string id;
                                float cp_no;
                                strin >> id >> cp_no;
                                
                                if (cp_number.count(id) == 0)
                                                           cp_number[id] = cp_no;
                                else cerr << "Warning, dup id : " << id << endl;
                                
                                count ++;
                                }
    
    in_cp_number.close();
    in_cp_number.clear();
    
    return count;
    }


int Load_ID(const char * idfilename, vector <string> & ID){
    
    ifstream infile(idfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open ID file: " << idfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          stringstream strin(buffer);
                          string id;
                          strin >> id;
                          
                          ID.push_back(id);
                          
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list){
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp);
                                                    
                          list.push_back(temp);
                          
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list, string prefix){ //with prefix
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp);
                                                    
                          list.push_back(prefix + temp);
                          
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids){//with id
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          vector <string> str_temp;
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp)
                                      str_temp.push_back(temp);
                          
                          if (str_temp.size() >= 2){
                                              
                                              ids.push_back(str_temp[0]);
                                              list.push_back(str_temp[1]);
                                              
                                              }
                          
                          else{                               
                               list.push_back(buffer);   
                               int name_begin = buffer.find_first_of('/');
                               int name_end = buffer.find_last_of('/');
                               if (name_begin == name_end)
                                              name_begin = 0;
                               else name_begin = buffer.rfind('/', name_end-1)+1;
                               ids.push_back(buffer.substr(name_begin, name_end - name_begin));
                               }
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids, string prefix){//with id
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          vector <string> str_temp;
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp)
                                      str_temp.push_back(temp);
                          
                          if (str_temp.size() >= 2){
                                              
                                              ids.push_back(str_temp[0]);
                                              list.push_back(prefix + str_temp[1]);
                                              
                                              }
                          
                          else{                               
                               list.push_back(prefix + buffer);   
                               int name_begin = buffer.find_first_of('/');
                               int name_end = buffer.find_last_of('/');
                               if (name_begin == name_end)
                                              name_begin = 0;
                               else name_begin = buffer.rfind('/', name_end-1)+1;
                               ids.push_back(buffer.substr(name_begin, name_end - name_begin));
                               }
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

/*
int Load_List_With_Id(const char * listfilename, vector <string> & list, vector <string> & ids){//with id
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          vector <string> str_temp;
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp)
                                      str_temp.push_back(temp);
                          
                          if (str_temp.size() >= 2){
                                              
                                              ids.push_back(str_temp[0]);
                                              list.push_back(str_temp[1]);
                                              
                                              }
                          
                          else{                               
                               
                               cerr << "Error: Please input the ID at column 1" << endl;
                               return 0;
                               }
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }
*/

void Make_list(const char * listname, const char * outpathname, vector <string> ids, int mode){ //0: classification.txt; 1: functions.txt; 2: id
    
    ofstream outfile(listname, ofstream::out);
    if (!outfile){
                  cerr << "Error: Cannot open output file : " << listname << endl;
                  return;
                  }
    
    for (int i = 1; i < ids.size(); i ++) // id with label
        switch(mode){
                     case 0: outfile << ids[i] << "\t" << outpathname << "/" << ids[i] << "/classification.txt" << endl; break;
                     case 1: outfile << ids[i] << "\t" << outpathname << "/" << ids[i] << "/functions.txt" << endl; break;
                     case 2: outfile << ids[i] << endl; break;
                     case 3: outfile << ids[i] << "\t" << outpathname << "/" << ids[i] << "/classification.rare.txt" << endl; break;
                     default: break;
                     }
    
    outfile.close();
    outfile.clear();
    }

void Add_list_prefix(const char * inlistname, const char * prefix, const char * outlistname){ //for prefix
     
     ifstream infile(inlistname, ios::in);
     if(!infile){
                 cerr << "Error: Cannot open input file: " << inlistname << endl;
                 return;
                 }
     
     ofstream outfile(outlistname, ios::out);
     if (!outfile){
                   cerr << "Error: Cannot open output file: " << outlistname << endl;
                   return;
                   }   
     
     string buffer;
     while(getline(infile, buffer)){
                           
                           stringstream strin(buffer);
                           string id;
                           string path;
                           strin >> id >> path;
                           path = prefix + path;
                           outfile << id << "\t" << path << endl;
                           }
     infile.close();
     infile.clear();
     
     outfile.close();
     outfile.clear();
     }
#endif
