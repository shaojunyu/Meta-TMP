// Updated at May 19, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>

#include "hash.h"

#ifndef TAXONOMY_H
#define TAXONOMY_H

#define LEVEL 7 //(0: Kingdom, 5: Genus, 7: OTU)

using namespace std;

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
    
unsigned int Parse_Taxonomy(const char * infilename, const char * outfilename, hash_map<string, string, std_string_hash> * table, hash_map<string, unsigned int, std_string_hash> * taxa_table, int Q, int * a_diver, bool is_paired){
    
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
    string last_seq_id = "";
    
    set <string> a_diver_table [LEVEL];
        
    while(getline(infile, buffer)){
                          if (buffer.size()==0) continue;
                          stringstream strin(buffer);
                          string seq_id, database_id, flag;
                          int mapq;
                          strin >> seq_id >> flag >> database_id >> mapq >> mapq;
                          
                          total_count ++;
                          if (is_paired)
                             getline(infile, buffer);
                                                       
                          if (database_id[0] == '*') continue;
                          if (mapq < Q) continue;
                          
                          if (seq_id == last_seq_id) continue; //Duplications
                          else last_seq_id = seq_id;
                          
                          //Generate out the abundance count
                          if (taxa_table->count(database_id) > 0) 
                                                             (*taxa_table)[database_id] ++;
                                                             
                          else
                                             
                              (*taxa_table)[database_id] = 1;
                          
                          outfile << seq_id << "\t" << database_id << "\t" << flag << "\t" << mapq << "\t";
                              
                          //Output the taxonomy annotation
                          if (table->count(database_id) != 0){
                          
                             string taxa_buffer = (*table)[database_id];
                             outfile << taxa_buffer << endl;
                             
                             //Alpha diversity
                             int start = 0;
                             int end = 0;
                               
                             for(int i = 0; i < LEVEL - 1; i ++){
                                       end = taxa_buffer.find(';', start + 1);
                                       string taxa = taxa_buffer.substr(start, end - start);
                                       if (taxa.find("otu_") != string::npos) break;
                                       a_diver_table[i].insert(taxa);
                                       start = end + 2;
                                       if (start >= taxa_buffer.size()) break;
                                       }
                                                          
                             }
                          else 
                             outfile << "Unclassified" << endl;                                                            
                          
                          a_diver_table [LEVEL - 1].insert(database_id); //OTU
                          
                          match_count ++;                                                    
                          }
                          
    
    for (int i = 0; i < LEVEL; i ++)
        a_diver[i] = a_diver_table[i].size();
    
    infile.close();
    infile.clear();
    outfile.close();
    outfile.clear();
    
    cout << endl << total_count << " sequences are loaded " << endl;
    
    return match_count;
    
    }

unsigned int Out_Taxonomy(const char * infilename, string this_path, string database_path, string out_path, int Q, int * a_diver, int is_paired){
    
    string taxonomyfilename = database_path + "taxonomy_annotation.txt"; 
    string outfilename = out_path + "/classification.txt";

    hash_map<string, string, std_string_hash> table;
    hash_map<string, unsigned int, std_string_hash> taxa_table;
    
    unsigned int taxonomy_count = 0;
    unsigned int match_count = 0;
    
    taxonomy_count = Read_Taxonomy(taxonomyfilename.c_str(), table);
    
    cout << endl << taxonomy_count << " taxonomy annotations are loaded" << endl;
    
    match_count = Parse_Taxonomy(infilename, outfilename.c_str(), &table, &taxa_table, 0, a_diver, is_paired);
    
    cout << endl << match_count << " taxonomy annotations are parsed out" << endl << endl;
    
    return match_count;
    
    }

#endif
/*

int main(int argc, char * argv[]){
    
    string out_path = argv[3];
    string database = argv[2];
    
    Out_Taxonomy(argv[1], database, out_path, 0);

    
    return 0;
    
    }
*/
