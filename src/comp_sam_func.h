// Updated at Apr 5, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>

#include "hash.h"
#include "utility.h"

#include "table_format.h"

#ifndef COMP_FUNC_H
#define COMP_FUNC_H

//using namespace std;

//#define GeneN 6909

class _Comp_Tree_Func{
      
      public:
             _Comp_Tree_Func(){                               
                               Database = "gg_13";
                               Database_path = Check_Env() + "/databases/" + Database + "/";                           
                               Id_file = Database_path + "KO/ko_id.tab";
                               GeneN = Load_Id();
                               }
             
             int Load_Gene_Count(const char * infilename, float * abd);
             int Load_Gene_Count(_Table_Format table, float * abd, int sample);
             float Calc_Dist_Cos(float * abd_m, float * abd_n);
             float Calc_Dist_E(float * abd_m, float * abd_n); //Weighted Euclidean dist
             
             int Get_GeneN(){
                 return GeneN;
                 }
             
      private:
              //void Init()
              int Load_Id();
              
              vector <string> Gene;
              
              string Database;
              string Database_path;
              string Id_file;
              
              int GeneN;              
      };

int _Comp_Tree_Func::Load_Gene_Count(const char * infilename, float * abd){
    
    //tr1::unordered_map<string, float, StrHash, StrCompare> table;
    hash_map<string, float, std_string_hash> table;
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open input file : " << infilename << endl;
                 exit(0);
                 }
    
    int count = 0;
    string buffer;
    getline(infile, buffer);//title
    while(getline(infile, buffer)){                          
                          stringstream strin(buffer);  
                          string gene;
                          string temp;
                          float gene_count;
                          
                          //strin >> gene >> temp >> temp >> gene_count;
                          strin >> gene >> gene_count;
                          /*
                          if (table.find(gene) == table.end())
                                                table[gene] = gene_count;
                          */
                          if (table.count(gene) == 0)
                                                table[gene] = gene_count;
                          else table[gene] += gene_count; 
                          count ++;
                          }
    
    for (int i = 0; i < GeneN; i ++){
        abd[i] = 0;
        if (table.count(Gene[i]) != 0)
           abd[i] = table[Gene[i]];
        }
                 
    infile.close();
    infile.clear();
    
    return count;
    }

int _Comp_Tree_Func::Load_Gene_Count(_Table_Format table, float * abd, int sample){
    
    int count = 0;
    
    for (int i = 0; i < GeneN; i ++){
        abd[i] = table.Get_Abd_By_Feature(sample, Gene[i]);
        if (abd[i] > 0) count ++;
        }
    return count;
    }

float _Comp_Tree_Func::Calc_Dist_Cos(float * abd_m, float * abd_n){
      
      float f = 0;
      
      float f_m_sum = 0;
      float f_n_sum = 0;
            
      for (int i = 0; i < GeneN; i ++){
          
          float fm = abd_m[i];
          float fn = abd_n[i];
          
          f += fm * fn;
          
          f_m_sum += fm * fm;
          f_n_sum += fn * fn;
                   
          }
      
      f = f / (sqrt(f_m_sum * f_n_sum));
      
      return (1 - f);
      
      }

float _Comp_Tree_Func::Calc_Dist_E(float * abd_m, float * abd_n){
      
      float abd_m_norm[GeneN];
      float abd_n_norm[GeneN];
      
      float sum_m = 0;
      float sum_n = 0;
      
      //Norm
      for (int i = 0; i < GeneN; i ++){
          
          abd_m_norm[i] = abd_m[i];
          abd_n_norm[i] = abd_n[i];
          
          sum_m += abd_m_norm[i];
          sum_n += abd_n_norm[i];          
          }
      
      for (int i = 0; i < GeneN; i ++){
          abd_m_norm[i] /= sum_m;
          abd_n_norm[i] /= sum_n;
          }
      
      //calc
      float f = 0;
      
      for (int i = 0; i < GeneN; i ++)
          f += pow(abd_m_norm[i] - abd_n_norm[i], 2);
      
      return sqrt(f);
      }

int _Comp_Tree_Func::Load_Id(){
    
    ifstream in_idfile(Id_file.c_str(), ifstream::in);
    if (!in_idfile){
                    cerr << "Error: Open KO ID file error : " << Id_file << endl;
                    return 0;
                    }
    string buffer;
    int count = 0;
    while(getline(in_idfile, buffer)){
                          if (buffer.size() == 0) continue;
                          Gene.push_back(buffer);
                          count ++;
                          }
    in_idfile.close();
    in_idfile.clear();
    
    return count;
    }

#endif
