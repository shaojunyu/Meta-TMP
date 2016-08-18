// Updated at Apr 5, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "hash.h"
#include "utility.h"
#include "version.h"

#include "table_format.h"

#ifndef COMP_H
#define COMP_H

/* gg_91
#define LeafN 22090
#define OrderN 22089
*/
 //gg_97
/*
#define LeafN 99322
#define OrderN 99321
*/
#define REG_SIZE 70

using namespace std;

class _Comp_Tree{
      
      public:
             _Comp_Tree(){
                           Database = "gg_13";
                           Database_path = Check_Env() + "/databases/" + Database + "/";
                           
                           Id_file = Database_path + "/tree/id.txt";
                           Order_file = Database_path + "/tree/order.txt";
                           
                           LeafN = 0;
                           OrderN = 0;
                           
                           /*
                           Dist_1 = new float[OrderN];
                           Dist_2 = new float[OrderN];
     
                           Order_1 = new int[OrderN];
                           Order_2 = new int[OrderN];
                           Order_d = new int[OrderN];
    
                           Id = new string[LeafN];
                           */
                           Init();
                          }
                          
            
             int Load_abd(const char * infilename, float * Abd, bool is_cp_correct);
             int Load_abd(const char * infilename, float * Abd);
             int Load_abd(_Table_Format table, float * Abd, int sample); //Load by table_format
             float Calc_sim(float * Abd_1, float * Abd_2);
             
             int Get_LeafN(){
                 return LeafN;
                 }
             
      private:
              string Database;
              string Database_path;
              string Id_file;
              string Order_file;

              /*
              float * Dist_1;
              float * Dist_2;
     
              int * Order_1;
              int * Order_2;
              int * Order_d;
    
              string * Id;
              */
              
              vector <float> Dist_1;
              vector <float> Dist_2;
              
              vector <int> Order_1;
              vector <int> Order_2;
              vector <int> Order_d;
              
              vector <string> Id; 
              
              hash_map <string, float, std_string_hash> Cp_number;
              
              void Init();
              int Load_id();
              int Load_order();   
              
              int LeafN;
              int OrderN;                         
              };

void _Comp_Tree::Init(){    
                    
     //load tree        
     LeafN = Load_id();
     OrderN = Load_order();
    
     //load cp number         
     Load_Copy_Number(Database_path, Cp_number);                    
     }

int _Comp_Tree::Load_id(){
     
     ifstream infile(Id_file.c_str(), ifstream::in);
     if (!infile){
                 cerr << "Error: Cannot open file : " << Id_file << endl;
                 return 0;
                 }
     
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){
                           if (buffer.size() == 0) continue;
                           Id.push_back(buffer);
                           count ++;
                           }
     
     infile.close();
     infile.clear();
     return count;
     }


int _Comp_Tree::Load_order(){
    ifstream infile(Order_file.c_str(), ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << Order_file << endl;
                 return 0;
                 }
    
     string buffer;
     int count = 0;
     while(getline(infile, buffer)){                           
                           if(buffer.size() == 0) continue;
                           stringstream strin(buffer);
                           int order_1 = 0;
                           int order_2 = 0;
                           int order_d = 0;
                           float dist_1 = 0;
                           float dist_2 = 0;
                           strin >> order_1 >> dist_1 >> order_2 >> dist_2 >> order_d;
                           Order_1.push_back(order_1);
                           Order_2.push_back(order_2);
                           Order_d.push_back(order_d);
                           Dist_1.push_back(dist_1);
                           Dist_2.push_back(dist_2);
                           count ++;                           
                           }
       
    infile.close();
    infile.clear(); 
    
    return count;
    }

int _Comp_Tree::Load_abd(const char * infilename, float * Abd, bool is_cp_correct){
    
    hash_map<string, float, std_string_hash> hash;
     
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr << "Error: Cannot open file : " << infilename << endl;
                 return 0;
                 }
    
    string buffer;
    int line_count = 0;    
    float total = 0;
    
    getline(infile, buffer); //title
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          line_count ++;
                                                    
                          string current_id;                                                                 
                          stringstream strin(buffer);
                          
                          strin >> current_id >> current_id;                          
                          if (hash.count(current_id) == 0)
                                                     hash[current_id] = 1;
                          else hash[current_id] ++;
                          }
        
    //cp_number_correct
    for (hash_map<string, float, std_string_hash>::iterator miter = hash.begin(); miter != hash.end(); miter ++){
        
        float cp_no = 1.0;
        
        if (is_cp_correct)
           if (Cp_number.count(miter->first) != 0)
                                          cp_no = Cp_number[miter->first];
        miter->second /= cp_no;
        total += miter->second;
        }
        
    //norm
    total /= 100.0;
    for (hash_map<string, float, std_string_hash>::iterator miter = hash.begin(); miter != hash.end(); miter ++)        
        miter->second /= total;
    
    for (int i = 0; i < LeafN; i++)
        if (hash.count(Id[i]) == 0)
           Abd[i] = 0;
        else Abd[i] = hash[Id[i]];        
    
    infile.close();
    infile.clear();
        
    return hash.size();
    }

int _Comp_Tree::Load_abd(const char * infilename, float * Abd){
    
    return Load_abd(infilename, Abd, true);
    }

int _Comp_Tree::Load_abd(_Table_Format table, float *Abd, int sample){ 
    
    int count = 0;
    
    for (int i = 0; i < LeafN; i ++){
        Abd[i] = table.Get_Abd_By_Feature(sample, "OTU_" + Id[i]) * 100.0;
        if (Abd[i] > 0) count ++;
        }
    return count;
    }

float _Comp_Tree::Calc_sim(float * Abd_1, float * Abd_2){
      
      float Reg_1[REG_SIZE];
      float Reg_2[REG_SIZE];
      
      float total = 0;
      
      for(int i = 0; i < OrderN; i++){
              
              int order_1 = Order_1[i];
              int order_2 = Order_2[i];
              int order_d = Order_d[i] + REG_SIZE;
              
              float dist_1 = 1- Dist_1[i];
              float dist_2 = 1- Dist_2[i];
              
              float c1_1;
              float c1_2;
              
              float c2_1;
              float c2_2;
                                
              if (order_1 >= 0){
                          
                          c1_1 = Abd_1[order_1];
                          c1_2 = Abd_2[order_1];
                          
                          }
              else {
                   c1_1 = Reg_1[order_1 + REG_SIZE];
                   c1_2 = Reg_2[order_1 + REG_SIZE];
                   }
              
              if (order_2 >= 0){
                          
                          c2_1 = Abd_1[order_2];
                          c2_2 = Abd_2[order_2];
                          
                          }
              else {
                   c2_1 = Reg_1[order_2 + REG_SIZE];
                   c2_2 = Reg_2[order_2 + REG_SIZE];
                   }
              //min
              float min_1 = (c1_1 < c1_2)?c1_1:c1_2;
              float min_2 = (c2_1 < c2_2)?c2_1:c2_2;
              
              total += min_1;
              total += min_2;
              
              //reduce
              Reg_1[order_d] = (c1_1 - min_1) * dist_1 + (c2_1 - min_2) * dist_2;
              Reg_2[order_d] = (c1_2 - min_1) * dist_1 + (c2_2 - min_2) * dist_2;
              
              }
      
      return total/100.0; //scale 0-1
      
      }

#endif
