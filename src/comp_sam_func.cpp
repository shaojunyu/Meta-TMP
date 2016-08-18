// Updated at May 13, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// version 3.1 or above with _Table_Format

#include <iostream>

#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "utility.h"
#include "version.h"
#include "comp_sam_func.h"
#include "table_format.h"

//using namespace std;

string database;

string Listfilename;
string Listprefix;

string Queryfile1;
string Queryfile2;

string Tablefilename;
string Outfilename;

int Coren = 0;

int Dist_metric = 0; //0: cos 1: eu;

bool Is_sim; //true: sim, false: dist;

int Mode = 0; //0: single, 1: multi, 2: multi_table
bool Is_heatmap;
int Cluster = 2;

int printhelp(){
    
    cout << "Comp_sam_func Version : " << Version << endl;
    cout << "Usage : " << endl;
    cout << "comp-sam-func [-option] value" << endl;
    cout << "\toption : " << endl;

    cout << "\t-i Two samples path for single sample comparison [Conflict with -l and -T]" << endl;
    cout << "\t-l Sample name list table for multi-sample comparison [Conflict with -i and -T]" << endl;
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    cout << "\t-T Sample KO table (*.Abd) for multi-sample comparison [Conflict with -i and -l]" << endl;
    cout << "\t-o Result output file, default is to output on screen" << endl;
    cout << "\t-d Output format, similarity (F) or distance (T), default is F" << endl;
    cout << "\t-D Distance Metric, 0: Cosine; 1: Euclidean, default is 0" << endl;
    cout << "\t-P Print heatmap and clusters, T(rue) or F(alse), default is F" << endl;
    cout << "\t-c Cluter number, default is 2" << endl;
    cout << "\t-t Cpu core number, default is auto" << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    
    return 0;
    
    };
    
int Parse_Para(int argc, char * argv[]){
    
    database = "gg_13"; 
    
    Coren = 0;
    Mode = true; //default is single;
    
    Is_sim = true;
    Is_heatmap = false;
    
    int i = 1;
      
      if (argc ==1) 
		printhelp();
      
      while(i<argc){
         if (argv[i][0] != '-') {
                           cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'i': Queryfile1 = argv[i+1]; Queryfile2 = argv[i+2]; i++; Mode = 0; break;
                            case 'l': Listfilename = argv[i+1]; Mode = 1; break;
                            case 'p': Listprefix = argv[i+1]; break;
                            case 'T': Tablefilename = argv[i+1]; Mode = 2; break;
                            case 'o': Outfilename = argv[i+1]; break;
                                                        
                            case 'd': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_sim = false; break;
                            case 'D': Dist_metric = atoi(argv[i+1]); break;
                            case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_heatmap = true; break;
                            case 'c': Cluster = atoi(argv[i+1]); break;
                                 
                            case 't': Coren = atoi(argv[i+1]); break;         
                            case 'h': printhelp(); break;
                            default : cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break; 
                            }
         i+=2;
         }
    
    if ((Dist_metric > 1) || (Dist_metric < 0)){
                     cerr << "Warning: Distance metric (-D) must be 0 or 1, change to default (0)" << endl;
                     Dist_metric = 0;
                     }     
    int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((Coren <= 0) || (Coren > max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    Coren = max_core_number;
                    }  
    if (Cluster <= 0){
                cerr << "Warning: cluster number must be larger than 0, change to default (2)" << endl;
                }
    }

void Output_Matrix(const char * outfilename, int n, float *sim_matrix, bool is_sim, vector <string> sam_name){
     
     ofstream outfile(outfilename, ofstream::out);
     if (!outfile){
                   cerr << "Error: Cannot open output file : " << outfilename << endl;
                   return; 
                   }
     
     //Label
     for(int i = 0; i < n; i ++)
             outfile << "\t" << sam_name[i];
     outfile << endl;
     
     for(int i = 0; i < n; i ++){
             outfile << sam_name[i];
             
             for (int j = 0; j < n; j ++)                 
                 if (is_sim){                 
                    if (i == j) outfile << "\t" << 1.0;
                    else outfile << "\t" << 1.0- sim_matrix[i * n + j];
                    }
                 else {
                      if (i == j) outfile << "\t" << 0.0;
                      else outfile << "\t" << sim_matrix[i * n + j];
                      }
                 
             outfile << endl;
             }
                 
     outfile.close();
     outfile.clear();
     }

void Single_Comp(){
     
     _Comp_Tree_Func comp_tree_func;
    
    float * abd_1 = new float [comp_tree_func.Get_GeneN()];
    float * abd_2 = new float [comp_tree_func.Get_GeneN()];
    
    cout << comp_tree_func.Load_Gene_Count(Queryfile1.c_str(), abd_1) << " loaded" << endl;
    cout << comp_tree_func.Load_Gene_Count(Queryfile2.c_str(), abd_2) << " loaded" << endl;
    
    float dist = 0;
    if (Dist_metric == 0)
                    dist = comp_tree_func.Calc_Dist_Cos(abd_1, abd_2);
    else dist = comp_tree_func.Calc_Dist_E(abd_1, abd_2);
    
    if (Is_sim)
       cout << (1.0 - dist) << endl;
    else cout << dist << endl;
    
    }

void Multi_Comp(){
    
    _Comp_Tree_Func comp_tree_func;
         
     //load list
    vector <string> sam_name;
    vector <string> file_list;
    
    int file_count = 0;
    
    file_count = Load_List(Listfilename.c_str(), file_list, sam_name, Listprefix);
            
    //load abd
    float **Abd = new float * [file_count];
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [comp_tree_func.Get_GeneN()];
        cout << comp_tree_func.Load_Gene_Count(file_list[i].c_str(), Abd[i]) << " KOs in file " << i + 1 << endl;
        }
    
    cout << file_count << " files loaded" << endl;
    
    //make order
    int * order_m = new int [file_count * (file_count - 1) / 2];
    int * order_n = new int [file_count * (file_count - 1) / 2];
    int iter = 0;
    
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m[iter] = i;
            order_n[iter] = j;
            iter ++;
            }
    
    //openmp;
    
    float * sim_matrix = new float [file_count * file_count];
    memset(sim_matrix, 0, file_count * file_count * sizeof(float));
        
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < iter; i ++){
        
        int m = order_m[i];
        int n = order_n[i];
        
        if (Dist_metric == 0)
                    sim_matrix[m * file_count + n] = comp_tree_func.Calc_Dist_Cos(Abd[m], Abd[n]);
        else sim_matrix[m * file_count + n] = comp_tree_func.Calc_Dist_E(Abd[m], Abd[n]);
        //sim_matrix[m * file_count + n] = comp_tree_func.Calc_Dist_Cos(Abd[m], Abd[n]);
        sim_matrix[n * file_count + m] = sim_matrix[m * file_count + n];
        }
    
    Output_Matrix(Outfilename.c_str(), file_count, sim_matrix, Is_sim, sam_name);
    
    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];
     
     if (Is_heatmap){
                    char command[BUFFER_SIZE];
                    sprintf(command, "Rscript %s/Rscript/PM_Heatmap.R -d %s -o %s", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".heatmap.pdf").c_str());
                    system(command);
                    sprintf(command, "Rscript %s/Rscript/PM_Hcluster.R -d %s -o %s -c %d", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".clusters.pdf").c_str(), Cluster);
                    system(command);
                    }        
     };

void Multi_Comp_Table(){
          
     _Comp_Tree_Func comp_tree_func;
     _Table_Format abd_table(Tablefilename.c_str());
         
    int file_count = abd_table.Get_Sample_Size();
    
    //load abd
    float **Abd = new float * [file_count];
    for (int i = 0; i < file_count; i ++){
        Abd[i] = new float [comp_tree_func.Get_GeneN()];
        cout << comp_tree_func.Load_Gene_Count(abd_table, Abd[i], i) << " KOs in file " << i + 1 << endl;
        }
    
    cout << file_count << " files loaded" << endl;
    
    //make order
    int * order_m = new int [file_count * (file_count - 1) / 2];
    int * order_n = new int [file_count * (file_count - 1) / 2];
    int iter = 0;
    
    for (int i = 0; i < file_count - 1; i ++)
        for (int j = i + 1; j < file_count; j ++){            
            order_m[iter] = i;
            order_n[iter] = j;
            iter ++;
            }
    
    //openmp;
    
    float * sim_matrix = new float [file_count * file_count];
    memset(sim_matrix, 0, file_count * file_count * sizeof(float));
        
    omp_set_num_threads(Coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < iter; i ++){
        
        int m = order_m[i];
        int n = order_n[i];
        
        if (Dist_metric == 0)
                    sim_matrix[m * file_count + n] = comp_tree_func.Calc_Dist_Cos(Abd[m], Abd[n]);
        else sim_matrix[m * file_count + n] = comp_tree_func.Calc_Dist_E(Abd[m], Abd[n]);
        //sim_matrix[m * file_count + n] = comp_tree_func.Calc_Dist_Cos(Abd[m], Abd[n]);
        sim_matrix[n * file_count + m] = sim_matrix[m * file_count + n];
        }
    
    Output_Matrix(Outfilename.c_str(), file_count, sim_matrix, Is_sim, abd_table.Get_Sample_Names());
    
    for (int i = 0; i < file_count; i ++)
        delete [] Abd[i];
     
     if (Is_heatmap){
                    char command[BUFFER_SIZE];
                    sprintf(command, "Rscript %s/Rscript/PM_Heatmap.R -d %s -o %s", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".heatmap.pdf").c_str());
                    system(command);
                    sprintf(command, "Rscript %s/Rscript/PM_Hcluster.R -d %s -o %s -c %d", Check_Env().c_str(), Outfilename.c_str(), (Outfilename + ".clusters.pdf").c_str(), Cluster);
                    system(command);
                    }        
     }

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);         
                  
    switch (Mode) {
           case 0: Single_Comp(); break;
           case 1: Multi_Comp(); break;
           case 2: Multi_Comp_Table(); break;
           default: break;
           }
    
    return 0;
    }
