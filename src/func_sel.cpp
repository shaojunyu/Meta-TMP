// Updated at May 13, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>

#include "class_func.h"

using namespace std;

// Parameters Def

string out_file = "functions_category";

string database = "gg_13";
string database_path;

string Listfilename;
string Listprefix;

vector <string> Infilename;
vector <string> Sam_name;

int Level = 2; //category level, default is 2

float Max_abd = 0;
float Min_abd = 0;

bool Is_print = false;

string Func_level[4] = {"Pathway level 1", "Pathway level 2", "Pathway level 3", "KO"};

int printhelp(){
    
    cout << "Function selection version " << Version << endl;
    cout << "Usage : func_sel [Options] Value" << endl;
    cout << "Options : " << endl;
    cout << "\t-l Input filename list [Required]" << endl;    
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    cout << "\t-o Output file" << endl;
    cout << "\t  To appoint the result output path, default is \"functions_category\"" << endl;
    cout << "\t-L KEGG Pathway level" << endl;
    cout << "\t  Level 1, 2, 3 or 4 (KO number), default is 2" << endl;
    cout << "\t-P Print distribution barchart, T(rue) or F(alse), default is F" << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){
    
    database_path =  Check_Env() + "/databases/" + database + "/";
    
    if (argc ==1) 
		printhelp();
    
    int i = 1;
    int sam_num = 0;

    Level = 2;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){                            
                            case 'l': Listfilename = argv[i+1]; break;  
                            case 'p': Listprefix = argv[i+1]; break;                                                                  
                            case 'o': out_file = argv[i+1]; break;                                                       
                            case 'L': Level = atoi(argv[i+1]); break; //by cate     
                            case 'P': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Is_print = true; break;                        
                            case 'h': printhelp(); break;

                            default : printf("Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }
    
    //check
    if ((Level > 4) || (Level < 1)){               
               cerr << "Warning: Level must between 1 and 4, change to default (2)" << endl;
               Level = 2;
               }
    //load
    sam_num = Load_List(Listfilename.c_str(), Infilename, Sam_name, Listprefix);      
                
    return sam_num;
    }

int main(int argc, char * argv[]){
    
    int sam_num = Parse_Para(argc, argv);
    
    _KO_OTU_Table_All KOs(database_path, sam_num, 1); // sim mode
    
    cout << endl << "Functional Selection Starts" << endl;
    cout << endl << "Level is " << Func_level[Level-1] << endl;
       
    cout << "Total Sample Number is " << sam_num << endl;
        
    for (int i = 0; i< sam_num; i++){

        cout << KOs.Load_Sample_By_Table(Infilename[i].c_str(), Sam_name[i], i) << " genes are loaded" << endl;

    }
            
    if (Level > 0)
       cout << endl << KOs.Output_By_Category(out_file.c_str(), Level - 1, Max_abd, Min_abd) << " Category have been parsed out" << endl;;        
        
    if (Is_print){
                  char command[BUFFER_SIZE];
                  sprintf(command, "Rscript %s/Rscript/PM_Distribution.R -i %s -o %s", Check_Env().c_str(), (out_file + ".Abd").c_str(), (out_file + ".Abd.distribution.pdf").c_str());
                  system(command);
                  }
                  
    cout << endl << "Functional Selection Finished"<< endl;

    return 0;
    }
