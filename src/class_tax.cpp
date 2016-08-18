// Updated at May 13, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <stdlib.h>

#include <sys/dir.h>
#include <sys/stat.h>

#include <dirent.h>
#include <unistd.h>

#include "class_tax.h"

using namespace std;

// Parameters Def
string Out_path = "Result_Plot";
string Html_file = "taxonomy.html";
string Txt_file = "taxonomy.txt";
string Svg_file = "taxonomy.svg";

string database = "gg_13";
string database_path;

int Depth = 7;

string Listfile;
string Listprefix;

vector <string> Infilename;
vector <string> Sam_name;

int Mode = -1;

int printhelp(){
    
    cout << "Class-tax version " << Version << endl;
    cout << "Usage : class-tax [Options] Value" << endl;
    cout << "Options : " << endl;
    cout << "\t-i Input filename [Conflict with -l]" << endl;
    cout << "\t-l Input filename list [Conflict with -i]" << endl;
    cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
    
    cout << "\t-o Output Path, default is \"Result_Plot\"" << endl;
    cout << "\t-h Help" << endl;
    
    exit(0);
    return 0;
    }

int Parse_Para(int argc, char * argv[]){
        
    database_path = Check_Env() + "/databases/" + database + "/";
    
    if (argc ==1) 
		printhelp();
    
    Mode = -1;
    Listprefix = "";
    
    int i = 1;    
    int sam_num = 0;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            case 'i':                                  
                                 if (Mode != -1) {
                                             
                                             cerr << "Error: -i conflicts with -l" << endl;
                                             exit(0);
                                             
                                             }
                                  Infilename.push_back(argv[i+1]);
                                  Sam_name.push_back("Sample");
                                  Mode = 0;                                      
                                  break;
                                                      
                            case 'l':                                      
                                      if (sam_num != 0){
                                                  
                                                  cerr << "Error: -l conflicts with -i" << endl;
                                                  exit(0);
                                                  
                                                  }
                                      
                                      Listfile = argv[i+1];
                                      Mode = 1;
                                      break;                                      
                            
                            case 'p': Listprefix = argv[i+1]; break;               
                            case 'o': Out_path = argv[i+1]; break;                            
                            case 'd': Depth = atoi(argv[i+1]); break;                            
                            case 'h': printhelp(); break;

                            default : printf("Error: Unrec argument %s\n", argv[i]); printhelp(); break; 
                            }
         i+=2;
         }
    
    Html_file = Out_path + "/" + Html_file;
    Txt_file = Out_path + "/" + Txt_file;
    Svg_file = Out_path + "/" + Svg_file;
    
    if (Depth <=0) Depth = 7;
    
    if (Mode == 0)
       sam_num = 1;
    else 
         sam_num = Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix);         
    
    Check_Path(Out_path.c_str(), 1);

    return sam_num;
    }

int Copy_JS(const char * path){
    
    string js_copy ="cp "; 

    js_copy += Check_Env();
    js_copy += "/html/* ";
    js_copy += path;
    
    system(js_copy.c_str());
    
    return 0;
    
    }


int main(int argc, char * argv[]){
    
    int sam_num = Parse_Para(argc, argv);
    
    TNode::Init(database_path, sam_num);
                
    TNode root;
    
    cout << endl << "Taxonomic Classification Starts" << endl;
        
    for (int i = 0; i < sam_num; i++){
        
        cout << root.Read_file(Infilename[i].c_str(), Sam_name[i], i) << " sequences are loaded" << endl;
    }
    
    if (sam_num > 1)
    
       root.Out_Tree_SVG(Svg_file.c_str(), Depth);
    
    root.Out_Tree_Html(Html_file.c_str(), Depth);
    //root.Out_Tree_Distribution(txt_file.c_str(), depth);
    
    Copy_JS(Out_path.c_str());
    
    cout << endl << "Taxonomic Classification Finished"<< endl;
    cout << "Plot Finished"<< endl;

    //system("pause");
    
    return 0;
    
    }

