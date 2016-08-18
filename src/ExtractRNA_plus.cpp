// Updated at Jan 14, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <stdlib.h>

#include <sys/dir.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>

#include "ExtractRNA.h"
#include "fastq.h"
#include "utility.h"
#include "version.h"


using namespace std;

int Print_ExtractRNA_Help(){
    
    cout << "ExtractRNA version " << Version << endl;
    cout << "Usage : extract-rna [Options] Value" << endl;
    cout << "\t-m Input single sequence file (Shotgun)" <<endl;
    cout << "\t   To appoint the inputfile metagenomic shotgun file." << endl;
    cout << "\t-o Output Path" <<endl;
    cout << "\t   To appoint the result output path, default is \"Extract_RNA\"." << endl;
   	cout << "\t-l integer value" << endl;
    cout << "\t   To assign the 16S rRNA length threshold of rRNA extraction, conflict with -r. 0 is disabled, default is 0." << endl;
    //cout << "\t-b 'B' or 'E'" << endl;
    //cout << "\t   To appoint the domain of the sample, B is for Bacteria(16S rRNA) and E is for Eukaryota(18S rRNA), default is B" << endl;
    cout << "\t-t Integer value" << endl;
	cout << "\t   To assign the core number of CPU. Default is auto." << endl;
    cout << "\t-h Help" << endl;
    exit(0);
    return 0;
    }

int main(int argc, char * argv[]){
    
    string this_path = Check_Env();
    
    //string modelfilename = "$ParallelMETA/models/bac_ssu.hmm" ;
    int boundary = 0;
    string infilename = "meta.fasta";
    string out_path = "Extract_RNA";
    
    int length_filter = 0;
    int format = 0; //0:fasta 1:fastq
    
    int coren = 0;
    
    if (argc ==1) 
		Print_ExtractRNA_Help();
    
    int i = 1;
    
    while(i<argc){
         if (argv[i][0] != '-') {
                           printf("Argument # %d Error : Arguments must start with -\n", i);
                           exit(0);
                           };           
         switch(argv[i][1]){
                            
                            case 'm' : infilename = argv[i+1]; break; 
                            case 'o' : out_path = argv[i+1]; break; 
                            case 'l' : length_filter = atoi(argv[i+1]); break;
                            //case 'b' : if ((argv[i+1][0] == 'E') || (argv[i+1][0] == 'e')) boundary = 1; break;
                            case 't' : coren = atoi(argv[i+1]); break;
                            case 'h' : Print_ExtractRNA_Help(); break;

                            default : printf("Error: Unrec argument %s\n", argv[i]); Print_ExtractRNA_Help(); break; 
                            }
         i+=2;
         }
    
    int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
    if ((coren <= 0) || (coren > Max_core_number)){
                    //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
                    coren = Max_core_number;
                    }
    
    format = Check_Format(infilename.c_str());
     
    Check_Path(out_path.c_str(), 0);
    Check_Path((out_path + "/tmp").c_str(), 0);
    
    if (format == 1){//If fastq
                  
        string tempfilename = out_path + "/meta.fasta";
        cout << "Pre-computation for Fastq Starts" << endl;
        cout << endl << Fastq_2_Fasta(infilename.c_str(), tempfilename.c_str()) << " sequences have been pre-computed" << endl << endl;
        infilename = tempfilename;
        }
    
    ExtractRNA(boundary, infilename, out_path, length_filter, this_path, coren);
    
    string command = "rm -rf " + out_path + "/tmp";
    system(command.c_str());
    

    return 0;
    
    }
