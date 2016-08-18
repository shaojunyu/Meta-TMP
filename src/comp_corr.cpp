// Updated at Nov 26, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
//version 3.1 or above with _Table_Format

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>

#include "version.h"
#include "utility.h"
#include "table_format_corr.h"

using namespace std;
string Infilename;
string Infilenamecorr;
string Outfilenameself = "self_corr_matrix.out";
string Outfilenamecorr = "corr_matrix.out";
string OutCytoscape_filename = "Cytoscape.out";
int Coeff = 0; //0: S 1: P
string Selotu;
int Flag = 0;
int Tn = 0;
float co_occurrence = 0.8;//Correlation based co-occurrence
bool is_output_cytoscape = true;
bool Is_network = false;
float Network_t = 0.5;

void Print_Help() {

    cout << "Comp-corr version " << Version << endl;
    cout << "Usage : comp-corr [Options] Value" << endl;
    cout << "Options : " << endl;
    cout << "\t-i Feature (Taxa/OTU/Function etc.) table file [Required]" << endl;
    cout << "\t-m Correlation file name" << endl;
    cout << "\t-o Output prefix, default is \"corr_matrix\"" << endl;
    cout << "\t-c Selected feature, separated by \",\"" << endl;
    cout << "\t-f 0:(Spearman) or 1:(Pearson) metrics,default is 0" << endl;
    cout << "\t-C Correlation based co-occurrence, default is 0.8 and 0 means no co-occurrence" <<endl;
    cout << "\t-s output Cytoscape data, default T(rue) or F(alse), default is T"<< endl;
    cout << "\t-N Network based co-occurrence analysis, T(rue) or F(alse), default is F" << endl;
    cout << "\t-T Netowrk analysis threshold, and Cytoscape output threshold, default is 0.5" << endl;
    cout << "\t-t Cpu core number, default is auto" << endl;
    cout << "\t-h Help" << endl;

    exit(0);
}

void Para(int argc, char *argv[]) {
    int i = 1;
    if (argc == 1)
        Print_Help();
    while (i < argc) {
        if (argv[i][0] != '-') {
            cout << "Argument #" << i << " Error : Arguments must start with -\n" << endl;
            exit(0);
        }
        switch (argv[i][1]) {
        case 'i':
            Infilename = argv[i + 1];
            break;
        case 'm':
            Infilenamecorr = argv[i + 1];
            break;
        case 'o': {
            Outfilenameself = argv[i + 1];
            Outfilenameself += ".self_matrix.out";
            Outfilenamecorr = argv[i + 1];
            Outfilenamecorr += ".corr_matrix.out";
            OutCytoscape_filename = argv[i + 1];
            OutCytoscape_filename += ".Cytoscape.out";
            break;
        }
        case 'f':
            Coeff = atoi(argv[i + 1]);
            break;
        case 'C':
            co_occurrence = atof(argv[i + 1]);
            if (co_occurrence > 1 || co_occurrence < 0)
            {
                cout << "Argument #" << i << " Error : co_occurrence must between 0 and 1" << endl;
                exit(0);
            }
            break;
        case 's':
            if ((argv[i + 1][0] == 't') || (argv[i + 1][0] == 'T'))
                is_output_cytoscape = true;
            break;
        case 'c': {
            Selotu = argv[i + 1];
            Flag = 1;
            break;
        }
        case 'N':
            if ((argv[i + 1][0] == 't') || (argv[i + 1][0] == 'T'))
                Is_network = true;
            break;
        case 'T': {
            Network_t = atof(argv[i + 1]);
            break;
        }
        case 't': {
            Tn = atoi(argv[i + 1]);
            break;
        }
        default:
            cout << "Error: Unrec arguments" << argv[i] << endl;
            Print_Help();
            break;
        }
        i += 2;
    }

    int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);

    if ((Tn <= 0) || (Tn > Max_core_number)) {
        //cerr << "Core number must be larger than 0, change to automatic mode" << endl;
        Tn = Max_core_number;
    }
}

int main(int argc, char *argv[]) {

    Para(argc, argv);

    _Table_Format taxa(Infilename.c_str(), 0);
    if (co_occurrence != 0)
    {
        taxa.Calc_Corr_Matrix_based_on_co_occurrence(Outfilenameself.c_str(), Coeff, Tn, co_occurrence, is_output_cytoscape, Network_t);
    }else{
        taxa.Calc_Corr_Matrix(Outfilenameself.c_str(), Coeff, Tn, is_output_cytoscape, Network_t);
    }
    
    // if (is_output_cytoscape)
    // {
    //     cout<<"adasd"<<endl;
    //     taxa.Output_Cytoscape(OutCytoscape_filename.c_str(), Network_t);
    // }
    // cout <<"dadasdadada" <<endl;
    if (Infilenamecorr.size() != 0) {
        _Table_Format_Meta Corr_With_Meta;
        if (Flag == 0)
            Corr_With_Meta.Calc_Corr_Meta_Matrix(Infilename.c_str(), Infilenamecorr.c_str(), Outfilenamecorr.c_str(), Coeff, Tn);
        else
            Corr_With_Meta.Calc_Corr_Meta_Matrix_S(Infilename.c_str(), Infilenamecorr.c_str(), Outfilenamecorr.c_str(), Selotu, Coeff, Tn);
    }

    if (Is_network) {
        char command[BUFFER_SIZE];
        sprintf(command, "Rscript %s/Rscript/PM_Network.R -i %s -o %s -t %f", Check_Env().c_str(), Outfilenameself.c_str(), (Outfilenameself + ".network.pdf").c_str(), Network_t);
        system(command);

        sprintf(command, "Rscript %s/Rscript/PM_MCL.R -i %s -o %s", Check_Env().c_str(), (Outfilenameself + "_to_Cytoscape.out").c_str(), (Outfilenameself + ".MCL_network.html").c_str());
        system(command);

        sprintf(command, "Rscript %s/Rscript/PM_MCODE.R -i %s -o %s", Check_Env().c_str(), (Outfilenameself + "_to_Cytoscape.out").c_str(), (Outfilenameself + ".MCODE_network.html").c_str());
        system(command);

    }

    return 0;
}
