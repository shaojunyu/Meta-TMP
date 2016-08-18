// Updated at May 13, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>

#include "class_func.h"

using namespace std;

// Parameters Def
string Out_path = "Results_Func/";

string database = "gg_13";
string database_path;

string Listfile;
string Listprefix;

vector<string> Infilename;
vector<string> Sam_name;
vector<string> Outfilename;

int Mode;

int Coren = 0;

int printhelp() {

  cout << "Class-func version " << Version << endl;
  cout << "Usage : class-func [Options] Value" << endl;
  cout << "Options : " << endl;
  cout << "\t-i Input single filename [Conflict with -l]" << endl;
  cout << "\t-l Input filename list list [Conflict with -i]" << endl;
  cout << "\t-p List file path prefix for '-l' [Optional]" << endl;
  cout << "\t-o Output path, default is \"Results_Func\"" << endl;
  cout << "\t-t Cpu core number, default is auto" << endl;
  cout << "\t-h Help" << endl;

  exit(0);
  return 0;
}

int Parse_Para(int argc, char *argv[]) {

  database_path = Check_Env() + "/databases/" + database + "/";

  if (argc == 1)
    printhelp();

  int i = 1;
  int sam_num = 0;
  Mode = -1;

  while (i < argc) {
    if (argv[i][0] != '-') {
      printf("Argument # %d Error : Arguments must start with -\n", i);
      exit(0);
    };
    switch (argv[i][1]) {
    case 'i':
      if (Mode != -1) {

        cerr << "Error: -i conflicts with -l" << endl;
        exit(0);
      }

      Infilename.push_back(argv[i + 1]);
      Sam_name.push_back("Sample");
      Mode = 0;
      break;

    case 'l':
      if (Mode != -1) {

        cerr << "Error: -l conflicts with -i" << endl;
        exit(0);
      }
      Listfile = argv[i + 1];
      Mode = 1;
      break;

    case 'p':
      Listprefix = argv[i + 1];
      break;
    case 'o':
      Out_path = argv[i + 1];
      break;
    case 't':
      Coren = atoi(argv[i + 1]);
      break;

    case 'h':
      printhelp();
      break;

    default:
      printf("Error: Unrec argument %s\n", argv[i]);
      printhelp();
      break;
    }
    i += 2;
  }

  Check_Path(Out_path.c_str(), 1);

  if (Mode == 0) { // single sample
    Outfilename.push_back(Out_path + "/functions.txt");
    sam_num = 1;
  }

  else if (Mode == 1) {

    sam_num = Load_List(Listfile.c_str(), Infilename, Sam_name, Listprefix);

    for (int i = 0; i < sam_num; i++) {
      Check_Path((Out_path + "/" + Sam_name[i]).c_str(), 1);
      Outfilename.push_back(Out_path + "/" + Sam_name[i] + "/functions.txt");
    }
  }

  int max_core_number = sysconf(_SC_NPROCESSORS_CONF);
  max_core_number = (max_core_number > sam_num) ? sam_num : max_core_number;

  if ((Coren <= 0) || (Coren > max_core_number)) {
    // cerr << "Core number must be larger than 0, change to automatic mode" <<
    // endl;
    Coren = max_core_number;
  }

  return sam_num;
}

int main(int argc, char *argv[]) {

  int sam_num = Parse_Para(argc, argv);

  cout << endl
       << "Functional Annotation Starts" << endl;

  _KO_OTU_Table_All KOs(database_path, sam_num, 0);

  KOs.Load_Sample_Multi(Infilename, Sam_name, Coren);

  cout << KOs.Output_Multi(Outfilename) << " KOs have been parsed out" << endl;

  cout << endl
       << "Functional Annotation Finished" << endl;

  return 0;
}
