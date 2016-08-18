// Updated at Nov 25, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <ctime>
#include <omp.h>
#include <math.h>
#include <unistd.h>

#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdlib.h>

#include "utility.h"
#include "version.h"

using namespace std;

string Infilename;
string Outdirectory = "result";
string Resultname = "out"; 
int Bot=20;
int Threadnum=0;
int line_max=0;
bool Annotation = false;

int printhelp(){
    cout << "Rare-curv Version : " << Version << endl;
    cout << "Usage : " << endl;
    cout << "rare-curv [-option] value" << endl;
    cout << "\toption : " << endl;

	cout << "\t-i or -T Input feature table with Absolute Count (*.Count) [Required]" << endl;
	cout << "\t-o Output file directory, default is \"result\"" << endl;
	cout << "\t-p Prefix name of output, default is \"out\""<<endl;
	cout << "\t-l Rarefaction curve label, T is enable and F is disable, default is F"<<endl;
	cout << "\t-b The bootstrap value, default is 20" << endl;
	cout << "\t-t Cpu core number, default is auto" <<endl;	
	cout << "\t-h Help" << endl;

	exit(0);

	return 0;

};

void Parse_Para(int argc, char * argv[]){
	if(argc==1){
		printhelp();
	}

	int i=1;
	while(i<argc){
		if(argv[i][0]!='-'){
			cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
			exit(0);
		}
		switch(argv[i][1]){
		case 'i': 
        case 'T': Infilename = argv[i+1];break;
		case 'o': Outdirectory = argv[i+1];break;
        case 'p': Resultname = argv[i+1];break;
		case 't': Threadnum = atoi(argv[i+1]);break;
		case 'b': Bot = atoi(argv[i+1]);break;
		case 'l': if ((argv[i+1][0] == 't') || (argv[i+1][0] == 'T')) Annotation = false; break;
		default : cerr << "Error: Unrec argument " << argv[i] << endl;
			printhelp();
			break;
		}
		i=i+2;

	}

   int Max_core_number = sysconf(_SC_NPROCESSORS_CONF);
    
   if ((Threadnum <= 0) || (Threadnum > Max_core_number)){
                    Threadnum = Max_core_number;
                    }  
   
   if (Bot <= 0){
         cerr << "Warning: The bootstrap value should be larger than 0, change to default (20)";
         Bot = 20;
         }
                    
} 

int main(int argc, char *argv[])
{
	Parse_Para(argc, argv);
	Check_Path(Outdirectory.c_str(), 1);
	ifstream infile(Infilename.c_str(),ios::in);
	if(!infile)
	{
		cerr<<"Error: Cannot open input file:  "<<Infilename<<endl;
		return 0;
	}
	string Outfilename_Shannon=Outdirectory+"/Shannon.txt";
	string Outfilename_Shannon_max=Outdirectory+"/Shannon_max.txt";
	string Outfilename_Observe_otu=Outdirectory+"/Observe_otu.txt";
	string Outfilename_Observe_otu_max=Outdirectory+"/Observe_otu_max.txt";
	ofstream Outfile_Shannon(Outfilename_Shannon.c_str(),ios::out);
	if (!Outfile_Shannon)
	{
		cerr<<"Error: Cannot open file: "<<Outfilename_Shannon<<endl;
		return 0;
	}
	Outfile_Shannon<<"SampleID\tx\ty"<<endl;
	ofstream Outfile_Observe_otu(Outfilename_Observe_otu.c_str(),ios::out);
	if (!Outfile_Observe_otu)
	{
		cerr<<"Error: Cannot open file: "<<Outfilename_Observe_otu<<endl;
		return 0;
	}
	Outfile_Observe_otu<<"SampleID\tx\ty"<<endl;
	fstream Outfile_Shannon_max(Outfilename_Shannon_max.c_str(),ios::out);
	if (!Outfile_Shannon_max)
	{
		cerr<<"Error: Cannot open file: " << Outfilename_Shannon_max<<endl;
		return 0;
	}
	Outfile_Shannon_max<<"SampleID\txm\tym"<<endl;
	ofstream Outfile_Observe_otu_max(Outfilename_Observe_otu_max.c_str(),ios::out);
	if (!Outfile_Observe_otu_max)
	{
		cerr<<"Error: Cannot open file: " <<Outfilename_Observe_otu_max<<endl;
		return 0;
	}
	Outfile_Observe_otu_max<<"SampleID\txm\tym"<<endl;
	string buffer;
	string Otu_id;
	int line_count=0;
	vector<string> Sample;
	vector< vector<int> > Otu_num;
	vector<string> Otu_name;
	vector< vector<string> > Otu_table;
	map< string,vector<int> > map_table;
	//vector< vector<double> > Observe_otu;
	//vector< vector<double> > Shannon_vec;
	getline(infile,buffer);
	stringstream strin(buffer);
	strin>>Otu_id;
	while(strin>>Otu_id)
	{
		Otu_name.push_back(Otu_id);
	}
	while (getline(infile,buffer))
	{
		vector<int> line_taxa;
		vector<string> Otu_count;
		stringstream Str_Otu_count(buffer);
		string SampleID;
		Str_Otu_count>>SampleID;
		Sample.push_back(SampleID);
		string Taxa_count;
		int i=0;
		while(Str_Otu_count>>Taxa_count)
		{
			if (atoi(Taxa_count.c_str())!=0)
			{
				for (int j=0;j<atoi(Taxa_count.c_str());j++)
				{
					Otu_count.push_back(Otu_name[i]);
				}
			} 
			i++;
		}
		Otu_table.push_back(Otu_count);
	}
	
	//cout << "Loaded" << endl; //debug
    //cout << Threadnum << endl;//debug

    int line_max=Otu_table[0].size();
	 for(int i=0;i<Otu_table.size();i++)
        {
                if(Otu_table[i].size()>=line_max)
                        line_max=Otu_table[i].size();
        }
	int **Observe_otu=new int*[Otu_table.size()];
	double **Shannon_vec=new double*[Otu_table.size()];
	for(int i=0;i<Otu_table.size();i++)
	{
		Observe_otu[i]=new int[line_max];
		Shannon_vec[i]=new double[line_max];
	}
	
    omp_set_num_threads(Threadnum);
    #pragma omp parallel for schedule(dynamic, 1)
    //#pragma omp parallel for num_threads(Threadnum)	
	for (int i=0;i<Otu_table.size();i++)
	{
		int times=Otu_table[i].size();
		for(int k=0;k<Bot;k++)
		{
			vector<int> line_taxa;
			vector<double> line_shannon;
			int ran_num;
			string temp;
			for(int j=0;j<times;j++)
			{
				srand((unsigned)time(0)*unsigned(j)*unsigned(i+1)*unsigned(k+1));
				ran_num=rand()%times;
				swap(Otu_table[i][j],Otu_table[i][ran_num]);
			}
			set<string> Taxa;
			map<string,int> Otu_map;
			vector<int> Otu_map_vec;
			map<string,int> Map_shannon;
			for (int j=0;j<times;j++)
			{
				Taxa.insert(Otu_table[i][j]);
				if (Map_shannon.find(Otu_table[i][j])!=Map_shannon.end())
				{
					Map_shannon[Otu_table[i][j]]++;
				}
				else
					Map_shannon[Otu_table[i][j]]=1;

				//if ((j+1)%step==0&&(j+1)/step>=1)
				//{
				    double Shannon=0.0;
					line_taxa.push_back(Taxa.size());
					for (map<string,int>::iterator Map_shannon_it=Map_shannon.begin();Map_shannon_it!=Map_shannon.end();Map_shannon_it++)
					{
						Shannon+=(float(Map_shannon_it->second)/(float(j+1)))*(log((long double)(float(Map_shannon_it->second)/float(j+1))));
					}
					line_shannon.push_back(Shannon);
				//}
				
			}
			if (k==0)
			{
				for(int ii=0;ii<line_taxa.size();ii++)
				{
				        Observe_otu[i][ii]=line_taxa[ii];
				        Shannon_vec[i][ii]=line_shannon[ii];
                }
			}
			else
			{
				for(int j=0;j<line_taxa.size();j++)
				{
					Observe_otu[i][j]+=line_taxa[j];
					Shannon_vec[i][j]+=line_shannon[j];
				}

			}
				
		}
	}
		
		for(int j=0;j<Sample.size();j++){
			for (int i=0;i<Otu_table[j].size();i++)
			{
				Outfile_Observe_otu<<Sample[j]<<"\t"<<i<<"\t"<<float(Observe_otu[j][i])/float(Bot)<<endl;
				Outfile_Shannon<<Sample[j]<<"\t"<<i<<"\t"<<-float(Shannon_vec[j][i])/float(Bot)<<endl;
			}
			Outfile_Observe_otu_max<<Sample[j]<<"\t"<<Otu_table[j].size()<<"\t"<<float(Observe_otu[j][Otu_table[j].size()-1])/float(Bot)<<endl;
			Outfile_Shannon_max<<Sample[j]<<"\t"<<Otu_table[j].size()<<"\t"<<-float(Shannon_vec[j][Otu_table[j].size()-1])/float(Bot)<<endl;
		}
	Outfile_Observe_otu.close();
	Outfile_Observe_otu_max.close();
	Outfile_Shannon.close();
	Outfile_Shannon_max.close();

	string Radd=Check_Env()+"/Rscript/rarefaction.R ";
	string Opt;
	if(Annotation)
	 Opt="Rscript "+Radd+"-o "+Outdirectory+" -p "+Resultname+" -a T";
	else
	 Opt="Rscript "+Radd+"-o "+Outdirectory+" -p "+Resultname+" -a F";
	system(Opt.c_str()); 
		return 0;
	//system("pause");
}
