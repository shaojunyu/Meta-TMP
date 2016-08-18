// Updated at June 3, 2016
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include "pipeline.h"

using namespace std;

int main(int argc, char * argv[]){
    
    Parse_Para(argc, argv);
    
    char command[BUFFER_SIZE];
    //Parallel for R
    vector <string> command_parallel_scripts; 
    vector <string> command_parallel_titles;
            
    Check_Path(Out_path.c_str(), 1);
        
    //open script file
    ofstream outscript((Out_path + "/scripts.sh").c_str(), ofstream::out);
    if (!outscript){
                    cerr << "Warning: Cannot open output file : " << Out_path << "/scripts.sh" << endl;                    
                    }
    
    cout << "Parallel-META Pipeline Version: " << Version << endl;
    outscript << "#Parallel-META Pipeline Version: " << Version << endl;
    
    //check metadata
    if (Load_ID(Meta_file.c_str(), Ids) == 0){
        cerr << "Error: Please check the Meta data file (-m): at least contains 1 columns of Sample ID" << endl;
        return 0;
        }
    
    if (!Check_Ids(Ids)){
        cerr << "Error: Sample ID cannot be started by number" << endl;
        return 0;
    }
    
    Input_sam_num = Ids.size() -  1;
    
    switch (Step){
        
    //Step 0: Parallel-META
    case 0:
               Check_Path(Singlesample_dir.c_str(), 1);
               Check_Path(Singlesamplelist_dir.c_str(), 1);
         
               
               if (Load_List(Seq_list_file.c_str(), Seq_files, List_prefix) == 0){
                                              cerr << "Error: Please check the sequence list file (-i) or the list path prefix (-p)" << endl;
                                              return 0;
                                              };
               
               //is_pair_end
               if (Seq_files.size() == (Ids.size() - 1) * 2){
                                    Is_paired_seq = true;
                                    if (Seq_type != 'r'){
                                                 cerr << "Error: Pair-end sequences only support 16S rRNA (-M F)" << endl;
                                                 return 0;
                                                 }
                                    }
               else if (Seq_files.size() != (Ids.size() - 1)){
                                    cerr << "Error: Sequence files (pairs) and meta data should have the same sample number and order" << endl;
                                    //cout << Seq_files.size() << "\t" << Ids.size() << endl; //debug
                                    return 0;
                                    }
            
               //format_seq
               /*
               if (Is_format_check){
               cout << endl << "Sequence format check" << endl;    
               outscript << endl << "#Sequence format check" << endl;                             
               
               if (List_prefix.size() > 0) // add prefix
                  sprintf(command, "%s/format-seq -l %s -p %s", Bin_path.c_str(), Seq_list_file.c_str(), List_prefix.c_str());          
               else
                  sprintf(command, "%s/format-seq -l %s", Bin_path.c_str(), Seq_list_file.c_str());                               
               Run_With_Error(command, "format-seq", Error_file.c_str());
               outscript << command << endl;
               }
               */
               
               //profiling
               cout << endl << "Microbial Community profiling" << endl;    
               outscript << endl << "#Microbial Community profiling" << endl;   
               
               //taxa               
               for (int i = 0; i < Ids.size() - 1; i ++){
                               
                               cout << endl << "Processing sample " << i + 1 << " of " << Ids.size() - 1 << endl;   
                               if (Is_paired_seq) //pair -end                            
                                  sprintf(command, "%s/parallel-meta -r %s -R %s -o %s -t %d -e %s -f F -P %s -k %c", Bin_path.c_str(), Seq_files[i * 2].c_str(), Seq_files[i * 2 + 1].c_str(), (Singlesample_dir + "/" + Ids[i+1]).c_str(), Coren, Align_mode.c_str(), Paired_mode.c_str(), Is_format_check);
                               else
                                   sprintf(command, "%s/parallel-meta -%c %s -o %s -t %d -e %s -f F -L %d -k %c", Bin_path.c_str(), Seq_type, Seq_files[i].c_str(), (Singlesample_dir + "/" + Ids[i+1]).c_str(), Coren, Align_mode.c_str(), Length_t, Is_format_check);
                               Run_With_Error(command, "parallel-meta", Error_file.c_str());
                               outscript << command << endl;
                               }                              
              //taxa list
              Taxa_list_file = Singlesamplelist_dir + "/taxa.list";
              Make_list(Taxa_list_file.c_str(), Singlesample_dir.c_str(), Ids, 0);
              
              List_prefix = "";
              //Step 0 finished
            
    //Step 1: chek list
    case 1: //check list first
           if (Taxa_list_file.size() == 0){
                              cerr << "Error: Please check the taxa list (-l)" << endl;
                              return 0;
                              }
                                 
           //If prefix
           if (List_prefix.size() > 0){ // Add prefix to list
                                  Check_Path(Singlesamplelist_dir.c_str(), 1);
                                  Add_list_prefix(Taxa_list_file.c_str(), List_prefix.c_str(), (Singlesamplelist_dir + "/taxa.list").c_str());
                                  Taxa_list_file = Singlesamplelist_dir + "/taxa.list";
                                  }
           //Ids reload from taxa list file
           Ids.clear();
           Ids.push_back("Sample");
           Load_ID(Taxa_list_file.c_str(), Ids);
            
           //check output dir
           Check_Path(Abd_dir.c_str(), 1);
           Check_Path(Dist_dir.c_str(), 1);
           Check_Path(Clust_dir.c_str(), 1);
           Check_Path(Marker_dir.c_str(), 1);
           Check_Path(Network_dir.c_str(), 1);
           Check_Path(Alpha_dir.c_str(), 1);
           Check_Path(Beta_dir.c_str(), 1);
           Check_Path(Sampleview_dir.c_str(), 1);
           Check_Path(Temp_dir.c_str(), 1);
    
    //Step 2: Rarefaction & Normalization
    if (Is_rare){
            cout << endl << "Sequence Normalization" << endl;
            outscript << endl << "#Sequence Normalization" << endl;
            //rand-rare
            
            Check_Path(Singlesamplerare_dir.c_str(), 1);
            
            sprintf(command, "%s/rand-rare -l %s -o %s -b %d -s %d", Bin_path.c_str(), Taxa_list_file.c_str(), Singlesamplerare_dir.c_str(), Bootstrap, Rare_depth);
            Run_With_Error(command, "rand-rare", Error_file.c_str());
            outscript << command << endl;
            
            //Update ids
            Check_Ids_By_Path(Singlesamplerare_dir.c_str(), Ids);
        
            //taxa list rare
            Check_Path(Singlesamplelist_dir.c_str(), 1);
            Taxa_list_file = Singlesamplelist_dir + "/taxa_rare.list";
            Make_list(Taxa_list_file.c_str(), Singlesamplerare_dir.c_str(), Ids, 3); //classification.rare.txt
                
            Singlesample_dir = Singlesamplerare_dir;
            }
    
    //Update metadata
    if (!(Check_Metadata_By_Ids(Ids, Meta_file, Out_path + "/meta.txt")))
                return 0;
    
    Filter_sam_num = Ids.size() - 1;
    
    if (Filter_sam_num < MIN_SAM_NUM){
        cerr << "Error: " << Filter_sam_num << " sample(s) passed the filtering." << MIN_SAM_NUM << endl;
        cerr << "The required minimum sample number is " << MIN_SAM_NUM << endl;
        return 0;
            }
            
    //Step 3: Function Prediction
    if (Is_func){
            cout << endl << "Function Prediction" << endl;
            outscript << endl << "#Function Prediction" << endl;
            
            //Parse func
            Check_Path(Singlesample_dir.c_str(), 1);
            sprintf(command, "%s/class-func -l %s -o %s -t %d", Bin_path.c_str(), Taxa_list_file.c_str(), Singlesample_dir.c_str(), Coren);
            Run_With_Error(command, "class-func", Error_file.c_str());
            outscript << command << endl;
                
            //func list
            Check_Path(Singlesamplelist_dir.c_str(), 1);
            Func_list_file = Singlesamplelist_dir + "/func.list";
            Make_list(Func_list_file.c_str(), Singlesample_dir.c_str(), Ids, 1);
            }
    
    //Step 3: Sample visualization
    cout << endl << "Sample visualization" << endl;
    outscript << endl << "#Sample visualization" << endl;
          //class-tax
          if (Is_taxa){
          sprintf(command, "%s/class-tax -l %s -o %s/", Bin_path.c_str(), Taxa_list_file.c_str(), Sampleview_dir.c_str());
          Run_With_Error(command, "class-tax", Error_file.c_str());
          outscript << command << endl;    
          }          
    
    //Step 4: Feature select
    cout << endl << "Feature Selection" << endl;
    outscript << endl << "#Feature Selection" << endl;
           //taxa-sel
           if (Is_taxa){              
           for (int i = 0; i < TLevN; i ++) // no OTU
               if (TLevel_Set[i]){             
               sprintf(command, "%s/taxa-sel -l %s -o %s/taxa.%s -L %d -m 0 -n 0 -P T", Bin_path.c_str(), Taxa_list_file.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), i + 1);
               Run_With_Error(command, "taxa-sel", Error_file.c_str());
               outscript << command << endl;
               }
           }
           //OTU table
           if (!TLevel_Set[6])
                              {
                               sprintf(command, "%s/taxa-sel -l %s -o %s/taxa.%s -L %d -m 0 -n 0 -z 0 -v 0", Bin_path.c_str(), Taxa_list_file.c_str(), Abd_dir.c_str(), TLevel[6].c_str(), 7);
                               Run_With_Error(command, "taxa-sel", Error_file.c_str());
                               outscript << command << endl;
                               //TLevel_Set[6] = false; //disable OTU
                               }
                                 
           //func-sel
           if (Is_func){          
           for (int i = 0; i < FLevN - 1; i ++)
               if (FLevel_Set[i]){
               sprintf(command, "%s/func-sel -l %s -o %s/func.%s -L %d -P T", Bin_path.c_str(), Func_list_file.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), i + 1);
               Run_With_Error(command, "func-sel", Error_file.c_str());
               outscript << command << endl;
               }
           //KO table
           if (!FLevel_Set[4]) {
                               sprintf(command, "%s/func-sel -l %s -o %s/func.%s -L %d", Bin_path.c_str(), Func_list_file.c_str(), Abd_dir.c_str(), FLevel[3].c_str(), 4);
                               Run_With_Error(command, "func-sel", Error_file.c_str());
                               outscript << command << endl;
                               }
           //calc func-nsti           
            sprintf(command, "%s/class-func-nsti -l %s -o %s", Bin_path.c_str(), Taxa_list_file.c_str(), (Abd_dir + "/func.nsti").c_str());                               
            Run_With_Error(command, "class-func-nsti", Error_file.c_str());
            outscript << command << endl;           
           }
           
    //Step 5: dist
    cout << endl << "Dist Calculation" << endl;
    outscript << endl << "#Dist Calculation" << endl;
           //calc taxa dist
           if (Is_taxa){
           sprintf(command, "%s/comp-sam -T %s/taxa.OTU.Abd -o %s/taxa.dist -t %d -d T -P T -c %d", Bin_path.c_str(), Abd_dir.c_str(), Dist_dir.c_str(), Coren, Cluster);
           Run_With_Error(command, "comp-sam", Error_file.c_str());
           outscript << command << endl;
           }
           
           if (Is_func){      
           //calc func dist
           sprintf(command, "%s/comp-sam-func -T %s/func.KO.Abd -o %s/func.dist -t %d -d T -P T -c %d", Bin_path.c_str(), Abd_dir.c_str(), Dist_dir.c_str(), Coren, Cluster);
           Run_With_Error(command, "comp-sam-func", Error_file.c_str());
           outscript << command << endl;
           }
           
    
           
    //Step 6: Co-occur network
    cout << endl << "Correlation Calculation" << endl;
    outscript << endl << "#Correlation Calculation" << endl;
           //taxa-corr                 
           if (Is_taxa){
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){               
               sprintf(command, "%s/comp-corr -i %s/taxa.%s.Abd -o %s/taxa.%s -f 1 -N T -T %f", Bin_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Network_dir.c_str(), TLevel[i].c_str(), Network_t);
               Run_With_Error(command, "comp-corr", Error_file.c_str());
               outscript << command << endl;
               }
           }
           
    //Step 7: Rarefaction curve
    //might be slow
    if (Is_rare_curve){
    cout << endl << "Rarefaction Analysis" << endl;
    outscript << endl << "#Rarefaction Analysis" << endl;
           //OTU based
           if (Is_taxa){            
           sprintf(command, "%s/rare-curv -i %s/taxa.OTU.Count -o %s -p taxa.OTU -b 10", Bin_path.c_str(), Abd_dir.c_str(), Alpha_dir.c_str());
           Run_With_Error(command, "rare-curv", Error_file.c_str());
           outscript << command << endl;
           }
    }     
    
    cout << endl << "Statistical Analysis" << endl;  
    //Step 8: PCoA
    //cout << endl << "PCoA Calculation" << endl;
    outscript << endl << "#PCoA Calculation" << endl;
           //taxa PCoA
           if (Is_taxa){
           sprintf(command, "Rscript %s/PM_Pcoa.R -m %s -d %s/taxa.dist -o %s/taxa.pcoa.pdf", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Clust_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Pcoa.R");
           outscript << command << endl;
           }
           
           //func PCoA
           if (Is_func){
           sprintf(command, "Rscript %s/PM_Pcoa.R -m %s -d %s/func.dist -o %s/func.pcoa.pdf", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Clust_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Pcoa.R");
           outscript << command << endl;
           }
           
    //Step 9: PCA
    //cout << endl << "PCA Calculation" << endl;
    outscript << endl << "#PCA Calculation" << endl;
            
           //taxa PCA
           if (Is_taxa){
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){           
               sprintf(command, "Rscript %s/PM_Pca.R -m %s -i %s/taxa.%s.Abd -o %s/taxa.%s.pca.pdf", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Clust_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Pca.R");
               outscript << command << endl;
               }
           }
           
           //func PCA
           if (Is_func){
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){                        
                sprintf(command, "Rscript %s/PM_Pca.R -m %s -i %s/func.%s.Abd -o %s/func.%s.pca.pdf", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Clust_dir.c_str(), FLevel[i].c_str());
                command_parallel_scripts.push_back(command);
                command_parallel_titles.push_back("PM_Pca.R");
                outscript << command << endl;
               }
           }
           
    //Step 10: Multivariate Statistical Analysis
    //cout << endl << "Multivariate Statistical Analysis" << endl;
    outscript << endl << "#Multivariate Statistical Analysis" << endl;           
    
           //Parallel  
           //Alpha: taxa         
           if (Is_taxa){                        
               for (int i = 0; i < TLevN; i ++) //contain OTU
               if ((TLevel_Set[i]) || (i == TLevN - 1)){ // OTU                  
               sprintf(command, "Rscript %s/PM_Adiversity.R -m %s -i %s/taxa.%s.Abd -o %s -p taxa.%s", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Alpha_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Adiversity.R");
               outscript << command << endl;
               }                      
            
           //Beta: taxa dist
           sprintf(command, "Rscript %s/PM_Bdiversity.R -m %s -d %s/taxa.dist -o %s -p taxa.dist -n Meta-Storms", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Beta_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Bdiversity.R");
           outscript << command << endl;               
           }
           
           //Alpha: func
           if (Is_func){
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){               
               sprintf(command, "Rscript %s/PM_Adiversity.R -m %s -i %s/func.%s.Abd -o %s -p func.%s", R_path.c_str(), Meta_file.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Alpha_dir.c_str(), FLevel[i].c_str()); 
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Adiversity.R");
               outscript << command << endl;           
               }                     
                       
           //Beta: func dist
           sprintf(command, "Rscript %s/PM_Bdiversity.R -m %s -d %s/func.dist -o %s -p func.dist -n Meta-Storms-func", R_path.c_str(), Meta_file.c_str(), Dist_dir.c_str(), Beta_dir.c_str());
           command_parallel_scripts.push_back(command);
           command_parallel_titles.push_back("PM_Bdiversity.R");
           outscript << command << endl;           
           }
                   
    //Step 11: Test
    //cout << endl << "Marker Analysis" << endl;
    outscript << endl << "#Marker Analysis" << endl;
    
           //Parallel
           //taxa-test
           if (Is_taxa){ 
           //test          
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){          
               sprintf(command, "Rscript %s/PM_Marker.R -i %s/taxa.%s.Abd -m %s -o %s -p taxa.%s -P %c", R_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), TLevel[i].c_str(), Is_pair);
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker.R");
               outscript << command << endl;           
               }     
           //RF
           for (int i = 0; i < TLevN; i ++)
               if (TLevel_Set[i]){          
               sprintf(command, "Rscript %s/PM_RFscore.R -i %s/taxa.%s.Abd -m %s -o %s -p taxa.%s", R_path.c_str(), Abd_dir.c_str(), TLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), TLevel[i].c_str());
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_RFscore.R");
               outscript << command << endl;           
               }                                             
           }
           
           //func-test
           if (Is_func){   
           //test        
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){                  
               sprintf(command, "Rscript %s/PM_Marker.R -i %s/func.%s.Abd -m %s -o %s -p func.%s -P %c", R_path.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), FLevel[i].c_str(), Is_pair);
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_Marker.R");
               outscript << command << endl;
               }         
           //RF
           for (int i = 0; i < FLevN; i ++)
               if (FLevel_Set[i]){                  
               sprintf(command, "Rscript %s/PM_RFscore.R -i %s/func.%s.Abd -m %s -o %s -p func.%s", R_path.c_str(), Abd_dir.c_str(), FLevel[i].c_str(), Meta_file.c_str(), Marker_dir.c_str(), FLevel[i].c_str());               
               command_parallel_scripts.push_back(command);
               command_parallel_titles.push_back("PM_RFscore.R");
               outscript << command << endl;
               }                                           
           }
    
    //Parallel R
     omp_set_num_threads(Min(command_parallel_scripts.size(), Coren));    
     #pragma omp parallel for schedule(dynamic, 1)
     for (int i = 0; i < command_parallel_scripts.size(); i ++){
         char error_parallel[BUFFER_SIZE];
         sprintf(error_parallel, "%s/%d.log", Temp_dir.c_str(), i);
         Run_With_Error(command_parallel_scripts[i].c_str(), command_parallel_titles[i].c_str(),error_parallel);
         }
     
     //Combine error
     for (int i = 0; i < command_parallel_scripts.size(); i ++){
         sprintf(command, "cat %s/%d.log >> %s", Temp_dir.c_str(), i, Error_file.c_str());
         system(command);
         }
                                                                                    
    default: break;
    }
    
    if (rmdir(Temp_dir.c_str()) < 0){
        sprintf(command, "rm -rf %s", Temp_dir.c_str());
        system(command);
        }
    
    if (outscript){
                   outscript.close();
                   outscript.clear();
                   }
    
    Print_Report(Report_file.c_str());
    
    cout << endl << "Parallel-META Pipeline Finished" << endl;
    outscript << endl << "#Parallel-META Pipeline Finished" << endl;
    cout << "Please check the analysis results and report at " << Out_path <<endl;
    
    return 0;
    }

