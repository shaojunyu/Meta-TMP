// Updated at Nov 25, 2015
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <iomanip>
#include <stdlib.h>

#include "utility.h"
#include "hash.h"
#include "version.h"

using namespace std;

#ifndef _CLASS_TAX_H
#define _CLASS_TAX_H

#define Text_Length 120
#define Node_Width (Text_Length + Sam_num * 8 + 35)
#define X_Offset 54
#define Color_Count 10

string Default_color[10] =  {"255,0,0","0,0,255","0,255,0", "255,255,0", "151,173,172", "173,137,118", "255,94,72", "38,188,213", "158,157,131", "60,79,57"};

class TNode {
       public:
              TNode(){
                      Count = new float [Sam_num];
                      for (int i = 0; i< Sam_num; i++)
                          Count[i] = 0;
                     
                     Width = 0;
                     Upper = 0;
                     Lower = 0;
                     Number = 0;                
                     }
               
               TNode(string _name){
                            Count = new float [Sam_num];
                            for (int i = 0; i< Sam_num; i++)
                                  Count[i] = 0;
                   
                            Name = _name;  
                            Width = 0; 
                            Upper = 0;
                            Lower = 0;    
                            Number = 0;           
                            }
       
                 int Align(int depth, int x);
                 int Get_Y();
                 
                 int Read_file(const char * infilename, string sam_name, int sam_n);                                                                              
                 void Out_Tree_Html(const char * outfilename, int depth);
                 void Out_Tree_Distribution(const char * outfilename, int depth);       
                 int Out_Tree_SVG(const char * outfilename, int depth);
                 static void Init(string database_path, int sam_num);                    
                                  
       private:
               int Add_Node(string * branch, string id, int sam_n, float copy_number);
               
               string Name;
               int Width;
               int Upper;
               int Lower;
               int Number;
       
               map<string, TNode *> Table; //taxa to lower level
               map<string, int *> Id_table; //only on leaf nodes, for otu

               float * Count;
               
               static void Generate_Html(ostream & out);
               static void Show_Tree_Html(ostream & out, TNode * root, int depth, int i);
               static void Show_Tree_Distribution(ostream & out, TNode * root, int depth, int i, string taxa);
               static void Make_Legend_SVG(ostream & out);
               static void Show_Tree_SVG(ostream & out, TNode * root, int depth, int x);
               
               static int Sam_num;
               static float * Sam_count;
               static int N_count; //counter for number
               static int Y; //y-axis
               static hash_map <string, float, std_string_hash> Cp_number;
               static string * Sam_name;
                                                                   
       };

int TNode::Sam_num = 0; 
float * TNode::Sam_count = NULL;
int TNode::N_count = 1;
int TNode::Y = 0;
hash_map <string, float, std_string_hash> TNode::Cp_number = hash_map <string, float, std_string_hash>();
string * TNode::Sam_name = NULL;

int TNode::Align(int depth, int x){
    
     if ((x >= depth)||(Name == "Unclassified" )){
           Upper = Y;
           Lower = Y;
           Width = 1;
           Y++;
           return 1;
           } 

    int temp_align = 0; 
    //int init_y = y;
    //int final_y = y;      
    
    map <string, TNode *>::iterator map_iter = Table.begin();
    while (map_iter != Table.end()){
          
          map_iter->second->Number = N_count;
          
          N_count ++;
          
          temp_align += map_iter->second->Align(depth, x+1);
          
          Width += map_iter->second->Width;
          
          //if (map_iter != Table.begin()) y++;
                    
          map_iter ++;
    
          }
    map_iter --;
    Upper = Table.begin()->second->Get_Y();
    Lower = map_iter->second->Get_Y();
    return temp_align;
    
    }

int TNode::Get_Y(){
    
    return (Upper + Lower)/2;
    
    }

int TNode::Add_Node(string * branch, string id, int sam_n, float copy_number){

    TNode * pt = this;
    
          for (int i = 0; i< 7; i++){
              
              //cin >> buffer;
        
              if (pt->Table.count(branch[i]) == 0){
                                    pt->Table[branch[i]] = new TNode(branch[i]);
                                    pt->Table[branch[i]]->Count[sam_n] = 1.0 / copy_number;
                                    }
              else pt->Table[branch[i]]->Count[sam_n] += (1.0 / copy_number);
        
              pt = pt->Table[branch[i]];
              
              if ((branch[i] == "Unclassified")||(i == 6)) {
                                               if (pt->Id_table.count(id) == 0){
                                                  //Init the id_table
                                                  pt->Id_table[id] = new int [Sam_num];
                                                  for (int j = 0; j< Sam_num; j++)
                                                      (pt->Id_table[id])[j] = 0;
                                                  (pt->Id_table[id])[sam_n] = 1;
                                                  }
                                               else (pt->Id_table[id])[sam_n] ++;
                                               break;
                                               }
              } 
    //pt = root;
    return 0;
}

int TNode::Read_file(const char * infilename, string sam_name, int sam_n){
    
    ifstream infile(infilename, ifstream::in);
    if (!infile){
                 cerr <<"Error: Open infile error :" << infilename << endl;
                 //system("pause");
                 exit(0);
                 }
    
    int line_count = 0;
    float count = 0;
    string buffer;
    getline(infile, buffer); //Lable
    while(getline(infile, buffer)){
                          if (buffer.size()==0) continue;
                          int bio_level = 0;

                          int str_length = buffer.size();
                          string temp = "";
                          string branch[7];
                          
                          // To get Id
                          int id_begin = buffer.find('\t')+1;
                          int id_end = buffer.find('\t', id_begin);
                          string id = buffer.substr(id_begin, id_end - id_begin).c_str();
                          
                          // To get copy number
                          float copy_number = 1.0;
                          if (Cp_number.count(id) != 0)
                             copy_number = Cp_number[id];
                          
                          int str_iter = buffer.find_last_of('\t') + 1;//To locate last of '\t'
                          
                          while ((bio_level < 7) && (str_iter < str_length)){
                                
                                if ((buffer[str_iter] == ';')||(str_iter == str_length - 1)){

                                      if (temp.find("otu_")!=string::npos) temp = "Unclassified";
                                      if (temp.find("other")!=string::npos) temp = "Unclassified";
                                      if (temp.find("unclassified")!=string::npos) temp = "Unclassified";
                                      
                                      branch[bio_level] = temp;
                                                    
                                      if (temp == "Unclassified") break;
                                      temp = "";
                                      
                                      bio_level ++;
                                      str_iter ++;
                                      
                                      }
                                else 
                                     temp += buffer[str_iter];
                                
                                str_iter ++;
                                }

                          Add_Node(branch, id, sam_n, copy_number);
                          
                          count += 1.0 / copy_number;
                          line_count ++;
                          }
                          
    infile.close();
    infile.clear();
    
    Sam_count[sam_n] = count;
    Sam_name[sam_n] = sam_name;
    return line_count;
    }

void TNode::Init(string database_path, int sam_num){
     
     Load_Copy_Number(database_path, Cp_number);
     Sam_num = sam_num;
     Sam_count = new float [sam_num];
     for (int i = 0; i < sam_num; i ++)
         Sam_count[i] = 0;
     Sam_name = new string [sam_num];
     return;
     }

void TNode::Generate_Html(ostream & out){
    
    out << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" <<endl;
    out << "<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"e\" lang=\"en\">" <<endl;
    out <<"\t"<< "<head>"<< endl;
    out <<"\t\t"<< "<meta charset=\"utf-8\"/>"<<endl;
    out <<"\t\t"<< "<style>"<<endl;
    out <<"\t\t\t"<< "body {margin:0;}"<<endl;
    out <<"\t\t"<< "</style>"<<endl;
    out <<"\t\t"<< "</head>"<<endl;
    out << "\t\t"<< "<body style=\"padding:0;position:relative\">" <<endl;
    out <<"\t\t\t"<< "<div id=\"options\" style=\"position:absolute;left:0;top:0;\">" << endl;
    out <<"\t\t\t"<< "</div>" << endl;
    out <<"\t\t\t"<< "<div id=\"details\" style=\"position:absolute;top:1%;right:2%;text-align:right;\">" << endl;
    out <<"\t\t\t"<< "</div>" << endl;
    out <<"\t\t\t"<< "<canvas id=\"canvas\" width=\"100%\" height=\"100%\">" << endl;
    out <<"\t\t\t"<< "This browser does not support HTML5." << endl;
    out <<"\t\t\t"<< "</canvas>" << endl;
    out <<"\t\t\t"<< "<img id=\"hiddenImage\" src=\"hidden.png\" visibility=\"hidden\"/>" << endl;
    out <<"\t\t\t"<< "<script name=\"meta_viewer\" src=\"meta_viewer.js\"></script>" << endl;
    out <<"\t\t"<< "</body>" << endl;
    out <<"\t\t"<< "<data>" << endl;
    out <<"\t\t\t"<< "<options collapse=\"false\" key=\"true\"></options>" << endl;
	out <<"\t\t\t"<< "<magnitude attribute=\"magnitude\"></magnitude>" << endl;
    out <<"\t\t\t"<< "<attributes magnitude=\"Total\"></attributes>" << endl;
    
    return;
    }

void TNode::Show_Tree_Html(ostream & out, TNode * root, int depth, int i){
    
    if (i >= depth) return;
           
    map<string, TNode *>::iterator map_iter = root->Table.begin();
    while (map_iter != root->Table.end()){
          
          out << "<node name=\"" <<map_iter->first << "\" magnitude=\""; 
          
          for (int j = 0; j< Sam_num; j++){
              out <<map_iter->second->Count[j] ;
              if (j != Sam_num-1) out << ",";
              }
              
          out <<"\">"<< endl;
          
          if (map_iter->first != "Unclassified") 
          
             Show_Tree_Html(out, map_iter->second, depth, i+1);
          
          out <<"</node>"<<endl;
          map_iter ++;
          
          }
    return;
    
    }

void TNode:: Show_Tree_Distribution(ostream & out, TNode * root, int depth, int i, string taxa){
    
    if (i >= depth) 
                   return;
    
    taxa += root->Name;
    if (taxa.size()>0)
        taxa+='\t';
                       
    map<string, TNode *>::iterator map_iter = root->Table.begin();
    
    while (map_iter != root->Table.end()){                    
                            
          if (!map_iter->second->Id_table.empty()){
             for (map<string, int *>::iterator id_iter = map_iter->second->Id_table.begin(); id_iter!=map_iter->second->Id_table.end();id_iter++){
                                  
                                  float copy_number = 1;
                                  if (Cp_number.count(id_iter->first) != 0)
                                                                      copy_number = Cp_number[id_iter->first];
                                                                                                        
                                  out << id_iter->first<<"\t";
                                      
                                  for (int j = 0; j < Sam_num; j++){
                                      //out << (id_iter->second)[j]<<"\t"<<map_iter->second->count[j]<<"\t";
                                      out << setprecision(5) << (float)(id_iter->second)[j] / copy_number / Sam_count[j] * 100 << "\t" << map_iter->second->Count[j] / Sam_count[j] * 100 << "\t";
                                                      
                                      }
                                  out << taxa + map_iter->first << endl;
                                  }
             }
                                      
          if (map_iter->first != "Unclassified") 
          
             Show_Tree_Distribution(out, map_iter->second, depth, i + 1, taxa);
               
          map_iter ++;
          
          }
    return;
    
    }

void TNode::Out_Tree_Html(const char * outfilename, int depth){
    
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                 cerr <<"Error: Open outfile error :" << outfilename << endl;
                 //system("pause");
                 exit(0);
                 }
                     
    Generate_Html(outfile);
    
    outfile <<"\t\t\t"<<"<datasets names=\"";
    for (int i = 0; i< Sam_num; i++){
        outfile <<Sam_name[i];
        if (i != Sam_num -1) outfile << ",";
        }
    outfile <<"\"></datasets>"<<endl;    
    
    outfile << "<node name=\"Root\" magnitude=\"";
    
    for (int i = 0; i< Sam_num; i++){
        outfile << Sam_count[i];
        if (i != Sam_num-1) outfile << ",";
        }
        
    outfile<<"\">"<< endl;
    
    Show_Tree_Html(outfile, this, depth, 0);
    
    outfile<<"</node>"<<endl;
    outfile <<"\t\t" << "</data>" << endl;
    outfile <<"</html>" << endl;
    
    outfile.close();
    outfile.clear();
    
    return;
    }

void TNode::Out_Tree_Distribution(const char * outfilename, int depth){

    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                 cerr <<"Error: Open outfile error :" << outfilename << endl;
                 //system("pause");
                 exit(0);
                 }
    //Header
    //outfile << "Database Id\t";
    for (int i = 0; i< Sam_num; i++)
        outfile << "OTU\t" << Sam_name[i]<< "_OTU_Abundance(%)\t" << Sam_name[i] << "_Taxa_Abundance(%)\t";
    outfile <<"Kingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies"<<endl;
    
    string taxa = "";
    Show_Tree_Distribution(outfile, this, depth, 0, taxa);
    
    outfile.close();
    outfile.clear();
    
    return;
    }

void TNode::Make_Legend_SVG(ostream & out){
    
    int l_x = Node_Width * 7 + 100;
    int l_y = 100;
    
    int item_width = 30;
    
    out << "<rect x=\""<< l_x <<"\" y=\""<< l_y <<"\" width=\"150\" height=\""<<item_width * (Sam_num+1) <<"\" fill = \"none\" stroke=\"black\"/>" << endl; 
    
    int text_x = l_x + 20;
    int text_y = l_y + 20;
    
    out << "<text x=\""<< text_x+25 <<"\" y=\""<< text_y <<"""\" font-family=\"Verdana\" font-size=\"16\" style=\"fill:rgb(0,0,0)\">"<<"Legend"<<"</text>"<<endl;
     
    text_y += item_width;
    
    for (int i = 0; i< Sam_num; i++){
        
        float text_zoom = 1.0;
        if (Sam_name[i].size()>13)
           text_zoom = 14.0 / (float)Sam_name[i].size();
        
        out <<"<a xlink:href=\"taxonomy.html?dataset="<<i<<"\">" << endl;
        
        out <<"<text transform=\"matrix("<< text_zoom <<" 0 0 1 "<< text_x <<" "<< text_y <<")\" font-family=\"Verdana\" font-size=\"16\" style=\"fill:rgb(0,0,0)\">"<<Sam_name[i]<<"</text>"<<endl;
        
        out << "<rect x=\""<< text_x-10 <<"\" y=\""<< text_y-16 <<"\" width=\"8\" height=\""<< 16<<"\" style=\"fill:rgb("<<Default_color[i % Color_Count] << ");\"/>" << endl; 
        
        out <<"</a>"<<endl;
        
        text_y += item_width;
        }
    
    return;
    
    }

void TNode::Show_Tree_SVG(ostream & out, TNode * root, int depth, int x){
    
    if (x >= depth) return;
    
    string align_perfix = "";
          
    for (int i = 0; i< x; i++)
        align_perfix += '\t';
    
    //V-line
    
     if (!root->Table.empty()){
            int v_line_x = Node_Width * x -20 + 2 + X_Offset;
            int v_line_init_y = root->Upper * 50 + 39 + 8;
            int v_line_final_y = root->Lower * 50 + 39 + 8;
            
            out << align_perfix << "<line x1=\""<< v_line_x <<"\" y1=\""<< v_line_init_y <<"\" x2=\""<< v_line_x <<"\" y2=\""<< v_line_final_y <<"\" style=\"stroke:rgb(0,0,0);stroke-width:2.5\"/>" << endl;
            }            
           
    map<string, TNode *>::iterator map_iter = root->Table.begin();    
                      
    while (map_iter != root->Table.end()){
                              
          int text_x = Node_Width * x + X_Offset + 5 ;
          int text_y = map_iter->second->Get_Y() * 50 + 39 + 16;
          
          int figure_x = text_x + Text_Length + 5;
          int figure_y = text_y - 40;
          
          int line_x = text_x - 18 -5 ;
          int line_y = text_y -8;
          
          float text_zoom = 1.0;
          if (map_iter->first.size() > 10)
             text_zoom = 11.0 / (float)map_iter->first.size();
          
          if (map_iter->first == "Unclassified")
          //Without link out
             out << align_perfix << "<text transform=\"matrix("<< text_zoom <<" 0 0 1 "<< text_x <<" "<< text_y <<")\" font-family=\"Verdana\" font-size=\"16\" style=\"fill:rgb(0,0,0)\">"<<map_iter->first<<"</text>"<<endl;
          else {
         //With link out  
             out << align_perfix << "<a xlink:href=\"http://www.computationalbioenergy.org/web_service/parallel-meta/visualization_redirect.cgi?target="<<map_iter->first<<"\" target=\"_blank\">" <<endl;    
             out << align_perfix << "<text transform=\"matrix("<< text_zoom <<" 0 0 1 "<< text_x <<" "<< text_y <<")\" font-family=\"Verdana\" font-size=\"16\" style=\"fill:rgb(0,0,0)\">"<<map_iter->first<<"</text>"<<endl;
             out << align_perfix << "</a>" << endl;
               }
          
          //Branch Line
          out << align_perfix << "<line x1=\""<< line_x <<"\" y1=\""<< line_y <<"\" x2=\""<< line_x + 20 <<"\" y2=\""<< line_y <<"\" style=\"stroke:rgb(0,0,0);stroke-width:2.5\"/>" << endl;
          
          //Define Bar Graph
          out << align_perfix << "<defs>" << endl;
          
          out << align_perfix << "<g id=\"Node_"<<map_iter->second->Number<<"\" stroke=\"none\">" <<endl;

          //Get Max
          int max = 0;
          
          for (int i = 1; i< Sam_num; i++)
              if ((float)map_iter->second->Count[i]/(float)Sam_count[i] > (float)map_iter->second->Count[max]/(float)Sam_count[max])
                                                          max = i;
          
          for (int i = 0; i< Sam_num; i++){
              
              int bar_x = i * 8;
              
              float bar_height = (((float)map_iter->second->Count[i] * 40 /(float)Sam_count[i]) / ((float)map_iter->second->Count[max]/(float)Sam_count[max]));
    
              int bar_y = 40- (int) bar_height;
              
              if (bar_height - (int) bar_height > 0) bar_y --;
              
              
              out << align_perfix << "<rect x=\""<< bar_x <<"\" y=\""<< bar_y <<"\" width=\"8\" height=\""<< 40 - bar_y <<"\" style=\"fill:rgb("<<Default_color[i % Color_Count]<<");\"/>" << endl;
              }

          out << align_perfix << "</g>"<<endl;
          out << align_perfix << "</defs>"<<endl;
          
          //coordinate axis
          out << align_perfix << "<line x1=\""<< figure_x - 4 <<"\" y1=\""<< figure_y+40 <<"\" x2=\""<< figure_x + 8 * Sam_num + 4 <<"\" y2=\""<< figure_y+40 <<"\" style=\"stroke:rgb(0,0,0);stroke-width:1.5\"/>" << endl;

          
          for (int i = 0; i<= Sam_num; i++)
              out << align_perfix << "<line x1=\""<< figure_x + i * 8  <<"\" y1=\""<< figure_y+40 -4 <<"\" x2=\""<< figure_x + i * 8  <<"\" y2=\""<< figure_y+40 <<"\" style=\"stroke:rgb(0,0,0);stroke-width:1.5\"/>" << endl;
                        
          //Use Bar Graph
          out << align_perfix << "<a xlink:href=\"taxonomy.html?node=" <<map_iter->second->Number<<"\">"<< endl;
          out << align_perfix << "<use xlink:href=\"#Node_"<<map_iter->second->Number<<"\" x=\""<< figure_x <<"\" y=\""<< figure_y <<"\"/>" << endl;
          out << align_perfix << "</a>" << endl;
          
          //Recursion
          if (map_iter->first != "Unclassified") 
          
             Show_Tree_SVG(out, map_iter->second, depth, x+1);
          
          map_iter ++;

          
          }
    return;
    }

int TNode::Out_Tree_SVG(const char * outfilename, int depth){
    
    ofstream outfile(outfilename, ofstream::out);
    if (!outfile){
                 cerr <<"Error: Open outfile error :" << outfilename << endl;
                 system("pause");
                 exit(0);
                 }
        
    int item_count = Align(depth,0);
    
    //Header
    
    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"" << endl;
    outfile << "\"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl;
    outfile << "<svg width=\""<< Node_Width * 9 <<"\" height=\""<< 50 * (item_count+1) <<"\" xmlns=\"http://www.w3.org/2000/svg\"" << endl;
    outfile << "xmlns:xlink=\"http://www.w3.org/1999/xlink\">" << endl;
    outfile << "<desc>Reusing items</desc>" << endl;
    
    Make_Legend_SVG(outfile);

    Show_Tree_SVG(outfile, this, depth, 0);
    
    outfile << "</svg>" << endl;
    
    outfile.close();
    outfile.clear();
    
    return item_count;
    }

#endif


