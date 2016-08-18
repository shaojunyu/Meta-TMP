#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <cmath>
using namespace std;

int main(int argc, char const *argv[]) {
  ifstream inFile(argv[1]);
  ofstream outHtmlFile("index.html");
  map<string, float> name_map;
  string line;
  string A;
  string B;
  float corr;

  if (inFile.is_open())
  {
    outHtmlFile<<"<!DOCTYPE html>"<<endl;
    outHtmlFile<<"<html><head><title>Basic sigma.js example</title>"<<endl;
    outHtmlFile<<"<style type=\"text/css\">body {margin: 0;}#container {position: absolute;width: 100%; height: 100%;}</style>"<<endl;
    outHtmlFile<<"</head>"<<endl;
    outHtmlFile<<"<body>"<<endl;
    outHtmlFile<<"<div id=\"container\"></div>"<<endl;
    outHtmlFile<<"<script src=\"./dist/sigma.min.js\"></script>"<<endl;
    outHtmlFile<<"<script>"<<endl;
    outHtmlFile<<"g = {nodes: [],edges: []};"<<endl;

    int id = 1;
    char color[20];
    // sprintf(color, "rgb(");
    while(inFile>>A>>B>>corr){
      if (corr >= 0)
      {
        sprintf(color, "rgb(%d,0,0)", (int)(corr*255));
      }else{
        sprintf(color, "rgb(0,%d,0)", (int)(-corr*255));
      }
      cout<<color<<endl;
      if (name_map.find(A) == name_map.end())
      {
        name_map[A] = corr;
        outHtmlFile<<"g.nodes.push({id:'"<<A<<"',label:'"<<A<<"',x: Math.random(),y: Math.random(),size: Math.random(),color: '#666'});"<<endl;
      }
      if (name_map.find(B) == name_map.end())
      {
        name_map[B] = corr;
        outHtmlFile<<"g.nodes.push({id:'"<<B<<"',label:'"<<B<<"',x: Math.random(),y: Math.random(),size: Math.random(),color: '#666'});"<<endl;
      }
      outHtmlFile<<"g.edges.push({id:"<<id<<",label:'"<<corr<<"',source:'"<<A<<"',target:'"<<B<<"',size:"<<fabs(corr)<<",color: '"<<color<<"',type:'line'});"<<endl;
      id++;
    }

    outHtmlFile<<"var s = new sigma({graph:g,renderer:{container: document.getElementById('container'),type:'canvas'},"<<endl;
    outHtmlFile<<"settings: {edgeLabelSize: 'proportional',edgeLabelThreshold: 3,minEdgeSize: 0.1}"<<endl;
    outHtmlFile<<"});"<<endl;
    outHtmlFile<<"s.refresh();"<<endl;
    outHtmlFile<<"</script>"<<endl;
    outHtmlFile<<"</body>"<<endl;
    outHtmlFile<<"</html>"<<endl;
  }else{
    cout<< argv[1];
  }
  return 0;
}