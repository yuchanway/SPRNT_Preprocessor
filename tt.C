/* 
  read the csv file and build the graph
*/

#include <stdio.h>
#include <vector>
#include "smap.h"
#include "tgraph.h"
#include "nc_wrapper.h"

#define COMMENT_TOKEN '#'
#define TOKEN_DELIM ", \t"
//#define steady_flag 0

// name of the netlist file
const char* outfile="Junction_Test_1.spt";

void printcontrol(FILE *fp) {
  fprintf(fp, "\n######\n");      // options
  fprintf(fp, "\n######\n");
  fprintf(fp, "def options metric=1 verbose=1 end\n"); // properly no need to change 
  fprintf(fp, "def options prtq=1 prta=1 prtdepth=1 end\n");
  fprintf(fp, "def options stoptime=1488.0 stoptimeunit=hour end\n");     // 24x30 hours = 1 month
  fprintf(fp, "def options prtinterval=2.0 prtintervalunit=hour end\n"); // printing interval
  //fprintf(fp, "def options timestep=5 timestepunit=minute end\n"); // 240 sec time step
  fprintf(fp, "def options ssfile=%s.ssf end\n", outfile);
  fprintf(fp, "def options epoch=2010-01-01T00:00:00Z end\n"); // time is in UTC
  fprintf(fp, "\n######\n");
  fprintf(fp, "\n###### end #########\n\n");
}

void printheader(FILE *fp, const char *conn_fname, const char *flow_fname, const char *root) {
  fprintf(fp,"########################################################################\n");
  fprintf(fp,"###### generated using connectivity file: \"%s\"\n", conn_fname);
  fprintf(fp,"######                 flow file: \"%s\"\n", flow_fname);
  fprintf(fp,"######                 root node: \"%s\"\n", root);
  fprintf(fp,"########################################################################\n\n");
}

int main(int argc, char* argv[]) {
  const int MAX_LINE_LENGTH = 999999;
  char   tmpbuf[MAX_LINE_LENGTH];
  char   brch_name[MAX_LINE_LENGTH], frm_name[MAX_LINE_LENGTH], to_name[MAX_LINE_LENGTH];
  char  *pt;
  int   fnum, eid, up_id, dn_id, lnum, tp;
  double len, maf, zr, slope, lon, lat, Manning_n;
  int steady_flag = 0;
  // ***Revised by Cheng-Wei *** New variables
  int shape;
  double wid, SWS;
  // ***Revised by Cheng-Wei *** New vector for saving x and y coordinates
  vector<double> xData, yData;
  vector<double> AAData, PPData, YYData, WWData; //P (perimeter), y(depth) and w(free surface width) all depend on A;
  FILE *F;
  tgraph  TG;
  NC      runoff;
  NC      QSOURCE;

  if (argc<4) {
    printf("Usage: %s <NHD+ connectivity csv file> <netcdf lateral file> <netcdf qsource file> [root name] -s(optional)\n", argv[0]);
    return (1);
  }
  
  if ((argc== 6) && (strcmp(argv[5],"-s") == 0)){
    steady_flag= 1;
  }
  else if ((argc == 6) && (strcmp(argv[5],"-s") != 0)){
    printf("[EE] Please use valid steady state flag -s \n");
    printf("Usage: %s <NHD+ connectivity csv file> <netcdf file> <netcdf qsource file> [root name] -s(optional)\n", argv[0]);
    return (1);
  }
  //fprintf(stdout,"steady_flag is %i \n", steady_flag);

  F = fopen(argv[1],"r"); //open connectivity file in csv 
  if (!F) {
    fprintf(stdout,"bummer: unable to open file \"%s\"\n", argv[1]);
    return (-1);
  }
  
  lnum = 0;
  while ( (fgets(tmpbuf, MAX_LINE_LENGTH, F)) != NULL ) {
    if ( tmpbuf[0] == COMMENT_TOKEN ) continue;
    pt = strtok(tmpbuf, TOKEN_DELIM);
    fnum = 0;
    //***Revised by Cheng-Wei*** New variable
    shape = -1;
    // ***Revised by Cheng-Wei*** Declare tp as interger
    tp = 1.0;
    len = -1.0; 
    maf = -1.0;
    slope = -1.0; 
    lon = -1.0;
    lat = -1.0;
    Manning_n = -1.0;     
      
    // ***Revised by Cheng-Wei*** New variables SWS = side wall slope; wid = bottom width
    wid = 1.0;
    SWS = 1.0;

    xData.clear();
    yData.clear();
    AAData.clear();
    PPData.clear();
    YYData.clear();
    WWData.clear();
     
    
   
    while ( pt != NULL ) {
          if (fnum <= 11) {
              switch (fnum) {
                     case 0: 
                          strncpy(brch_name, pt, MAX_LINE_LENGTH-1);
                          //	eid  = TG.MakeEdge(pt);
                          break;
                     case 1:
                          strncpy(frm_name, pt, MAX_LINE_LENGTH-1);
                          //	up_id = TG.MakeVertex(pt);
                          break;
                     case 2:
                          strncpy(to_name, pt, MAX_LINE_LENGTH-1);
                          //	dn_id = TG.MakeVertex(pt);
                          break;
                     case 3:
                          tp = atoi(pt);
                          break;
                     case 4:
                          len = atof(pt);
                          break;
                     case 5:
                          maf = atof(pt);
                          break;
                     case 6:
                          slope = atof(pt);
                          break;
                     case 7:
                          zr = atof(pt);
                          break;
                     case 8:
                          Manning_n = atof(pt);
                          break;
                     case 9:
                          lat = atof(pt);
                          break;
                     // ***Revised by Cheng-Wei	*** 9th column data is shape
                     case 10:
                          lon = atof(pt);
                          break;
                     case 11:
                     	  shape = atoi(pt);
                     	  break;
                     default: break;
              }
          }
          // ***Revised by Cheng-Wei *** For 11th or behind columns,
          // If shape = 2, read the following AA(Wetarea), PP(perimeter), YY(depth), WW(free surface width) --> instrinsic XS
          // If shape = 1, read the following x and y coordinates.     --> XY data XS
          // If shape = 0, read only 2 colmuns and save as wid and SWS --> Trapezoidal XS.
          else if (shape == 0){
               if (fnum == 12) {wid = atof(pt);}
               if (fnum == 13) {SWS = atof(pt);}
          }
          else if (shape == 1){
               if (fnum % 2 == 0) {xData.push_back(atof(pt));}
               else {yData.push_back(atof(pt));}
          }
          else if (shape == 2){
              if (fnum % 4 ==0) {AAData.push_back(atof(pt));}
              if (fnum % 4 ==1) {PPData.push_back(atof(pt));}
              if (fnum % 4 ==2) {YYData.push_back(atof(pt));}
              if (fnum % 4 ==3) {WWData.push_back(atof(pt));}
          }
               
          pt = strtok(NULL, TOKEN_DELIM);
          fnum++;
    }
              
    
    if (tp >= 2 ) {
      //      fprintf(stderr,"minor branch %s ignored\n", brch_name);
    } else {
      eid = TG.MakeEdge(brch_name);
      up_id = TG.MakeVertex(frm_name);
      dn_id = TG.MakeVertex(to_name);
      TG.AssignLength(eid, len);
      TG.AssignMAF(eid, maf);
      TG.AssignZr(eid, zr);
      TG.AssignSlope(eid, slope);
      TG.AssignLon(eid, lon);
      TG.AssignLat(eid, lat);
      TG.AssignManning(eid, Manning_n);
      
      // ***Revised by Cheng-Wei *** New assign fuctions
      TG.AssignShape(eid, shape);
      // ***Revised by Cheng-Wei *** If shape = 0,
      // assign only wid and SWS to TG(class tgraph)
      if (shape == 0){
            TG.AssignWid(eid, wid);
            TG.AssignSWS(eid, SWS);
      }
      // ***Revised by Cheng-Wei *** If shape = 1,
      // assign the XY vector to TG(class tgraph)
      if (shape == 1){
          for (int i = 0; i < xData.size(); i++){
              TG.AssignX(eid, xData[i]);
              TG.AssignY(eid, yData[i]);
          }
      }
        if (shape ==2){
            for (int i= 0; i < AAData.size(); i++){
                TG.AssignAA(eid, AAData[i]);
                TG.AssignPP(eid, PPData[i]);
                TG.AssignYY(eid, YYData[i]);
                TG.AssignWW(eid, WWData[i]);
            }
        }
          
      
      TG.ConnectEdge(eid, up_id, dn_id);
    }
    lnum++;
  }
  if (F) fclose(F);

  TG.QuickCheck();
  
  TG.DFS();  // DFS to look at connectivities

  printf("Reading runoff data from %s\n", argv[2]);

  int rc = runoff.Init(argv[2]);
  if ( rc < 0 ) {//return (-1);
    printf("control point of lateral flow file open.");
    return (-1);
  }

  int rc_2 = QSOURCE.Init(argv[3]);
  if ( rc_2 < 0 ) {//return (-1);
    printf("control point of qsource flow file open.");
    return (-1);
    }

#if 0
  printf("Dump:\n");
  printf("%d %d\n", runoff.NumComid(), runoff.NumTime());
  runoff.Dump(stdout, 3);
#endif

  if ((argc>=5) && (steady_flag==0)) {
    FILE *fp=fopen(outfile,"w");
    printf("[II] Netlist will be saved as %s\n", outfile);
    if (fp) printheader(fp, argv[1], argv[2],argv[3]);
    TG.DFSupstream( argv[4], &runoff, &QSOURCE, fp, steady_flag);
    printcontrol(fp);
    fclose(fp);
  } 
  else if((argc>=5) && (steady_flag==1)) {
    FILE *fp=fopen(outfile,"w");
    printf("[II] Steady State Netlist will be saved as %s\n", outfile);
    if (fp) printheader(fp, argv[1], argv[2],argv[3]);
    TG.DFSupstream( argv[4], &runoff, &QSOURCE, fp, steady_flag);
    printcontrol(fp);
    fclose(fp);
  }

  return 0;
}

// Local Variables:
// mode: c++
// End:

