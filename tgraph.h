/* 
  a simple graph package
*/

#ifndef _T_GRAPH_H
#define _T_GRAPH_H

#include "smap.h"
#include "util.h"
#include "nc_wrapper.h"
#include <vector>

typedef int edge_id;
typedef int vertex_id;

/*
  Each edge should have one upstream vertex and one downstream vertex
  
  Each vertex can have multiple fan-in edges from upstream, and multiple fan-out
  edges to downstream

 */
class edge;

class vertex {
protected:
  int                _id;
  int                _flag;
  int                _num_fanin;
  int                _num_fanout;

  FlexVec<edge*,3>   _fanins;  // assuming each node has up to 3 fan-in
  FlexVec<edge*,3>   _fanouts;

  char               _name[16];  

public:
  vertex():_id(-1),_flag(0),_num_fanin(0),_num_fanout(0) { _name[0]='\0';}
  ~vertex() {}

  /* methods */
  int& Id() { return _id; }
  int& Flag() { return _flag;}
  int GetNumFanin() const { return _num_fanin; }
  int GetNumFanout() const { return _num_fanout; }
  edge* GetIthFanin(int j) { return _fanins[j]; }
  edge* GetIthFanout(int j) { return _fanouts[j]; }
  void SwapFanin(int j, int k) {
    assert( j<_num_fanin && k<_num_fanin);
    edge* tmp = _fanins[j];
    _fanins[j] = _fanins[k];
    _fanins[k] = tmp;
  }
  char* Name() { return _name; }
  void SetName(const char *n) { strncpy(_name, n,15); _name[15]='\0';}
  
  int AddEdgeOut(edge *e) {
    _fanouts.grow(_num_fanout+1);
    _fanouts[_num_fanout] = e;
    return (_num_fanout++);
  }

  int AddEdgeIn(edge *e) {
    _fanins.grow(_num_fanin+1);
    _fanins[_num_fanin] = e;
    return (_num_fanin++);
  }

    };

class edge {
protected:
  int       _id;
  int       _flag;
  double    _length;
  double    _maf;
  double    _zr;
  double    _slope;
  double    _long;
  double    _lat;
  
  // ***Revised by Cheng-Wei *** New variables
  double    _manning;
  int       _shape;
  double    _wid;
  double    _SWS;
  // ***Revised by Cheng-Wei *** New vectors for saving vectors for x, y, AA, PP,......
  std::vector<double> _x;
  std::vector<double> _y;
  std::vector<double> _AA;
  std::vector<double> _PP;
  std::vector<double> _WW;
  std::vector<double> _YY;


  vertex*   _up;
  vertex*   _dn;

  char      _name[16];

public:
  edge():_id(-1),_flag(0),_length(0.0),_up(NULL),_dn(NULL) { _name[0]='\0';}
  ~edge() {}
  
  /* methods */
  
  // REFERENCE MAKES PROTECTED VARIABLES UNPROTECTED...
  int& Id() { return _id; }
  int& Flag() { return _flag; }
  char* Name() { return _name; }
  void SetName(const char *n) { strncpy(_name, n,15); _name[15]='\0';}
  double& Length() { return _length; }
  double& MAF() { return _maf; }
  double& Zr() { return _zr; }
  double& Slope() { return _slope; }
  double& Lon() { return _long; }
  double& Lat() { return _lat; }
  double& Manning() { return _manning; }
  
  // ***Revised by Cheng-Wei *** New save functions (avoid using reference)
  void Save_x(double data) {_x.push_back(data); }
  void Save_y(double data) {_y.push_back(data); }
    //*** Revised by Cheng-Wei *** intrinsic cross section
  void Save_AA(double data) {_AA.push_back(data); }
  void Save_PP(double data) {_PP.push_back(data); }
  void Save_WW(double data) {_WW.push_back(data); }
  void Save_YY(double data) {_YY.push_back(data); }
  void Save_shape(int data) {_shape = data; }
  void Save_wid(double data) {_wid = data; }
  void Save_SWS(double data) {_SWS = data; }
  
  // ***Revised by  Cheng-Wei *** New callout functions, these callout will be used in printing .spt
  double Read_x(int num) {return _x[num]; }
  double Read_y(int num) {return _y[num]; }
  double Read_AA(int num) {return _AA[num];}
  double Read_PP(int num) {return _PP[num];}
  double Read_YY(int num) {return _YY[num];}
  double Read_WW(int num) {return _WW[num];}
  int Length_xy() {return _x.size(); }
  int Length_AA() {return _AA.size(); }
  int Read_Shape() {return _shape; }
  double Read_Wid() {return _wid; }
  double Read_SWS() {return _SWS; }
   
  
  vertex* Up() { return _up; }
  vertex* Dn() { return _dn; }
  void SetUpVert(vertex *up) { _up = up; }
  void SetDnVert(vertex *dn) { _dn = dn; }

};

// tgraph holds all vertices and edges
class tgraph {
protected:

  SMap                _v_map;
  SMap                _e_map;
  
  int                 _num_vertices;
  int                 _num_edges;
  int                 steady_f;
  FlexVec<vertex*>    _Vertices;
  FlexVec<edge*>      _Edges;

public:
  tgraph() { _num_vertices=0; _num_edges=0;}
  ~tgraph();

  vertex_id MakeVertex(const char* n);
  edge_id   MakeEdge(const char *n);
  void AssignLength(edge_id eid, double l) { _Edges[eid]->Length()=l; }
  void AssignMAF(edge_id eid, double m) { _Edges[eid]->MAF()=m; }

  void AssignZr(edge_id eid, double l) { _Edges[eid]->Zr()=l; }

  void AssignSlope(edge_id eid, double l) { _Edges[eid]->Slope()=l; }
  void AssignLon(edge_id eid, double l) { _Edges[eid]->Lon()=l; }
  void AssignLat(edge_id eid, double l) { _Edges[eid]->Lat()=l; }
  void AssignManning (edge_id eid, double l) { _Edges[eid] -> Manning()=l; }

  // ***Revised by Cheng-Wei ***  New assign functions
  void AssignX (edge_id eid, double l) {_Edges[eid]->Save_x(l); }
  void AssignY (edge_id eid, double l) {_Edges[eid]->Save_y(l); }
  void AssignAA (edge_id eid, double l) {_Edges[eid]->Save_AA(l);}
  void AssignPP (edge_id eid, double l) {_Edges[eid]->Save_PP(l);}
  void AssignYY (edge_id eid, double l) {_Edges[eid]->Save_YY(l);}
  void AssignWW (edge_id eid, double l) {_Edges[eid]->Save_WW(l);}
  void AssignShape (edge_id eid, int l) {_Edges[eid]->Save_shape(l); }
  void AssignWid (edge_id eid, double l) {_Edges[eid]->Save_wid(l); }
  void AssignSWS (edge_id eid, double l) {_Edges[eid]->Save_SWS(l); }

  void ConnectEdge(edge_id eid, vertex_id up_id, vertex_id dn_id);

  void QuickCheck();  // some simple printing
  int DFS();  // DFS traversal, returns the number of connected regions
  int DFSupstream(const char *root_name, NC *runoff, NC *QSOURCE, FILE* fp=NULL, int steady_f=0); // DSF traversal to create a netlist

};


#endif

// Local Variables:
// mode: c++
// End:
