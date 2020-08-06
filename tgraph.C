/* 
   methods of tgraph 
*/

#include <stdio.h>
#include <math.h>

#include "tgraph.h"

#define PRINT_CRAP 0


// minimal length to split
/* This MAGIC number is the distance interval (delta x) beem used in the numerical method. The value of this MAGIC number may be changed for users' preference. */
#define MAGIC 200 // unit = meter


/* qsource is the upstream boundary condition Q. If user do not assign any hydrograph at the specific location for upstream BC, this preprocessor will automatically assign all the upstream BC = qsource. */
const double upstream_BC = 0.1;



/* Variable "downstream_BC" is used to define the downstream boundray condition. In here, the BC is depth and the unit is in METERs. SPRNT program also provide another boundrau condition option, which is wetted area. However, defining downstream BC by using wetted area requires extra computation since the wetted area is a function of river geometry and bathemetry. Using depth may be a more straight forward option. */
const double downstream_BC = 3.0; // unit = meter



// we either fetch the existing one
// or create a new one
vertex_id tgraph::MakeVertex(const char *n) {
  vertex *v;
  int idx = _v_map.Check(n);
  if ( idx == -1 ) { // new
    idx = _v_map.Create(n);
    v=new vertex;
    v->SetName(n);
    v->Id() = idx;
    _Vertices.grow( _num_vertices + 1 );
    _Vertices[ _num_vertices ] = v;
    return ( _num_vertices++);
  } else { // existing
    return (idx);
  }
}

// each edge can only appear once
edge_id tgraph::MakeEdge(const char *n) {
  int idx = _e_map.Check(n);
  if ( idx != -1 ) {
    fprintf(stdout, "bummer: duplicate edge w/ name %s\n", n);
    exit (-1);
  }

  idx = _e_map.Create(n);
  edge* ee=new edge;
  ee->SetName(n);
  ee->Id() = idx;
  _Edges.grow( _num_edges + 1 );
  _Edges[ _num_edges ] = ee;
  return ( _num_edges++);
}

// connect the edge to two vertices
void tgraph::ConnectEdge(edge_id eid, vertex_id up_id, vertex_id dn_id) {
  vertex *up, *dn;
  edge   *ee;

  up = _Vertices[up_id];
  dn = _Vertices[dn_id];
  ee = _Edges[eid];
  
  ee->SetUpVert( up );
  ee->SetDnVert( dn );

  up->AddEdgeOut( ee );
  dn->AddEdgeIn( ee );

}

tgraph::~tgraph() {
  for (int jj=0; jj<_num_vertices; jj++) {
    if ( _Vertices[jj] ) delete _Vertices[jj];
  }
  for (int jj=0; jj<_num_edges; jj++) {
    if ( _Edges[jj] ) delete _Edges[jj];
  }
}

// simple printing
void tgraph::QuickCheck() {
  printf("Total %d nodes.\n", _num_vertices);
#if PRINT_CRAP
  for (int jj=0; jj<_num_vertices; jj++) {
    if ( _Vertices[jj]->GetNumFanin() == 0 ) {
      printf("Found leaf nodes: %d - %s\n", jj, _Vertices[jj]->Name());
    }
  }
  for (int jj=0; jj<_num_vertices; jj++) {
    if ( _Vertices[jj]->GetNumFanout() == 0 ) {
      printf("Found root nodes: %d - %s\n", jj, _Vertices[jj]->Name());
    }
  }
  for (int jj=0; jj<_num_vertices; jj++) {
    if (_Vertices[jj]->GetNumFanout() > 1 ) {
      printf("Found multiple down branches: %d - %s\n",
	     jj, _Vertices[jj]->Name());
    }
  }
  for (int jj=0; jj<_num_vertices; jj++) {
    if (_Vertices[jj]->GetNumFanin() > 1 ) {
      printf("Found junction: %d - %s\n",
	     jj, _Vertices[jj]->Name());
    }
  }
  for (int jj=0; jj<_num_vertices; jj++) {
    if (_Vertices[jj]->GetNumFanin() > 1 && _Vertices[jj]->GetNumFanout() > 1 ) {
      printf("Found odd ball node: %d - %s\n",
	     jj, _Vertices[jj]->Name());
    }
  }
#endif
  int cnt=0;
  for (int jj=0; jj<_num_edges; jj++) {
    if (_Edges[jj]->Length() < MAGIC ) {
      _Edges[jj]->Length() = MAGIC;
      cnt++;
    }
  }
  printf("Fixed %d edges with length less than %.1f\n", cnt, (double)MAGIC);
}

// DFS traveral
// returns number of partitions
int tgraph::DFS() {
  MyStack<int>     S;
  int    n_roots, vid, rc, region_tag;
  edge    *e;

  // clear the flags
  for (int jj=0; jj<_num_vertices; jj++) {
    _Vertices[jj]->Flag() = -1;
  }

  // assign domain tags, starting from 0
  n_roots = 0;
  for (int jj=0; jj<_num_vertices; jj++) {
    if ( _Vertices[jj]->GetNumFanout() == 0 ) { 
      _Vertices[jj]->Flag() = n_roots++;
      S.Push( jj );  // they are also starting point
    }
  }
  
  // DFS to get them all
  while (1) {
    rc = S.Pop( vid );
    if ( rc == -1 ) break; 

    region_tag = _Vertices[ vid ]->Flag();

    // if there are more than 2 fanin, we sort them by MAF so that we only use the 2
    // largest ones
    if ( _Vertices[ vid ]->GetNumFanin() > 2 ) {
      for (int jj=1; jj<_Vertices[ vid ]->GetNumFanin(); jj++) {
	if ( _Vertices[ vid ]->GetIthFanin(0)->MAF() < _Vertices [vid]->GetIthFanin(jj)->MAF()) {
	  _Vertices[ vid ]->SwapFanin(0,jj);
	}
      }
      for (int jj=2; jj<_Vertices[ vid ]->GetNumFanin(); jj++) {
	if ( _Vertices[ vid ]->GetIthFanin(1)->MAF() < _Vertices [vid]->GetIthFanin(jj)->MAF()) {
	  _Vertices[ vid ]->SwapFanin(1,jj);
	}
      }
    }

    for (int jj=0; jj<_Vertices[ vid ]->GetNumFanin(); jj++) {
      e = _Vertices[ vid ]->GetIthFanin(jj);
      e->Up()->Flag() = region_tag;
      
      // we only push into the stack if there is any fanin
      if ( e->Up()->GetNumFanin() > 0 ) {
	S.Push( e->Up()->Id());
      }
    }
    
  }
  
  // do checking and printing
  for (int jj=0; jj<_num_vertices; jj++) {
    if ( _Vertices[jj]->Flag() == -1 ) {
      printf("orphan node %d\n", jj);
    }
  }

  // print the connected regions
#if PRINT_CRAP
  for (int kk=0; kk<n_roots; kk++) {
    printf("Region %d:\n", kk);
    for (int jj=0; jj<_num_vertices; jj++) {
      if ( _Vertices[jj]->Flag() == kk) {
	printf("%d: %d = %s %1s\n", 
	       kk, jj, _Vertices[jj]->Name(),
	       _Vertices[jj]->GetNumFanin()==0?
	       "L":(_Vertices[jj]->GetNumFanout()==0?
		    "R":" "));
      }
    }
    printf("\n\n");
  }
#endif
  for (int kk=0; kk<n_roots; kk++) {
    int cnt=0;
    char tmp[512];
    for (int jj=0; jj<_num_vertices; jj++) {
      if ( _Vertices[jj]->Flag() == kk ) {
	cnt++;
	if ( _Vertices[jj]->GetNumFanout()==0) strncpy(tmp, _Vertices[jj]->Name(), 511);
      }
    }
    printf("Region %d, root = %s, num vertices = %d\n", kk, tmp, cnt);
  }
  return n_roots;
}

// DFS from a starting point
// also creates the netlist
int tgraph::DFSupstream(const char* root, NC *runoff, NC *QSOURCE, FILE *fp, int steady_f) {
  MyStack<int>    S;
  int             root_int;
  int             npts = runoff->NumTime(); //this is total time step in the runoff file
  double         *buf = (double*)malloc( (npts+10)*sizeof(double));
  int             steady_flag = steady_f; // if steady mode on, print lateral flow at t=0 and t= 1

  root_int = _v_map.Check( root );
  if ( root_int == -1 ) {
    fprintf(stderr,"bummer: unable to locate node \"%s\"\n",
	    root);
    return (-1);
  }

  if ( _Vertices[ root_int ]->GetNumFanin() == 0 ) {
    fprintf(stderr,"Node \"%s\" (%d) has no upstream nodes\n",root, root_int);
    return (-1);
  }
  
  if ( _Vertices[ root_int ]->GetNumFanin() > 1 ) {
    fprintf(stderr,"Node \"%s\" (%d) has more than one upstream nodes\n",root, root_int);
    return (-1);
  }

  printf("Traversing upstream from root %s\n", root);

  S.Push( root_int );

  MyStack<edge*>  e_stack;

  int rc, vid; //cntr;
  int node_w_coord;
  edge  *e=NULL, *e_dn=NULL;
  int num;
  double len_step, wid;
  double total_lat, total_q;

  total_lat = 0.0;
  total_q = 0.0;
  //cntr = 0;
  
  while (1) {
    rc = S.Pop(vid);
    if ( rc == -1 ) break;
/*
#if PRINT_CRAP  
    printf("%d : %s\n", cntr++, _Vertices[ vid ]->Name());
#endif
*/
    for (int jj=0; jj<_Vertices[vid]->GetNumFanin(); jj++) {
      if ( jj>=2 ) break;

      e = _Vertices[vid]->GetIthFanin(jj);
      S.Push( e->Up()->Id() );
#if PRINT_CRAP
      printf("   %s %7.1f %7.1f %s %s\n", e->Name(), e->Length(),e->MAF(), e->Up()->Name(),
	     e->Length()<100?"S":" ");
#endif

      if ( fp ) {
	// we print the nodes first
	if ( e->Length() < MAGIC ) {
	  len_step = (MAGIC)/2.0;
	  num = 1;  // number of nodes to insert!
	} else {
	  num = (int)(e->Length())/(MAGIC);
	  len_step = e->Length()/(num+1);
	}

	node_w_coord = num/2;

	wid = 0.7*sqrt(e->MAF()); // geomorphology
	if (wid<0.2) wid=0.2;     // with minimal width
	
	// now we just print
	fprintf(fp,"def node id=%s_%s sR=%.1e n=%.3f zR=%.2e hR=0\n", e->Name(), e->Up()->Name(), e->Slope(), e->Manning(), e->Zr());
	// ***Revised by Cheng-Wei *** New print out
	if (e->Read_Shape() == 0){
        fprintf(fp,"   def trapezoidal bottomwidth=%.2f slope=%0.3f end\n", e->Read_Wid(), e->Read_SWS());
    }
	else if (e->Read_Shape() == 1){
        fprintf(fp,"   def xy \n");
    	for (int i = 0; i < e->Length_xy(); i++){
            fprintf(fp,"     x=%f y=%f\n", e->Read_x(i), e->Read_y(i));
        }
        fprintf(fp,"   end\n");
    }
    else if(e->Read_Shape()==2){
        fprintf(fp,"   def intrinsic\n");
        for (int i=0; i< e->Length_AA(); i++){
            fprintf(fp,"     A=%.3f P=%.3f Y=%.3f W=%.3f\n", e->Read_AA(i), e->Read_PP(i), e->Read_YY(i), e->Read_WW(i));
        }
        fprintf(fp,"   end\n");
    }
	fprintf(fp,"   end\n");

          
  

	for (int kk=0; kk<num; kk++) {
	  fprintf(fp,"def node id=%s_%d sR=%.1e n=%.3f zR=%.2e hR=0", e->Name(), kk,e->Slope(), e->Manning(), e->Zr());
	  if ( kk == node_w_coord ) {
	    fprintf(fp," xcoord=%.5f ycoord=%.5f\n", e->Lon(), e->Lat());
	  } 
      else {
	    fprintf(fp,"\n");
	  }
    // ***Revised by Cheng-Wei*** New print out
    if (e->Read_Shape() == 0){
        fprintf(fp,"   def trapezoidal bottomwidth=%.2f slope=%0.3f end\n", e->Read_Wid(), e->Read_SWS());
    }
    else if (e->Read_Shape() == 1){
        fprintf(fp,"   def xy \n");
    	for (int i = 0; i < e->Length_xy(); i++){
            fprintf(fp,"     x=%f y=%f\n", e->Read_x(i), e->Read_y(i));
        }
        fprintf(fp,"   end\n");
    }
    else if(e->Read_Shape()==2){
        fprintf(fp,"   def intrinsic\n");
        for (int i=0; i< e->Length_AA(); i++){
            fprintf(fp,"     A=%.3f P=%.3f Y=%.3f W=%.3f\n", e->Read_AA(i), e->Read_PP(i), e->Read_YY(i), e->Read_WW(i));
        }
        fprintf(fp,"   end\n");

    }



	  fprintf(fp,"end\n");
	}
	
	fprintf(fp,"def node id=%s_%s sR=%.1e n=%.3f zR=%.2e hR=0\n", e->Name(), e->Dn()->Name(), e->Slope(), e->Manning(), e->Zr());
	
    // ***Revised by Cheng-Wei *** New print out
	if (e->Read_Shape() == 0){
        fprintf(fp,"   def trapezoidal bottomwidth=%.2f slope=%0.3f end\n", e->Read_Wid(), e->Read_SWS());
    }
	else if (e->Read_Shape() == 1){
        fprintf(fp,"   def xy \n");
    	for (int i = 0; i < e->Length_xy(); i++){
            fprintf(fp,"     x=%f y=%f\n", e->Read_x(i), e->Read_y(i));
        }
        fprintf(fp,"   end\n");
    }
    else if(e->Read_Shape()==2){
        fprintf(fp,"   def intrinsic\n");
        for (int i=0; i< e->Length_AA(); i++){
            fprintf(fp,"     A=%.3f P=%.3f Y=%.3f W=%.3f\n", e->Read_AA(i), e->Read_PP(i), e->Read_YY(i), e->Read_WW(i));
        }
        fprintf(fp,"   end\n");
    }




    fprintf(fp,"   end\n");
	//fprintf(fp,"end\n");
	
	// print the segments
	int idx=0;
	fprintf(fp,"def segment up=%s_%s down=%s_%d length=%.1f end\n", 
	       e->Name(),e->Up()->Name(), e->Name(),idx, len_step);
	
	for (; idx<num-1;idx++) { // idx defined earlier
	  fprintf(fp,"def segment up=%s_%d down=%s_%d length=%.1f end\n", 
		 e->Name(),idx, e->Name(), idx+1,len_step);
	}
	fprintf(fp,"def segment up=%s_%d down=%s_%s length=%.1f end\n", 
	       e->Name(),idx, e->Name(), e->Dn()->Name(), len_step);


	// lateral sources
	int tmp=atoi(e->Name());
	rc = runoff->Query(tmp, buf);
	if ( rc && e->Up()->GetNumFanin() > 0 ) { // non-leaf node
	  fprintf(fp,"##### lateral flow %d ", tmp);
	  for (int ii=0; ii< (runoff->NumTime()>=4?4:runoff->NumTime()) ;ii++) {
	    fprintf(fp, "%.3e ", buf[ii]);
	  }
	  fprintf(fp,"\n");

	  total_lat += buf[0];
    double temp_lat = buf[0];
    double lat_tol = 0.00;
    int lateral_repeat_checker=0;
if (steady_flag == 0){
	  for (idx=0; idx<num; idx++) {
	    fprintf(fp, "def lateralsource\n");
	    fprintf(fp, "  location=%s_%d\n", e->Name(),idx);
	    fprintf(fp, "  def timeseries\n");
	    fprintf(fp, "    timeunit=hour\n");
 	    fprintf(fp, "    t=0 v=0.0\n"); //, buf[0]/(num));
	    for (int ti=0; ti<runoff->NumTime(); ti++){
        if (ti >=1 && ti < runoff->NumTime()-1){
          if (abs(buf[ti] - temp_lat) >= lat_tol){
            if (lateral_repeat_checker ==0){
              fprintf(fp,"    t=%d v=%.3e\n", ti, buf[ti-1]/(num));
              fprintf(fp,"    t=%d v=%.3e\n", ti+1, buf[ti]/(num));
              lateral_repeat_checker=1;}
            else if (lateral_repeat_checker ==1){fprintf(fp,"    t=%d v=%.3e\n", ti+1, buf[ti]/(num));}
            temp_lat = buf[ti];
            }
          else{lateral_repeat_checker = 0;}
        }
        else {
          fprintf(fp,"    t=%d v=%.3e\n", ti+1, buf[ti]/(num));
          lateral_repeat_checker =1;}
	    }
	    fprintf(fp,"  end\n");
	    fprintf(fp,"end\n");
	  }
    }

else if (steady_flag == 1){
    for (idx=0; idx<num; idx++){
      fprintf(fp, "def lateralsource\n");
      fprintf(fp, "   location=%s_%d\n", e->Name(), idx);
      fprintf(fp, "   def timeseries\n");
      fprintf(fp, "       t=0 v=%.3e\n", buf[0]/(num));
      fprintf(fp, "       t=1 v=%.3e\n", buf[0]/(num));
      fprintf(fp, "   end\n");
      fprintf(fp, "end\n");

  }
}

    }

    }
  }
}
  

  // do it again, this time we define the junctions
  S.Push( root_int );  // the implcit assumption is that there are no branches at the root!

  int num_ed;
  double maf1,maf0, ratio0, ratio1;
  while (1) {
    rc = S.Pop(vid);
    if ( rc == -1 ) break;
    num_ed = e_stack.Pop( e_dn );
    
#if PRINT_CRAP    
    printf("%d : %s\n", cntr++, _Vertices[ vid ]->Name());
#endif

    for (int jj=0; jj<_Vertices[ vid ]->GetNumFanin(); jj++) {
      if ( jj>=2 ) {
	printf("Extra upper edge at node %s discarded\n", _Vertices[vid]->Name());
	break;
      }
      e = _Vertices[vid]->GetIthFanin(jj);
      S.Push( e->Up()->Id() );
      e_stack.Push( e );
      //printf("e=%s, epush=%d\n", _Vertices[vid]->Name(), e->Up()->Id());
    }
    
    if ( fp ) {
      // root does not have a "down" edge!
      if ( strcmp(_Vertices[vid]->Name(),root)==0 ) { //0 == _Vertices[vid]->GetNumFanout() ) {
	fprintf(fp,"#### root node %s_%s\n", e->Name(), _Vertices[ vid ]->Name());
	fprintf(fp,"def boundarycondition\n");
	fprintf(fp,"  location=%s_%s type=depth\n", e->Name(), _Vertices[ vid ]->Name()); // need
	fprintf(fp,"  def timeseries\n");
	fprintf(fp,"    timeunit=hour\n");  // might need to change
	fprintf(fp,"    t=0 v=%.2f\n", downstream_BC);      // need to add
	fprintf(fp,"    t=100 v=%.2f\n",downstream_BC);      // need to add
	fprintf(fp,"  end\n");
	fprintf(fp,"end\n");
      }
      if ( num_ed > 0 ) {
	if ( 0 == _Vertices[vid]->GetNumFanin() ) {
	  e = _Vertices[vid]->GetIthFanout(0);      // we don't allow multiple fanout
						    // anyway

    int NVAR = runoff->Return_Nvar();
    fprintf(fp, "############################################ NVAR = %i\n\n\n\n\n", NVAR);


  // edit by Justin 20200324
	int tmp = atoi(e->Name());
	rc = QSOURCE->Query(tmp,buf);
	for (int ii=0; ii<QSOURCE->NumTime();ii++) {
	  if ( buf[ii]< upstream_BC ) buf[ii]= upstream_BC;
	  rc = -1;
	}
	fprintf(fp,"#### leaf node %s_%s %c\n", e_dn->Name(), _Vertices[ vid ]->Name(),
	  (rc<0?'F':'T'));
	fprintf(fp,"def qsource\n");
	fprintf(fp,"  location=%s_%s\n", e_dn->Name(), _Vertices[ vid ]->Name());
	fprintf(fp,"  def timeseries\n");
	fprintf(fp,"    timeunit=hour\n");  // might need to change

	total_q += buf[0];
  double qsource_tol = 0.0;
  double temp_qsource = buf[0];
  int qsource_repeat_checker = 0;
	fprintf(fp,"    t=0.0 v=%.3e\n",buf[0]);      // need to add
	for (int jj=0; jj<QSOURCE->NumTime(); jj++){
        if (jj >= 1 && jj < QSOURCE->NumTime()-1){
          if (abs(buf[jj] - temp_qsource) >= qsource_tol){
            if (qsource_repeat_checker == 0){
              fprintf(fp,"    t=%d v=%.3e\n", jj, buf[jj-1]);
              fprintf(fp,"    t=%d v=%.3e\n", jj+1, buf[jj]);
              qsource_repeat_checker = 1;
            }
            else if (qsource_repeat_checker ==1){fprintf(fp,"    t=%d v=%.3e\n", jj+1, buf[jj]);}
            temp_qsource = buf[jj];
          }
          else{qsource_repeat_checker = 0;}
        }
        else { 
          fprintf(fp,"    t=%d v=%.3e\n", jj+1, buf[jj]);
          qsource_repeat_checker = 1;}
  }
	fprintf(fp,"  end\n");
	fprintf(fp,"end\n");
}


	if ( 1 == _Vertices[vid]->GetNumFanin() ) 
  {
	  /*fprintf(fp, "def junction\n");
	  fprintf(fp, "  down=%s_%s up1=%s_%s coeff1=%.2f\n",
		  e_dn->Name(), _Vertices[ vid ]->Name(),
		  _Vertices[ vid ]->GetIthFanin(0)->Name(), _Vertices[ vid ]->Name(),
		  1.0
		  );
	  fprintf(fp, "end\n");
  */
  } 

  else if ( 2 == _Vertices[ vid ]->GetNumFanin() ) 
  {
	  fprintf(fp, "def junction\n");
	  maf0 = _Vertices[ vid ]->GetIthFanin(0)->MAF();
	  maf1 = _Vertices[ vid ]->GetIthFanin(1)->MAF();
	  ratio0 = maf0/(maf1+maf0);
	  ratio1 = maf1/(maf1+maf0);

	  if ( (maf0+maf1) > 2.0 ) {
	    double thresh= 0.01;//original value = 0.01
	    if ( ratio0 < thresh) {
	      ratio0=thresh;
	      ratio1=1.0-thresh;
	    } else if ( ratio1 < thresh) {
	      ratio1=thresh;
	      ratio0=1.0-thresh;
	    }
	  } else {
	    double thresh = 0.01;
	    if ( ratio0 < thresh) {
	      ratio0=thresh;
	      ratio1=1.0-thresh;
	    } else if ( ratio1 < thresh) {
	      ratio1=thresh;
	      ratio0=1.0-thresh;
	    }
	  }

	  fprintf(fp, "  down=%s_%s up1=%s_%s coeff1=%.3f up2=%s_%s, coeff2=%.3f\n",
		  e_dn->Name(), _Vertices[ vid ]->Name(),
		  _Vertices[ vid ]->GetIthFanin(0)->Name(), _Vertices[ vid ]->Name(),
		  ratio0,
		  _Vertices[ vid ]->GetIthFanin(1)->Name(), _Vertices[ vid ]->Name(),
		  ratio1
		  );
	  fprintf(fp, "end\n");
	  
	} else {
	}
}
}
}



  

  printf("Total Lateral = %.3e, Total Q = %.3e\n", total_lat, total_q);

  if (buf) free(buf);
  return 0;

}


// Local Variables:
// mode: c++
// End:
