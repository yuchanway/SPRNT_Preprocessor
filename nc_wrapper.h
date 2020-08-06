/* 
  wrapper for netcdf stored surface runoff values
*/

#ifndef _NC_WRAPPER_H
#define _NC_WRAPPER_H

class NC {
protected:
  int            _num_comid;
  int            _num_time;
  int 			     _num_vars;
  double         _delta_t;

  int           *_comids;
  float         *_runoff;
  
  double        *_buffer;

  int            _search_index(int x);

public:
  NC():_num_comid(0),_num_time(0),_delta_t(-1.0) {
    _comids = NULL;
    _runoff = NULL;
  }
  ~NC() {
    if (_comids) free (_comids);
    if (_runoff) free(_runoff);
    //if (_num_vars) free (_num_vars);
  }
  
  double& DeltaT() { return _delta_t;}
  int NumComid() const { return _num_comid; }
  int NumTime() const { return _num_time; }
  int Return_Nvar() { return _num_vars; }  

  int Init(const char *fname);

  int Query(int where, double *buffer);
  
  void Dump(FILE *F, int n_col);

};

#endif

// Local Variables:
// mode: c++
// End:
