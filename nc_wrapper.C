/* methods */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

#include "nc_wrapper.h"

#define ERRCODE 0
#define NETCDF_ERR(e) {printf("Error: %s\n", nc_strerror(e)); return(ERRCODE);}

#define PRT 1

/* search for the index of given x, return -1 if not found
   stupid linear search since the values are not sorted
   but (hopefully) the array is not too big */
int NC::_search_index(int x) {
  int jj;
  for (jj=0; jj<_num_comid; jj++) {
    if ( _comids[jj] == x ) return (jj);
  } 
  return (-1);
}

int NC::Init(const char *fname) {
  FILE *F;
  int   jj,rc;
  int   ncid, ndims, dimids[128];
  int   nvars, varids[128], vardim;
  size_t alen, len[128];

  F=fopen(fname,"r");
  if ( !F ) {
    fprintf(stderr, "Bummer: unable to open %s to read. \n", fname);
    return (-1);
  }
  fclose(F);

  /* read the netCDF file  */
  if ( (rc = nc_open(fname, NC_NOWRITE, &ncid)) != NC_NOERR) NETCDF_ERR(rc);

  
  if ( (rc = nc_inq_dimids(ncid, &ndims, dimids,0)) != NC_NOERR) NETCDF_ERR(rc);
#if PRT
  //int kk,where;
  //printf("dim = %d\n", ndims);
#endif
  for (jj=0; jj<ndims;jj++) { 
    if ( (rc = nc_inq_dimlen(ncid, dimids[jj],&alen)) != NC_NOERR) NETCDF_ERR(rc);
    len[jj] = alen;

#if PRT
    //printf("The %d th dimension's length is: %d \n", dimids[jj], (unsigned int)alen);
#endif
  }

  _num_comid = (int)len[0];
  _num_time = (int)len[1];
  
  if ( (rc = nc_inq_varids(ncid, &nvars, varids)) != NC_NOERR) NETCDF_ERR(rc);
#if PRT
  printf("var = %d\n", nvars);
#endif

  _num_vars = nvars;

  for (jj=0; jj<nvars;jj++) {
    if ( (rc = nc_inq_varndims(ncid, varids[jj],&vardim)) != NC_NOERR) NETCDF_ERR(rc);
#if PRT
    //printf("%d %d\n", varids[jj], vardim);
#endif
  }

#if PRT
  //printf("%d\n", (unsigned int)len[0]);
#endif



  _comids = (int*)malloc( len[0] *sizeof(int) ); //allocate memory for comids
  _runoff = (float*)malloc( len[0]*len[1] * sizeof(float));  //assign memory for runoff



  if (!_comids) return (-1);
  //assign comid to memory
  if ( ( rc = nc_get_var_int(ncid, varids[0], _comids)) != NC_NOERR) NETCDF_ERR(rc);
  
  if ( !_runoff ) return (-1);

  if ( (rc=nc_get_var_float(ncid, varids[1], _runoff)) != NC_NOERR) NETCDF_ERR(rc);

  if ( (rc=nc_close(ncid)) != NC_NOERR ) NETCDF_ERR(rc);

//printf("dimensions : %d , %d \n", len[0],len[1]);
//for (int kk = 0; kk<_num_comid; kk++){printf("%i \n", _comids[kk]);}
/*
#if PRT  
  for (kk=0; kk< len[0]; kk++) {
    where = _search_index(_comids[kk]);
    if ( where >= 0 ) {
      printf("%d: ", _comids[kk]);
      for ( jj=0; jj<len[1]; jj++) {
	printf("%g ", _runoff[jj*len[0]+where]); 
      }
      printf("\n");
      
    } else {
      printf("%d not found\n", _comids[kk]);
    }
  }
#endif
*/
#if 0
  *e_comids = comids;
  *e_sgs = bgs;
  *e_ss = ss;
#endif

  return rc;

}

void NC::Dump(FILE *F, int num_col) {
  int nn;
  int kk, where;

  nn = _num_time < num_col ? _num_time : num_col;

  for (kk=0; kk< _num_comid; kk++) {
    where = _search_index(_comids[kk]);
    if ( where >= 0 ) {
      printf("%d: ", _comids[kk]);
      for (int jj=0; jj<nn; jj++) {
	printf("%g ", _runoff[jj*_num_comid + where]);
      }
      printf("\n");
    } else {
      printf("%d not found\n", _comids[kk]);
    }
  }
  
}
int NC::Query(int comid, double *buffer) {
  int w;
  int rc=1;
  w = _search_index( comid );
  if ( w >=0 ) {
    for (int jj=0; jj<_num_time; jj++) buffer[jj] = _runoff[jj*_num_comid+w]; 
  } else {
    rc = -1;
    for (int jj=0; jj<_num_time; jj++) buffer[jj] = 0.0;
  }
  //printf("%i, %i\n", comid, w);
  return rc;
}



// Local Variables:
// mode: c++
// End:
