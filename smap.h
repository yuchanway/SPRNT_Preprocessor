/*
 a map to match strings to integers. Insertion only, no deletions
 plus a few inquiry

 anticipated number of unique entries 100 ~ 100,000

 use unordered_map 
  
 interface using char* to be backward compatible

 capacity: since we use int, the max number of entries is ~2,000,000,000

 */

#ifndef _SMAP_H
#define _SMAP_H

#include <string>
#include <unordered_map>

using namespace std;

class SMap {
private:
  int                             _cntr;
  unordered_map < string, int >   _M;

public:
  SMap():_cntr(0) { _M.rehash(15000000); } // slightly different hash function
  ~SMap() {}

  /* public methods */
  // return the size of the mapping table
  int Size() const { return ( _M.size() ); }

  // check, if the entry is found returns the key
  // if not, create a new entry and returns key
  int Create(const char *s) {
    string ns(s);
    unordered_map < string, int >::const_iterator loc;

    loc = _M.find(ns);
    if ( loc == _M.end() ) {
      _M[ns] = _cntr;
      return (_cntr++);
    } else {
      return _M[ ns ];
    }
  }

  // check only
  // if not found, return -1
  // if found, return key
  int Check(const char *s) {
    string ns(s);
    unordered_map < string, int >::const_iterator loc;

    loc = _M.find(ns);
    if ( loc == _M.end() ) {
      return (-1);
    } else {
      return _M[ ns ];
    }
  }
  
};

#endif

// Local Variables: 
// mode: c++
// End:


