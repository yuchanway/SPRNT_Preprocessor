/* 
   utilities: vector and stack

 */

#ifndef _DAG_H
#define _DAG_H

#include <string.h>
#include <memory.h>
#include <assert.h>

#define DEFAULT_APPSIZE 512
#define MAXFIXSZ 1023

/* 
   Vector with flexible length
   resize() will change the size, and wipe out the existing content
   grow() will grow, but also keep the existing content
*/

template<class T, int FIXSZ=MAXFIXSZ> 
class FlexVec {
protected:
  T   _fixmem[FIXSZ];
  int _asz;
  T*  _vec;
  
  //if not large enough allocates more than requested and recopies the old data into new memory
  int do_grow( int sz ) { 
    T* newvec = new T[sz + _asz]; // at least doubles the allocated size
    for( int k=0; k!=_asz; ++k) newvec[k] = _vec[k];
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = newvec;
    _asz += sz;                            
    return _asz;
  }
  
public:
  /* constructors and destructor */
  FlexVec() : _asz(FIXSZ), _vec(_fixmem) {};
  FlexVec( int sz ) { 
    if (sz <= FIXSZ) {
      _vec = _fixmem;
      _asz = FIXSZ;
    } else {
      _vec = new T[sz];
      _asz = sz;
    }
  };
  
 // disallow copy constructor
  FlexVec( const FlexVec<T,FIXSZ>& ) { assert(0); } 
  
  ~FlexVec() {
    if (_vec != (T*)_fixmem) delete[] _vec; 
  };
  
  /* access methods */
  operator T*() { return _vec; }
  T& operator[](int i) { assert(i<_asz); assert(i>=0); return _vec[i]; }; 
  
  // if not large enough, allocate exactly as required
  int size( int sz ) { 
    if ( sz <= _asz ) return _asz; 
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = new T[sz];                
    _asz = sz;                            
    return _asz;
  }
  
  // resize will grow the memory allocation, but will wipe out the content
  int resizep( int sz ) { 
    if ( sz <= _asz ) return _asz; 
    if (_vec != (T*)_fixmem) delete[] _vec;
    _vec = new T[sz + _asz]; 
    _asz += sz;                            
    return _asz;
  }
  
  inline int grow (int sz) { return (sz<_asz)?_asz:do_grow(sz); }
  
};

// a simple stack, for graph traversal
template<class T> 
class MyStack {
protected:
  T                               *_array;
  int                              _size;
  int                              _current;
  
public:
 MyStack():_size(DEFAULT_APPSIZE),_current(0) {
   _array = new T [DEFAULT_APPSIZE];
 }
  ~MyStack() {
    if (_array) delete [] _array;
  }
  
  void Push(T data) {
    int newsz;
    T *tmp;
    if ( _current >= _size ) {  // grow and copy
      newsz = 2*_size;
      tmp = new T[newsz];
      memcpy(tmp, _array, _size*sizeof(T));
      delete[] _array;
      _array = tmp;
      _size = newsz;
    }
    _array[_current++] = data;
  }
  
  // rc == 1, good
  // rc == -1, top of the stack
  int Pop(T& data) {  
    if ( _current > 0 ) {
      data = _array[--_current];
      return 1;
    } else {
      return -1;
    }
  }
  
  int StackSize(void) const { return _current; }

  // flush out whatever in the stack
  void Flush(void) { _current = 0; }
};

#endif

// Local Variables:
// mode: c++
// End:
