#ifndef __STRINGBUFFER__
#define __STRINGBUFFER__

#include "alignment.h"

//__________________________________________________________________________________________

class StringBuffer {

  char *sData;
  unsigned long sLength, saLength;

public:
  StringBuffer(void);
  ~StringBuffer(void);

  char *getString(void) const { return sData; }
  void appendChar(const char);
  void appendBuffer(const char *, const long = -1);
  void resetString(void);
  void swap(StringBuffer &);
  unsigned long length(void) const { return sLength; }
  void reset_length(unsigned long newL) {
    if (newL < sLength) {
      sLength = newL;
    }
  }

  char setChar(const long i, const char c) {
    char oc = sData[i];
    sData[i] = c;
    return oc;
  }
    
  char getChar(const long i) const { return sData[i]; }
  void flip (void);
  void detach (void) { sData = nullptr;}

  static long sbDefaultLength, sbDefaultBoost;
};

//__________________________________________________________________________________________

class VectorFP {

  cawlign_fp *vData;

  unsigned long vLength, vaLength;

public:
  VectorFP(void);
  ~VectorFP(void);

  void appendValue(const cawlign_fp);
  void appendValues(const cawlign_fp*, long);
  void storeValue(const cawlign_fp, const unsigned long);
  cawlign_fp value(const long idx) { return vData[idx]; }
  unsigned long length(void) const { return vLength; }
  const cawlign_fp * values (void) {return vData;}
  cawlign_fp * rvalues (void) {return vData;}

  static long vDefaultLength, vDefaultBoost;
};

//__________________________________________________________________________________________

class Vector {

  long *vData;

  unsigned long vLength, vaLength;

public:
  Vector(void);
  ~Vector(void);

  void appendValue(const long);
  void appendVector(const Vector &);
  long extractMin(VectorFP &);
  void resetVector(void);
  void remove(const unsigned long);
  void storeValue(const long, const unsigned long);
  void storeVector(const Vector &, const unsigned long);
  void sort(void);
  void swap(Vector &);
  long value(const long idx) const;
  unsigned long length(void) const { return vLength; }

  static long vDefaultLength, vDefaultBoost;
};

#endif
