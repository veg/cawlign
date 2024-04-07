#include "stringBuffer.h"
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <string.h>

#define SWAP(A, B, C)                                                          \
  ({                                                                           \
    (C) = (A);                                                                 \
    (A) = (B);                                                                 \
    (B) = (C);                                                                 \
  })

using namespace std;

long StringBuffer::sbDefaultLength = 16, StringBuffer::sbDefaultBoost = 16,
     Vector::vDefaultLength = 16, Vector::vDefaultBoost = 16,
     VectorFP::vDefaultLength = 16, VectorFP::vDefaultBoost = 16;

/*----------------------------------------------------------------------------------------------------
 */

StringBuffer::StringBuffer(void) {
  sLength = 0;
  saLength = StringBuffer::sbDefaultLength;
  sData = (char *)malloc(sizeof(char) * (saLength + 1));
  sData[0] = 0;
}

/*----------------------------------------------------------------------------------------------------
 */

void StringBuffer::swap(StringBuffer &src) {
  long t;
  SWAP(sLength, src.sLength, t);
  SWAP(saLength, src.saLength, t);
  char *tc;
  SWAP(sData, src.sData, tc);
}

/*----------------------------------------------------------------------------------------------------
 */

void StringBuffer::flip (void) {
    long half = sLength >> 1;
    const long lm1 = sLength - 1L;
    for (long i = 0L; i < half; i++) {
        char t = sData [i];
        sData [i] = sData [lm1 - i];
        sData [lm1 - i] = t;
    }
}

/*----------------------------------------------------------------------------------------------------
 */

StringBuffer::~StringBuffer(void) { free(sData); }

/*----------------------------------------------------------------------------------------------------
 */

void StringBuffer::appendChar(const char c) {
  long addThis;
  if (saLength == sLength) {
    addThis = saLength / 8;
    if (StringBuffer::sbDefaultBoost > addThis)
      addThis = StringBuffer::sbDefaultBoost;
    saLength += addThis;
    sData = (char *)realloc(sData, sizeof(char) * (saLength + 1));
  }
  sData[sLength] = c;
  sData[++sLength] = 0;
}

/*----------------------------------------------------------------------------------------------------
 */

void StringBuffer::appendBuffer(const char *buffer, const long length) {
  long addThis, pl = length > 0 ? length : strlen(buffer);

  if (pl > 0) {
    if (saLength < sLength + pl) {
      addThis = saLength / 8;

      if (StringBuffer::sbDefaultBoost > addThis)
        addThis = StringBuffer::sbDefaultBoost;
      if (addThis < pl)
        addThis = pl;

      saLength += addThis;
      sData = (char *)realloc(sData, sizeof(char) * (saLength + 1));
    }
    for (addThis = 0; addThis < pl; addThis++)
      sData[sLength++] = buffer[addThis];

    sData[sLength] = 0;
  }
}

/*----------------------------------------------------------------------------------------------------
 */

void StringBuffer::resetString(void) {
  sLength = 0;
  sData[sLength] = 0;
}

/*----------------------------------------------------------------------------------------------------
 */

Vector::Vector(void) {
  vLength = 0;
  vaLength = Vector::vDefaultLength;
  vData = (long *)calloc(vaLength, sizeof(long));
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::swap(Vector &src) {
  long t;
  SWAP(vLength, src.vLength, t);
  SWAP(vaLength, src.vaLength, t);
  long *tc;
  SWAP(vData, src.vData, tc);
}

/*----------------------------------------------------------------------------------------------------
 */

Vector::~Vector(void) {
    free(vData);
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::appendValue(const long l) {
  long addThis;
  if (vLength == vaLength) {
    addThis = vaLength / 8;
    if (Vector::vDefaultBoost > addThis)
      addThis = Vector::vDefaultBoost;
    vaLength += addThis;
    vData = (long *)realloc(vData, sizeof(long) * vaLength);
  }
  vData[vLength++] = l;
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::appendVector(const Vector &v) {
  for (unsigned long k = 0UL; k < v.vLength; k++) {
    appendValue(v.vData[k]);
  }
}
/*----------------------------------------------------------------------------------------------------
 */

void Vector::storeValue(const long v, const unsigned long l) {
  long addThis;
  if (l >= vaLength) {
    addThis = l - vaLength + 1;
    if (Vector::vDefaultBoost > addThis)
      addThis = Vector::vDefaultBoost;
    vaLength += addThis;
    vData = (long *)realloc(vData, sizeof(long) * vaLength);
    vLength = l + 1;
  }
  vData[l] = v;
}

/*----------------------------------------------------------------------------------------------------
 */

long Vector::value(const long idx) const {
  /*if (idx >= vLength) {
    cerr << "Indexing past the end of a vector" << endl;
    exit (1);
  }*/
  return vData[idx];
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::storeVector(const Vector &v, const unsigned long l) {
  if (l < vLength && vData[l])
    delete (Vector *)(vData[l]);

  storeValue((long)&v, l);
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::remove(const unsigned long l) {
  if (l < vLength) {
    for (unsigned long k = l + 1; k < vLength; k++) {
      vData[k - 1] = vData[k];
    }
    vLength--;
  }
}

/*----------------------------------------------------------------------------------------------------
 */

long Vector::extractMin(VectorFP &values) {
  cawlign_fp current_min = DBL_MAX;
  long best_index = -1;

  for (unsigned long i = 0; i < vLength; i++) {
    cawlign_fp try_value = values.value(vData[i]);
    if (try_value < current_min) {
      best_index = i;
      current_min = try_value;
    }
  }

  if (best_index >= 0) {
    current_min = vData[best_index];
    remove(best_index);
    return current_min;
  }

  return best_index;
}

/*----------------------------------------------------------------------------------------------------
 */

int long_comp(const void *v1, const void *v2) {
  long l1 = *(long *)v1, l2 = *(long *)v2;
  if (l1 < l2)
    return -1;
  if (l1 > l2)
    return 1;
  return 0;
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::sort(void) {
  qsort((void *)vData, vLength, sizeof(long), long_comp);
}

/*----------------------------------------------------------------------------------------------------
 */

void Vector::resetVector(void) { vLength = 0; }

/*----------------------------------------------------------------------------------------------------
 */

VectorFP::VectorFP(void) {
  vLength = 0;
  vaLength = Vector::vDefaultLength;
  vData = (cawlign_fp *)calloc(vaLength, sizeof(cawlign_fp));
}

/*----------------------------------------------------------------------------------------------------
 */

VectorFP::~VectorFP(void) { free(vData); }

/*----------------------------------------------------------------------------------------------------
 */

void VectorFP::appendValue(const cawlign_fp l) {
  long addThis;
  if (vLength == vaLength) {
    addThis = vaLength / 8;
      if (VectorFP::vDefaultBoost > addThis)
      addThis = VectorFP::vDefaultBoost;
    vaLength += addThis;
    vData = (cawlign_fp *)realloc(vData, sizeof(cawlign_fp) * vaLength);
  }
  vData[vLength++] = l;
}

/*----------------------------------------------------------------------------------------------------
 */

void VectorFP::appendValues(const cawlign_fp* l, long N) {
    long addThis;
    
    if (vLength + N > vaLength) {
        addThis =  vLength + N + 1 - vaLength;
        long def_inc = vaLength / 8;
        if (def_inc > addThis) {
            addThis = def_inc;
        }
        if (VectorFP::vDefaultBoost > addThis)
            addThis = VectorFP::vDefaultBoost;
        vaLength += addThis;
        vData = (cawlign_fp *)realloc(vData, sizeof(cawlign_fp) * vaLength);
    }
    
    memcpy (vData + vLength, l, N*sizeof (cawlign_fp));
    vLength += N;
}

/*----------------------------------------------------------------------------------------------------
 */

void VectorFP::storeValue(const cawlign_fp v, const unsigned long l) {
  long addThis;
  if (l >= vaLength) {
    addThis = l - vaLength + 1;
    if (VectorFP::vDefaultBoost > addThis)
      addThis = VectorFP::vDefaultBoost;
    vaLength += addThis;
    vData = (cawlign_fp *)realloc(vData, sizeof(cawlign_fp) * vaLength);
    vLength = l + 1;
  }
  vData[l] = v;
}
