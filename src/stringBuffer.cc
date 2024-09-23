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

/**
 * Initializes an empty `StringBuffer` with a default initial capacity.
 * This buffer dynamically grows as new characters are appended.
 */
StringBuffer::StringBuffer(void) {
  sLength = 0;
  saLength = StringBuffer::sbDefaultLength;
  sData = (char *)malloc(sizeof(char) * (saLength + 1));
  sData[0] = 0;
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Swaps the contents of this `StringBuffer` with another `StringBuffer`.
 *
 * This function exchanges the data, length, and capacity of two `StringBuffer` objects.
 *
 * @param src The `StringBuffer` object to swap with.
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

/**
 * Reverses the content of the `StringBuffer`.
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

/**
 * Frees the memory allocated for the string buffer.
 */
StringBuffer::~StringBuffer(void) { free(sData); }

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Appends a single character to the end of the buffer, growing the buffer if needed.
 *
 * @param c The character to append.
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

/**
 * Appends a string or a buffer of specified length to the `StringBuffer`.
 *
 * @param buffer The string or character buffer to append.
 * @param length The length of the buffer, if known. If not, the length is inferred using `strlen`.
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

/**
 * Resets the `StringBuffer` to an empty state.
 *
 * Clears the buffer content by resetting its length, but keeps the allocated memory.
 */
void StringBuffer::resetString(void) {
  sLength = 0;
  sData[sLength] = 0;
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Initializes an empty vector of `long` values with a default initial capacity.
 */
Vector::Vector(void) {
  vLength = 0;
  vaLength = Vector::vDefaultLength;
  vData = (long *)calloc(vaLength, sizeof(long));
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Swaps the contents of this `Vector` with another `Vector`.
 *
 * Exchanges the data, length, and capacity of two `Vector` objects.
 *
 * @param src The `Vector` object to swap with.
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

/**
 * Frees the memory allocated for the vector data.
 */
Vector::~Vector(void) {
    free(vData);
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Appends a value to the `Vector`.
 *
 * Adds a `long` value to the end of the vector, growing the vector if needed.
 *
 * @param l The `long` value to append.
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

/**
 * Appends the contents of another `Vector` to this `Vector`.
 *
 * @param v The source `Vector` whose contents are to be appended.
 */
void Vector::appendVector(const Vector &v) {
  for (unsigned long k = 0UL; k < v.vLength; k++) {
    appendValue(v.vData[k]);
  }
}
/*----------------------------------------------------------------------------------------------------
 */

/**
 * Stores a value at a specific index in the `Vector`.
 *
 * If the index is beyond the current capacity, the vector grows to accommodate the value.
 *
 * @param v The `long` value to store.
 * @param l The index at which to store the value.
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

/**
 * Retrieves a value from the `Vector` at the specified index.
 *
 * @param idx The index from which to retrieve the value.
 * @return The value stored at the specified index.
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

/**
 * Stores a pointer to a `Vector` at a specified index, growing the vector if necessary.
 *
 * @param v The `Vector` object to store.
 * @param l The index at which to store the vector.
 */
void Vector::storeVector(const Vector &v, const unsigned long l) {
  if (l < vLength && vData[l])
    delete (Vector *)(vData[l]);

  storeValue((long)&v, l);
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Removes the element at the given index and shifts the remaining elements to fill the gap.
 *
 * @param l The index of the element to remove.
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

/**
 * Extracts the minimum value from the vector based on a `VectorFP` of floating-point values.
 *
 * @param values The `VectorFP` of floating-point values to compare.
 * @return The index of the minimum value, or `-1` if the vector is empty.
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

/**
 * Used by `qsort` to compare two `long` values.
 *
 * @param v1 Pointer to the first value.
 * @param v2 Pointer to the second value.
 * @return -1 if `v1` is less than `v2`, 1 if greater, and 0 if equal.
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

/**
 * Sorts the `Vector` in ascending order.
 */
void Vector::sort(void) {
  qsort((void *)vData, vLength, sizeof(long), long_comp);
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Resets the `Vector` to an empty state.
 * 
 * Clears the vector by resetting its length, but keeps the allocated memory.
 */
void Vector::resetVector(void) { vLength = 0; }

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Initializes an empty vector of floating-point values with a default initial capacity.
 */
VectorFP::VectorFP(void) {
  vLength = 0;
  vaLength = Vector::vDefaultLength;
  vData = (cawlign_fp *)calloc(vaLength, sizeof(cawlign_fp));
}

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Frees the memory allocated for the floating-point vector data.
 */
VectorFP::~VectorFP(void) { free(vData); }

/*----------------------------------------------------------------------------------------------------
 */

/**
 * Appends a floating-point value to the `VectorFP`.
 *
 * Adds a floating-point value to the end of the vector, growing the vector if needed.
 *
 * @param l The floating-point value to append.
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

/**
 * Appends multiple floating-point values to the `VectorFP`.
 *
 * @param l The array of floating-point values to append.
 * @param N The number of values to append.
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

/**
 * Stores a floating-point value at a specific index in the `VectorFP`.
 *
 * If the index is beyond the current capacity, the vector grows to accommodate the value.
 *
 * @param v The floating-point value to store.
 * @param l The index at which to store the value.
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
