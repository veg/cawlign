#ifndef __STRINGBUFFER__
#define __STRINGBUFFER__

#include "alignment.h"

//__________________________________________________________________________________________

/**
 * @brief A simple string buffer class that allows for dynamic resizing.
 * 
 */
class StringBuffer {

  char *sData;
  unsigned long sLength, saLength;

public:
  /**
   * @brief Construct a new String Buffer object
   * 
   */
  StringBuffer(void);
  /**
   * @brief Destroy the String Buffer object
   * 
   */
  ~StringBuffer(void);

  /**
   * @brief Get the underlying C string.
   * 
   * @return char* A pointer to the underlying C string.
   */
  char *getString(void) const { return sData; }
  /**
   * @brief Append a character to the string buffer.
   * 
   * @param c The character to append.
   */
  void appendChar(const char);
  /**
   * @brief Append a C string to the string buffer.
   * 
   * @param buffer The C string to append.
   * @param length The length of the C string to append. If -1, the entire string is appended.
   */
  void appendBuffer(const char *, const long = -1);
  /**
   * @brief Reset the string buffer to an empty state.
   * 
   */
  void resetString(void);
  /**
   * @brief Swap the contents of this string buffer with another.
   * 
   * @param other The other string buffer to swap with.
   */
  void swap(StringBuffer &);
  /**
   * @brief Get the length of the string buffer.
   * 
   * @return unsigned long The length of the string buffer.
   */
  unsigned long length(void) const { return sLength; }
  /**
   * @brief Reset the length of the string buffer.
   * 
   * @param newL The new length of the string buffer.
   */
  void reset_length(unsigned long newL) {
    if (newL < sLength) {
      sLength = newL;
    }
  }

  /**
   * @brief Set a character at a specific index.
   * 
   * @param i The index to set the character at.
   * @param c The character to set.
   * @return char The character that was previously at the index.
   */
  char setChar(const long i, const char c) {
    char oc = sData[i];
    sData[i] = c;
    return oc;
  }
    
  /**
   * @brief Get a character at a specific index.
   * 
   * @param i The index to get the character from.
   * @return char The character at the index.
   */
  char getChar(const long i) const { return sData[i]; }
  /**
   * @brief Flip the string buffer in place.
   * 
   */
  void flip (void);
  /**
   * @brief Detach the underlying C string from the string buffer.
   * 
   */
  void detach (void) { sData = nullptr;}

  static long sbDefaultLength, sbDefaultBoost;
};

//__________________________________________________________________________________________

/**
 * @brief A simple vector class for floating point numbers that allows for dynamic resizing.
 * 
 */
class VectorFP {

  cawlign_fp *vData;

  unsigned long vLength, vaLength;

public:
  /**
   * @brief Construct a new VectorFP object
   * 
   */
  VectorFP(void);
  /**
   * @brief Destroy the VectorFP object
   * 
   */
  ~VectorFP(void);

  /**
   * @brief Append a value to the vector.
   * 
   * @param value The value to append.
   */
  void appendValue(const cawlign_fp);
  /**
   * @brief Append multiple values to the vector.
   * 
   * @param values The values to append.
   * @param count The number of values to append.
   */
  void appendValues(const cawlign_fp*, long);
  /**
   * @brief Store a value at a specific index.
   * 
   * @param value The value to store.
   * @param index The index to store the value at.
   */
  void storeValue(const cawlign_fp, const unsigned long);
  /**
   * @brief Get a value at a specific index.
   * 
   * @param idx The index to get the value from.
   * @return cawlign_fp The value at the index.
   */
  cawlign_fp value(const long idx) { return vData[idx]; }
  /**
   * @brief Get the length of the vector.
   * 
   * @return unsigned long The length of the vector.
   */
  unsigned long length(void) const { return vLength; }
  /**
   * @brief Get the underlying C array.
   * 
   * @return const cawlign_fp* A pointer to the underlying C array.
   */
  const cawlign_fp * values (void) {return vData;}
  /**
   * @brief Get the underlying C array.
   * 
   * @return cawlign_fp* A pointer to the underlying C array.
   */
  cawlign_fp * rvalues (void) {return vData;}

  static long vDefaultLength, vDefaultBoost;
};

//__________________________________________________________________________________________

/**
 * @brief A simple vector class for long integers that allows for dynamic resizing.
 * 
 */
class Vector {

  long *vData;

  unsigned long vLength, vaLength;

public:
  /**
   * @brief Construct a new Vector object
   * 
   */
  Vector(void);
  /**
   * @brief Destroy the Vector object
   * 
   */
  ~Vector(void);

  /**
   * @brief Append a value to the vector.
   * 
   * @param value The value to append.
   */
  void appendValue(const long);
  /**
   * @brief Append another vector to this vector.
   * 
   * @param other The other vector to append.
   */
  void appendVector(const Vector &);
  /**
   * @brief Extract the minimum value from a VectorFP object.
   * 
   * @param other The VectorFP object to extract the minimum value from.
   * @return long The minimum value.
   */
  long extractMin(VectorFP &);
  /**
   * @brief Reset the vector to an empty state.
   * 
   */
  void resetVector(void);
  /**
   * @brief Remove a value at a specific index.
   * 
   * @param index The index to remove the value at.
   */
  void remove(const unsigned long);
  /**
   * @brief Store a value at a specific index.
   * 
   * @param value The value to store.
   * @param index The index to store the value at.
   */
  void storeValue(const long, const unsigned long);
  /**
   * @brief Store another vector at a specific index.
   * 
   * @param other The other vector to store.
   * @param index The index to store the vector at.
   */
  void storeVector(const Vector &, const unsigned long);
  /**
   * @brief Sort the vector in place.
   * 
   */
  void sort(void);
  /**
   * @brief Swap the contents of this vector with another.
   * 
   * @param other The other vector to swap with.
   */
  void swap(Vector &);
  /**
   * @brief Get a value at a specific index.
   * 
   * @param idx The index to get the value from.
   * @return long The value at the index.
   */
  long value(const long idx) const;
  /**
   * @brief Get the length of the vector.
   * 
   * @return unsigned long The length of the vector.
   */
  unsigned long length(void) const { return vLength; }
    
  const long * rvalues (void) const {return vData;}

  static long vDefaultLength, vDefaultBoost;
};

#endif
