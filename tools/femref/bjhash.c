/*
  Taken from http://burtleburtle.net/bob/hash/evahash.html
  Hash Functions for Hash Table Lookup
  Robert J. Jenkins Jr., 1995-1997
*/
/* $Id: bjhash.c,v 1.5 2005/05/19 00:46:54 mstorti Exp $ */

#include <stdio.h>

typedef  unsigned long int  u4;   /* unsigned 4-byte type */
typedef  unsigned     char  u1;   /* unsigned 1-byte type */

/* The mixing step */
#define mix(a,b,c) \
{ \
  a=a-b;  a=a-c;  a=a^(c>>13); \
  b=b-c;  b=b-a;  b=b^(a<<8);  \
  c=c-a;  c=c-b;  c=c^(b>>13); \
  a=a-b;  a=a-c;  a=a^(c>>12); \
  b=b-c;  b=b-a;  b=b^(a<<16); \
  c=c-a;  c=c-b;  c=c^(b>>5);  \
  a=a-b;  a=a-c;  a=a^(c>>3);  \
  b=b-c;  b=b-a;  b=b^(a<<10); \
  c=c-a;  c=c-b;  c=c^(b>>15); \
}

/* The whole new hash function */
u4 bjhash( k, length, initval)
register u1 *k;        /* the key */
u4           length;   /* the length of the key in bytes */
u4           initval;  /* the previous hash, or an arbitrary value */
{
  register u4 a,b,c;  /* the internal state */
  u4          len;    /* how many key bytes still need mixing */

  /* Set up the internal state */
  len = length;
  a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
  c = initval;         /* variable initialization of internal state */

  /*---------------------------------------- handle most of the key */
  while (len >= 12)
    {
      a=a+(k[0]+((u4)k[1]<<8)+((u4)k[2]<<16) +((u4)k[3]<<24));
      b=b+(k[4]+((u4)k[5]<<8)+((u4)k[6]<<16) +((u4)k[7]<<24));
      c=c+(k[8]+((u4)k[9]<<8)+((u4)k[10]<<16)+((u4)k[11]<<24));
      mix(a,b,c);
      k = k+12; len = len-12;
    }

  /*------------------------------------- handle the last 11 bytes */
  c = c+length;
  switch(len)              /* all the case statements fall through */
    {
    case 11: c=c+((u4)k[10]<<24);
    case 10: c=c+((u4)k[9]<<16);
    case 9 : c=c+((u4)k[8]<<8);
      /* the first byte of c is reserved for the length */
    case 8 : b=b+((u4)k[7]<<24);
    case 7 : b=b+((u4)k[6]<<16);
    case 6 : b=b+((u4)k[5]<<8);
    case 5 : b=b+k[4];
    case 4 : a=a+((u4)k[3]<<24);
    case 3 : a=a+((u4)k[2]<<16);
    case 2 : a=a+((u4)k[1]<<8);
    case 1 : a=a+k[0];
      /* case 0: nothing left to add */
    }
  mix(a,b,c);
  /*-------------------------------------------- report the result */
  // printf("state %d, k %d, hash %d\n",initval,*k,c);
  return c;
}
