/*
 **********************************************************************
 ** md5driver.c -- sample routines to test                           **
 ** RSA Data Security, Inc. MD5 message digest algorithm.            **
 ** Created: 2/16/90 RLR                                             **
 ** Updated: 1/91 SRD                                                **
 **********************************************************************
 */

/*
 **********************************************************************
 ** Copyright (C) 1990, RSA Data Security, Inc. All rights reserved. **
 **                                                                  **
 ** RSA Data Security, Inc. makes no representations concerning      **
 ** either the merchantability of this software or the suitability   **
 ** of this software for any particular purpose.  It is provided "as **
 ** is" without express or implied warranty of any kind.             **
 **                                                                  **
 ** These notices must be retained in any copies of any part of this **
 ** documentation and/or software.                                   **
 **********************************************************************
 */

#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "./md5.h"

int hash(int *w,int n) {
  MD5_CTX ctx;
  MD5Init(&ctx);
  MD5Update(&ctx,w,n*sizeof(int));
  MD5Final(&ctx);
  int val;
  memcpy(&val,ctx.digest,sizeof(int));
  return val;
}

int main() {
  int v[] = {0,1,2,3,4,5,6,7};
  int val = hash(v,8);
  printf("hash %x\n",val);
}
