#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>

//This function was copied from stackoverflow
//users Carl Norum and Basic
void _mkdir(const char *dir) {
  char tmp[256];
  char *p = NULL;
  size_t len;
 
  snprintf(tmp, sizeof(tmp),"%s",dir);
  len = strlen(tmp);
  if(tmp[len - 1] == '/')
    tmp[len - 1] = 0;
     for(p = tmp + 1; *p; p++){
       if(*p == '/') {
         *p = 0;
         mkdir(tmp, S_IRWXU);
         *p = '/';
       }
       mkdir(tmp, S_IRWXU);
     }
}


