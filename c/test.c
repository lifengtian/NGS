#include <stdio.h>
#include <stdlib.h>

int main(int argc, char ** argv){
   if (argc<2) {
      fprintf(stderr, "usage: example <path-to-job>\n");
      return 1;
   } else {
	printf("%s\t%s\tsize:%d\t%d\n", argv[0],argv[1],sizeof(char),sizeof(char *));
	}

   return 1;
}
