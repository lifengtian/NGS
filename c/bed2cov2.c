#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* merge bed files
allows different count (not 0) so long as the pos is continuous
##chr1	14411	1
#chr1	14412	1
#chr1	14413	1
#chr1	14414	1
#chr1	14415	1
#chr1	14416	1

my $first=<>;
chomp($first);
my ($previous_chr,$previous_pos, $previous_count,@rest)=split/\t/,$first;
my $begin = $previous_pos;
my $end = $previous_pos;
my $count = $previous_count;

while (my $current = <>){
	chomp($current);
	my ($current_chr, $current_pos, $current_count,@rest) = split/\t/,$current;
	if ($current_chr == $previous_chr && $current_count == $previous_count && $current_pos == $previous_pos + 1 ) {
		$end = $current_pos;
		$count = $current_count;
		$previous_pos = $current_pos;
	} else {
		print join("\t",$previous_chr, $begin,$end,$count),"\n";
		$begin = $current_pos;
		$end = $current_pos;
		$count = $current_count;
		$previous_pos = $current_pos;
		$previous_count = $current_count;
		$previous_chr = $current_chr;
	}
}

print join("\t",$previous_chr, $begin,$end,$count),"\n";
*/


FILE *fp;
char *filename;
int begin, end, count = 0;
int previous_pos, current_pos, previous_count, current_count;

char current_chr[64];
char previous_chr[64];

char line[1024];

int main(int argc, char **argv){
if (argc != 2 ) {
	printf("Usage: bed2cov <bed file>\n");
	exit(-1);
}

filename = argv[1];

fp = fopen(filename, "r");
if ( fp == NULL ) {
	printf("can't open file: %s\n", filename);
	exit(-1);
}

if ( fgets(line, 1024, fp) == NULL ) {
	exit(-1);
}


 char * pch;
  	pch = strtok (line,"\t");
	strcpy(previous_chr,pch);

	pch = strtok (NULL, "\t");
	previous_pos = atoi(pch);
	pch = strtok (NULL, "\t");
	previous_count = atoi(pch);
	//printf("chr:%s pos:%d count:%d\n",previous_chr,previous_pos, previous_count);

begin = previous_pos;
end = previous_pos;
count = 0;

while ( fgets(line, 1024, fp) != NULL ) {
        pch = strtok (line,"\t");
        strcpy(current_chr,pch);

        pch = strtok (NULL, "\t");
        current_pos = atoi(pch);
        pch = strtok (NULL, "\t");
        current_count = atoi(pch);

	if ( strcmp(current_chr,previous_chr)==0 && (current_count != 0 ) && ( current_pos == previous_pos + 1) ) {
		end = current_pos;
		count += current_count;
		previous_pos = current_pos;
	} else {
		printf("%s\t%d\t%d\t%d\n",previous_chr, begin, end, count);
		begin = current_pos;
		end = current_pos;
		count = 0;
		previous_pos = current_pos;
		previous_count = current_count;
		strcpy(previous_chr, current_chr);
	}

}



	printf("%s\t%d\t%d\t%d\n",previous_chr, begin, end, count);

return 1;
}
