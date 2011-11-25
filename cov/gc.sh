SM=$1

#samtools depth -b ~/cov/SureSelect50mbclean.bed $SM  > $SM.samdepth
echo "Coverage "
perl ~/cov/gc.pl $SM.samdepth
