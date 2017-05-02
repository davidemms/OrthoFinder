#!/usr/bin/env bash

prefix="$1"
outfile="$2"

directory=`dirname $prefix`
fileprefix=`basename $prefix`


awk '
BEGIN{
  maxval = -1
}
{
  n = split($0,line,"\t")
  if (n == 3){
    if (line[0] > maxval){
      maxval = line[0]
    }
    if (line[1] > maxval){
      maxval = line[1]
    }
  }
}
END{
  printf "*Vertices %d\n*Arcs\n", maxval+1
}' `ls $directory | grep -e "^$fileprefix" | sed -e "s:^:$directory/&:"` > $outfile


