#!/usr/bin/env bash

# $1 -> louvain output
# $2 -> Solo nodes in the network
# $3 -> output

hierarchy="$1"
snodes="$2"
output="$3"

echo $hierarchy
louvain-hierarchy -n "$hierarchy"
louvain-hierarchy -n "$hierarchy" | tail -n1 > $hierarchy.highestlevel

level_id=`louvain-hierarchy -n $hierarchy | tail -n1 | cut -d\  -f2 | cut -d: -f1`
nclust=`louvain-hierarchy -n $hierarchy | tail -n1 | cut -d\  -f3`

 # Get the highest level cluster
nnodes=`louvain-hierarchy -l $level_id $hierarchy  | wc -l`
nsolo=`cat $snodes | sed '/^\s*$/d' | wc -l`

louvain-hierarchy -n $hierarchy
echo "Chosen level: $level_id"

 # Format into MCL output format
awk -v nnodes="$nnodes" -v nclust="$nclust" -v nsolo="$nsolo" '
 BEGIN{
   printf "#Louvain wrapper\n(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mclmatrix\nbegin\n", nnodes+nsolo, nclust+nsolo
   split("",idx)
   split("",sidx)
   sidx_i = 0
 }
 {

   if ( FNR == NR ) {
     split($0,line," ")
     node_id = line[1]
     cluster = line[2]

     if (cluster in idx) {
       idx[cluster] = idx[cluster] " " node_id
     } else {
       idx[cluster] = cluster " " node_id
     }
   } else {
     sidx[sidx_i] = $0
     sidx_i += 1
   }

 }
 END{
   for (i=0; i < nclust; i++) {
     printf "%s $\n", idx[i]
   }
   for(i=0; i < nsolo; i++) {
     printf "%d %s $\n", nclust+i, sidx[i]
   }
   printf " \n)"
 }' <(louvain-hierarchy -l $level_id $hierarchy) <(cat $snodes) > $output

