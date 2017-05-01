#!/usr/bin/env bash

# $1 -> louvain output
# $2 -> Solo nodes in the network (Only added if the gene ID is higher than the highest observed in the hierarchy)
# $3 -> output

hierarchy="$1"
snodes="$2"
output="$3"

echo $hierarchy
louvain-hierarchy -n "$hierarchy"
louvain-hierarchy -n "$hierarchy" | tail -n1 > $hierarchy.highestlevel

level_id=`louvain-hierarchy -n $hierarchy | tail -n1 | cut -d\  -f2 | cut -d: -f1`
nclust=`louvain-hierarchy -n $hierarchy | tail -n1 | cut -d\  -f3`

louvain-hierarchy -n $hierarchy
echo "Chosen level: $level_id"

 # Format into MCL output format
awk '
 BEGIN{
   split("",louvain_clusters)
   split("",louvain_ids)
   n_louvain_clusters = 0

   split("",solo_genes)
   n_solo_genes = 0

   max_gene_id = -1
   gcount = 0
 }
 {

   n = split($0,line," ")

   if ( n == 2 ) {
     node_id    = line[1]
     cluster_id = line[2]

     if ( node_id > max_gene_id ) {
       max_gene_id = node_id
     }

     gcount += 1
     if (cluster_id in louvain_clusters) {
       louvain_clusters[cluster_id] = louvain_clusters[cluster_id] " " node_id
     } else {
       louvain_clusters[cluster_id] = node_id
       louvain_ids[n_louvain_clusters] = cluster_id
       n_louvain_clusters += 1
     }
   } else if ( n == 1 && $0 > max_gene_id ) {
     gcount += 1
     solo_genes[n_solo_genes] = $0
     n_solo_genes += 1
   }

 }
 END{
   printf "#Louvain wrapper\n(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n(mclmatrix\nbegin\n", gcount, n_louvain_clusters+n_solo_genes
   for (i=0; i < n_louvain_clusters; i++) {
     printf "%d %s $\n", i, louvain_clusters[louvain_ids[i]]
   }
   for(i=0; i < n_solo_genes; i++) {
     printf "%d %s $\n", n_louvain_clusters+i, solo_genes[i]
   }
   printf " \n)"
 }' <(louvain-hierarchy -l $level_id $hierarchy) <(cat $snodes | sort | uniq) > $output

