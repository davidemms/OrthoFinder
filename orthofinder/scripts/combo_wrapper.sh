#!/usr/bin/env bash

# $1 -> combo output
# $2 -> Solo nodes in the network (Only added if the gene ID is higher than the highest observed in the hierarchy)
# $3 -> output

combo="$1"
snodes="$2"
output="$3"

awk '
 BEGIN{
   split("",combo_clusters)
   split("",combo_ids)
   n_combo_clusters = 0

   split("",solo_genes)
   n_solo_genes = 0

   max_gene_id = -1
   gcount = 0
 }
 {

   if ( NR == FNR ) {
     
     node_id    = FNR-1
     cluster_id = $0

     if ( node_id > max_gene_id ) {
       max_gene_id = node_id
     }

     gcount += 1
     if (cluster_id in combo_clusters) {
       combo_clusters[cluster_id] = combo_clusters[cluster_id] " " node_id
     } else {
       combo_clusters[cluster_id] = node_id
       combo_ids[n_combo_clusters] = cluster_id
       n_combo_clusters += 1
     }
   } else if ( n == 1 && $0 > max_gene_id ) {
     gcount += 1
     solo_genes[n_solo_genes] = $0
     n_solo_genes += 1
   }

 }
 END{
   printf "#Combo Wrapper\n(mclheader\nmcltypematrix\ndimensions%dx%d\n)\n(mclmatrix\nbegin\n", gcount, n_combo_clusters+n_solo_genes
   for (i=0; i < n_combo_clusters; i++) {
     printf "%d %s $\n", i, combo_clusters[combo_ids[i]]
   }
   for(i=0; i < n_solo_genes; i++) {
     printf "%d %s $\n", n_combo_clusters+i, solo_genes[i]
   }
   printf " \n)"
 }' $combo $snodes > $output


