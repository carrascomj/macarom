#!/bin/bash


for file in data/fasta/*fasta; do
  imodulon=`basename "${file%.fasta}"`
  meme "$file" -dna -oc "eval/meme_results/$imodulon" -nostatus -time 14400 -mod anr -nmotifs 1 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0
done


