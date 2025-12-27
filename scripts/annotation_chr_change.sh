#!/bin/bash

reflines=$(wc -l resources/chr.accesssion.list.txt | awk '{ print $1 }')
refcurrentline=1
while [ "$refcurrentline" -le "$reflines" ]
do
  refchr=$(sed "${refcurrentline}q;d" resources/chr.accesssion.list.txt | awk '{ print $1 }')
  refnc=$(sed "${refcurrentline}q;d" resources/chr.accesssion.list.txt | awk '{ print $2 }')
  sed -i "s/${refnc}/${refchr}/g" resources/EquCab3.0_annotation_chr_only.gff
  refcurrentline=$(($refcurrentline+1))
done
