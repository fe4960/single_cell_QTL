#!/bin/sh
#19_D003
#for s in 19_D005  19_D006  19_D007  19_D008  #19_D009  19_D010  19_D011  19_D019  D005_13  D009_13  D013_13  D017_13  D018_13  D019_13  D021_13  D026_13  D027_13  D028_13  D030_13
#for s in 19_D009  19_D010  19_D011  19_D019  D005_13
#for s in D009_13  D013_13  D017_13  D018_13  D019_13
#for s in D021_13  D026_13  D027_13  D028_13  D030_13
for s in 19_D003
do
s2=$s
IFS='_' read -a arr <<<"$s" 
s1=${arr[0]}${arr[1]}
for cell in Rod Cone ONBC OFFBC MG AC RGC HC  
do
#echo ${s1}
#echo $s2
#done
nohup sh sc_human_retina/scripts/ASE/1phaser.sh $s1 $s2 ${cell} >  sc_human_retina/scripts/ASE/1phaser_${s1}_${cell}.out 2> sc_human_retina/scripts/ASE/1phaser_${s1}_${cell}.err &
done
done 
