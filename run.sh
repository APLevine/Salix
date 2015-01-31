SCRIPTS=/home/zchads1/scripts
PEDIGREE=pedigree.pro #input pedigree file

#Generate framework
#Stage 1 (previously generate_hapdrop.py)
python $SCRIPTS/salix_stage1.py \
$PEDIGREE \
cm_lengths.txt \ #chromosome lengths in cM
0.1 31 8 #0.1 cM spacing, 31 alleles, 8 markers per tag


GENOTYPES=/home/zchads1/cluster/exome/annotate/AJ/AJ_moderate.gen #control sequence data
MERLIN=~/merlin/merlin
OUT=.
i=1

#Run Merlin
$MERLIN -p for_simulation.ped \
-d for_simulation.dat \
-m for_simulation.map \
-f for_simulation.frq \
--bits 0 --simulate -r ${i} \
--save --prefix $OUT/sim_${i} \
--swap --megabytes 1000 > $OUT/run_${i}.log

#Extract founder labelled flow
#Stage 2 (previously extract_from_simulated_multiple_v2.py)
python $SCRIPTS/salix_stage2.py \
$OUT/sim_${i}-replicate.ped \
$OUT/sim_${i}-replicate.dat \
8 \
$OUT/sim_flow_${i} >> $OUT/run_${i}.log

#Add sequence genotypes
#Stage 3 (previously hapdrop_v5_unphased.py)
#Requires mapping of sequence variants to nearest cM position (*.cm file) and sequence data ($GENOTYPES)
python $SCRIPTS/salix_stage3.py \
pedigree.pro \
../AJ_moderate.cm \
$OUT/sim_flow_${i}_maternal.txt \
$OUT/sim_flow_${i}_paternal.txt \
$GENOTYPES \
$OUT/${i} >> $OUT/run_${i}.log


