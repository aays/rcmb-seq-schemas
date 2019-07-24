
## 23/07/2019

to do -
1. create fastas from Liu et al tetrad data
2. parse through fastas, 'drawing' with different read pair schemas
3. estimate how many COs/gene conversions will be caught by each

eventually - simulate much larger amounts of 'gene converted' F1 sequence to run read pair code on

parental set H (2935 x 2936) - SRR5243250 -> SRR5243297

moving the gvcfs over:

```bash
mkdir -p data/tetrad-draws/vcfs
cd data/tetrad-draws/vcfs
ln -sv ../../../../../meiotic_mutation/Data\ /SRR*.g.vcf.gz .
ln -sv ../../../../../meiotic_mutation/Data\ /SRR*.g.vcf.gz.tbi .
```

let's start off with just chromosome 1:

```bash
time ./bin/vcf2fasta.py \
-v data/tetrad-draws/vcfs/SRR5243250.g.vcf.gz \
-r data/references/chlamy_reference.fasta \
-i chromosome_1:1-8033585 \
--min_GQ 30 > data/tetrad-draws/fastas/SRR5243250.chr01.fasta
```

took 4 min - there are 48 in total - should queue these up in parallel

```bash
parallel -j 20 -i sh -c 'time ./bin/vcf2fasta.py \
-v data/tetrad-draws/vcfs/SRR52432{}.g.vcf.gz \
-r data/references/chlamy_reference.fasta \
-i chromosome_1:1-8033585 \
--min_GQ 30 > data/tetrad-draws/fastas/SRR52432{}.chr01.fasta' -- {51..70}

parallel -j 20 -i sh -c 'time ./bin/vcf2fasta.py \
-v data/tetrad-draws/vcfs/SRR52432{}.g.vcf.gz \
-r data/references/chlamy_reference.fasta \
-i chromosome_1:1-8033585 \
--min_GQ 30 > data/tetrad-draws/fastas/SRR52432{}.chr01.fasta' -- {71..97}
```

## 24/07/2019

other things to also do -
1. create VCF containing just CC2935 and CC2936, chr1
2. create fastas for parental strains?

```bash
mkdir -p data/tetrad-draws/bams
cd data/tetrad-draws/bams
ln -sv /scratch/research/data/chlamydomonas/species_wide/Individual.InDelRealigned.BAMs/CC2935* .
ln -sv /scratch/research/data/chlamydomonas/species_wide/Individual.InDelRealigned.BAMs/CC2936* .
cd ../../../

# use GATK to create VCF
ln -sv /opt/gatk/GenomeAnalysisTK.jar bin/

# get reference dict
ln -sv /scratch/research/references/chlamydomonas/5.3_chlamy_w_organelles_mt_minus/chlamy.5.3.w_organelles_mtMinus.dict chlamy_reference.dict

java -jar bin/GenomeAnalysisTK.jar \
-R data/references/chlamy_reference.fasta \
-T HaplotypeCaller \
-I data/tetrad-draws/bams/CC2935.RG.bam \
-I data/tetrad-draws/bams/CC2936.RG.bam \
-L chromosome_1 \
-o data/tetrad-draws/vcfs/parental.vcf
# took 4 hours...

bgzip data/tetrad-draws/vcfs/parental.vcf
tabix -p vcf data/tetrad-draws/vcfs/parental.vcf.gz
```

filter for SNPs that differ b/w the two strains - either there are 2 ALTs 
(1 for each strain) or one of the strains has the REF allele and the other the ALT

```bash
time python3.5 analysis/tetrad-draws/parental_vcf_filter.py \
--vcf data/tetrad-draws/vcfs/parental.vcf.gz \
--out data/tetrad-draws/vcfs/parental_filtered.vcf
```

this retained just about half of the VCF:

```bash
Filtering complete.
105347 of 212163 records retained.
```

next up - run different read pair schemes through this file
to get upper bound of identifiable COs/NCOs





















