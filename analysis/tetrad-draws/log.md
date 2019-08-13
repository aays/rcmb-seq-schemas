
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

## 28/7/2019

today: read pair schemes to get upper bound of identifiable COs/NCOs

script should take in:
1. VCF (obv)
2. read length L
3. insert size b/w reads (I) (though technically 'insert size' includes reads...)
4. chromosome

from start of chromosome - 
1. count SNPs in window of length L AND window at L + I (left and right pair) - use tabix for this
2. move window forward by L/2 and repeat

do this over L = {50, 100, 150, 250}
and I = {50, 100, 150, 250, 2000, 3000, 4000, 5000} (short insert and mate pair)

```bash
mkdir -p data/tetrad-draws/snp-counts

# test run
time python3.5 analysis/tetrad-draws/read_snp_counts.py \
--vcf data/tetrad-draws/vcfs/parental_filtered.vcf.gz \
--read_length 150 \
--insert_size 50 \
--chrom chromosome_1 \
--out test.out
```

looks great, and only took ~3 minutes - let's queue this up for all of
the combinations above in parallel:

```bash
parallel -j 4 -i sh -c \
'for insert_size in 50 100 150 250 2000 3000 4000 5000; do \
time python3.5 analysis/tetrad-draws/read_snp_counts.py \
--vcf data/tetrad-draws/vcfs/parental_filtered.vcf.gz \
--read_length {} --insert_size ${insert_size} --chrom chromosome_1 \
--out data/tetrad-draws/snp-counts/r{}_i${insert_size}.tsv; done' -- 50 100 150 250
```

once these are done, will need to calculate SNP density on reads,
and assess how many pairs are 'callable' (ie at least two SNPs *per read*)
and how many are 'highly informative' (kind of arbitrary - say five SNPs per read?)

after that: creating a VCF with all the Liu samples and randomly drawing 'mate pairs'
using the best schema above to see how many are genuinely callable

update: longest took 10 min - genomewide run should take
about 1.5 hours

## 7/8/2019

to do:
- redo SNP density counts with more permissive parameters (ie at least 2 SNPs total)
- write script that detects COs and NCOs assuming all sites are callable
    - also simulate data where this is the case?
    - could also use AT or 12 system to designate parental haps - probably easier

```bash
mkdir -p analysis/rcmb-calls
mkdir -p data/rcmb-calls
```

creating 'all callable' sites:
- create two 10 kb fasta records
    - parent 1 - all As
    - parent 2 - all Ts
- create simulated 'tetrad' from cross b/w parent 1 and 2
    - walk along chromosome and determine crossover/NCO locations using Poisson process
    - update four tetrad strings accordingly
    - could be faster if done in 'long' format before transposing?
- do 'read draws' from single individuals in tetrad
    - what SNP combinations indicate one form of recombination over another?
    

creating two 10 kb 'parents':

```python
>>> outname = 'data/rcmb-calls/parents_all_callable.fasta'
>>> d = dict.fromkeys(['parent1', 'parent2'], '')
>>> d['parent1'] = ''.join(['A' for i in range(10000)])
>>> d['parent2'] = ''.join(['T' for i in range(10000)])
>>> with open(outname, 'w') as f:
  2     for id, seq in d.items():
  3         f.write('>' + id + '\n')
  4         f.write(seq + '\n')
```

## 9/8/2019

'recombining' the parents

some fixed rcmb events -
- COs
    - 200 kb
    - 800 kb - this could be CO-GC
- NCOs
    - 500 kb - 100 bp long
    - 700 kb - 50 bp long

script should encode 'instructions' for where phase changes happen instead
of storing large strings in memory

at start - create four 'offspring' - assign 2 to parent 1 (A) and
2 to parent 2 (T)

- if CO
    - draw one 'A' offspring and one 'T' offspring and switch to other phase
- if CO-GC
    - encode phase change _n_ bases downstream of initial 'switch point' _s_
    - in distance b/w _s_ and _n_ draw one offspring and switch phase for that length
- if NCO-GC
    - draw one offspring and temporarily switch phase for _n_ bases


















