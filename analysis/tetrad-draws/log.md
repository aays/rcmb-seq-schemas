
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
```
