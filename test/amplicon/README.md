# Data
Extract CHM13 chromosome 7.
```bash
samtools faidx /scratch/ucgd/lustre-labs/vollger/users/Keith/EldeLabRotation/data/annot/human/chm13v2.0.fa.gz chr7 \
    | bgzip > test/amplicon/chm13_chr7.fa.gz
```

Get low-complexity regions.
```bash
longdust test/amplicon/chm13_chr7.fa.gz > test/amplicon/chm13_chr7_lcr.bed
```

Create inversion. Allow events not overlapping low-complexity regions greater than 500 bp.
```bash
python test/amplicon/denovo_inv_create.py \
-f test/amplicon/chm13_chr7.fa.gz \
-s 42 \
-pl 50_000 \
-n <(awk '$3-$2>500' test/amplicon/chm13_chr7_lcr.bed) \
-o test/amplicon/chm13_chr7_sim_inv
```

Create deletion.
```bash
python test/amplicon/denovo_inv_create.py \
-f test/amplicon/chm13_chr7.fa.gz \
-s 120 \
-pl 50_000 \
-n <(awk '$3-$2>500' test/amplicon/chm13_chr7_lcr.bed) \
-e deletion \
-o test/amplicon/chm13_chr7_sim_del
```

Normal reads
```bash
badread simulate \
    --reference test/amplicon/chm13_chr7.fa.gz \
    --quantity 5x \
    --error_model pacbio2021 \
    --qscore_model pacbio2021 \
    --identity 30,3 \
    --seed 42 \
    | bgzip > test/amplicon/chm13_chr7_sim_norm.fq.gz
```

Simulated reads from event
```bash
for file in test/amplicon/chm13_chr7_sim_del_after_event.fa test/amplicon/chm13_chr7_sim_inv_after_event.fa; do
    badread simulate \
        --reference test/amplicon/chm13_chr7_sim_del_after_event.fa \
        --quantity 1x \
        --error_model pacbio2021 \
        --qscore_model pacbio2021 \
        --identity 30,3 \
        --seed 42 \
        | bgzip >> test/amplicon/chm13_chr7_sim_events.fq.gz
done
```

```bash
minimap2 -ax lr:hq test/amplicon/chm13_chr7.fa.gz test/amplicon/chm13_chr7_sim_events.fq.gz -t8 \
    | samtools sort -u - \
    | samtools view -bh -o test/amplicon/chm13_chr7_sim.bam

samtools index test/amplicon/chm13_chr7_sim.bam
```
