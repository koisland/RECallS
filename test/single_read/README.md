# Test data
From Schierup et al. 2025 and https://github.com/PeterSoerud/recombination_calling/tree/main/tables/classified_reads

## GCV
Gene-conversion event
```bash
samtools view -bh /uufs/chpc.utah.edu/common/home/u1643401/projects/sperm_ampl_sv/results/align/CT22.bam ENA_CBCUDK010000046_CBCUDK010000046.1:46,027,777-46,027,932 -o workflow/scripts/RECall/test/CT22_ENA_CBCUDK010000046_CBCUDK010000046.1_46027777-46027932.bam
```

## CO
Crossover event
```bash
samtools view -bh /uufs/chpc.utah.edu/common/home/u1643401/projects/sperm_ampl_sv/results/align/CT22.bam ENA_CBCUDK010000011_CBCUDK010000011.1:6,335,921-6,341,074 -o workflow/scripts/RECall/test/CT22_ENA_CBCUDK010000011_CBCUDK010000011.1_6335921-6341074.bam
```

## Visualizing mismatch PDF 
```bash
python test/single_read_vis_mismatch_dist.py \
-i test/amplicon/chm13_chr7_sim.bam \
-r chr7:6031981-6085629 \
-l 30580116-9da2-f84f-b0c8-46b140978a73 b919a3e8-c866-4793-60d4-c8d69792fe9f 
```
