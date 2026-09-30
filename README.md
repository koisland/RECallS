# RECallS
(R)ecombination (E)vent (call)er from (S)perm long-read sequencings data

## Usage
Align reads to donor-specific assembly.
```bash
minimap -ax lr:hq --eqx -Y  -I8g "${sample}.fa.gz" "${sample}.fq.gz" \
    | samtools sort -u - \
    | samtools view -F 4 -bh -o "${sample}.bam"
```

Compile (Requires rustc >=1.98.0)
```bash
cargo build --release
```

Then run:
```bash
# Requires minimap2 in PATH
./target/release/RECallS -i "${sample}.bam" -f "${sample}.fa.gz"
```

## Outputs
Within the output directory (`-o`):
```
recalls
├── calls_del.bed
├── calls_inv.bed
├── chrom_indel_stats.tsv
└── chrom_self_align.paf
```

Where:
|name|desc|
|-|-|
|calls_del.bed|Putative intrachromosomal deletion event breakpoints|
|calls_inv.bed|Putative intrachromosomal inversion event breakpoints|
|chrom_indel_stats.tsv|Chromosome read indel stats for primary and secondary alignments|
|chrom_self_align.paf|PAF file for self-alignment. Uses `minimap2` and params: `-PD -k19 -w19 -m200`|

## Docs
See [`docs/long_distance.md`](docs/long_distance.md) for overview.

## Test
To generate simulated data in CHM13 chr7.
```bash
# Install dependencies
pixi install
snakemake -c 8 -s test/amplicon/Snakefile -np
```

Then to run:
```bash
./target/release/RECallS \
-i test/amplicon/chm13_chr7_sim.bam \
-f test/amplicon/chm13_chr7.fa.gz \
-n "test/amplicon/chm13_chr7_lcr.bed"
```
