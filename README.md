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
./target/release/RECallS -i "${sample}.bam"
```

## Outputs
Within the output directory (`-o`):
```
recalls
├── calls.bed
└── chrom_indel_stats.tsv
```

Where:
|name|desc|
|-|-|
|calls.bed|Putative intrachromosomal inversion and deletion events|
|chrom_indel_stats.tsv|Chromosome read indel stats for primary and secondary alignments|

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
-n "test/amplicon/chm13_chr7_lcr.bed"
```
