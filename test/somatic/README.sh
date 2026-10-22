#!/bin/bash

set -euo pipefail

WD=/scratch/ucgd/lustre-labs/vollger/users/Keith/RECallS/test/somatic

# assembly
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz -P "${WD}"
# BAM
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/hifi_revio_pbmay24/hg002v1.1_hifi_revio_pbmay24.bam -P "${WD}"
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/mapping/hifi_revio_pbmay24/hg002v1.1_hifi_revio_pbmay24.bam.bai -P "${WD}"

# low-complexity regions
longdust "${WD}/hg002v1.1.fasta.gz" > "${WD}/hg002v1.1_lcr.bed"

./target/release/RECallS \
-i "${WD}/hg002v1.1_hifi_revio_pbmay24.bam" \
-n "${WD}/hg002v1.1_lcr.bed" \
-o "${WD}/recalls"
