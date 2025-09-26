# Small-scale somatic mutation detecting
Location: small_mutation_detection.sh

# Small-scale somatic mutation filtering
## 1. mutations present in the parental genomes
Location: mutation_present_in_parental_genomes.sh

## 2. mutations with VAFs greater than 0.25 in all samples 
Location: VAF_greater_than_0.25.sh

## 3. mutations present in more than one family
Location: mutation_present_in_more_than_one_family.sh

## 4. clustered mutations occurring within 150 bp of each other
Location: mutation_cluster_within_150bp.sh