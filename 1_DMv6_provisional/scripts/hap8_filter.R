
library(tidyverse)

vcf_file <- "../firstpass-filter.vcf.gz"

vcf_header <- system(paste9("zgrep '^#C' ", vcf_file), intern = T) %>% 
  str_remove("^#") %>% 
  str_split(pattern = "\t") %>% 
  pluck(1)

vcf <- read_tsv(vcf_file,
                col_names = vcf_header,
                comment = "#",
                show_col_types = F)

sample_vars <- colnames(vcf)[10:ncol(vcf)]
attr <- unique(vcf$FORMAT) %>% 
  str_split(pattern = ":") %>% 
  pluck(1)
attr

vcf_wide <- vcf %>% 
  separate_wider_delim(col = sample_vars, attr, delim = ":")

# filter RP depth >= 40 and <= 140
# show example discovery figures - 12.5% VAF is an appropriate cutoffs
# filter L1-specific if at least 4 of MF43, 63, 86, 93 and 105 have 12.5% VAF
# filter if each of MF11, 74, 75, 83, 84, 98 and 104 have <12.5% VAF

# write out candidate hap8 variants as BED3.
# need to do two things with this: intersect firstpass calls with this BED3
# then use the BED3 to define region input and the VCF to define allele input
# for genotyping the MF93 and BB populations
# with that VCF, move on to the asm and annotation workflows