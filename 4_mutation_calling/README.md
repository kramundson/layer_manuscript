Mutation calling steps:

1. Short and long read alignments to the RP assembly (variant discovery)
2. Variant filtering, applied in two steps
 a. first pass filter with BCFtools - good for fast rough cuts but not inter-sample comparisons
 b. second-pass filter in R - good for intersample comparisons
 c. multiple rounds of targeted variant calling as new samples were added
 d. generated raw data for manuscript figures
