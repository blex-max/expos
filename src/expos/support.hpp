#pragma once

#include <htslib/sam.h>
#include <htslib/vcf.h>

bool eval_support_snp (const bcf1_t* m, const bam_pileup1_t* p);

bool eval_support_mnp (const bcf1_t* m, const bam_pileup1_t* p);

bool eval_support_ins (const bcf1_t* m, const bam_pileup1_t* p);

bool eval_support_del (const bcf1_t* m, const bam_pileup1_t* p);
