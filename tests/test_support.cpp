#include <catch2/catch_test_macros.hpp>
#include <cstring>

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "expos/support.hpp"


// RAII wrapper for a bcf1_t with two alleles set and unpacked.
struct VarHelper {
    bcf_hdr_t* hdr;
    bcf1_t*    v;
    VarHelper(const char* ref, const char* alt) {
        hdr = bcf_hdr_init("w");
        (void)bcf_hdr_sync(hdr);
        v = bcf_init();
        const char* als[] = {ref, alt};
        bcf_update_alleles(hdr, v, als, 2);
        bcf_unpack(v, BCF_UN_STR);
    }
    ~VarHelper() { bcf_destroy(v); bcf_hdr_destroy(hdr); }
};


// RAII wrapper for a simple all-M bam1_t + bam_pileup1_t.
// seq must contain only uppercase A/C/G/T/N.
struct ReadHelper {
    bam1_t*       b;
    bam_pileup1_t pl;
    ReadHelper(const char* seq, int32_t qpos, int indel = 0) {
        b = bam_init1();
        auto l = strlen(seq);
        uint32_t cig = bam_cigar_gen(static_cast<int>(l), BAM_CMATCH);
        bam_set1(b, 1, "r", 0, 0, 0, 0, 1, &cig, 0, 0, 0, l, seq, nullptr, 0);
        pl      = {};
        pl.b    = b;
        pl.qpos = qpos;
        pl.indel = indel;
    }
    ~ReadHelper() { bam_destroy1(b); }
};


// RAII wrapper for a bam1_t with an explicit CIGAR + bam_pileup1_t.
// Used for insertion tests where the CIGAR includes BAM_CINS ops.
struct IndelReadHelper {
    bam1_t*       b;
    bam_pileup1_t pl;
    IndelReadHelper(
        const char*     seq,
        int             l_seq,
        const uint32_t* cigar,
        int             n_cigar,
        int32_t         qpos,
        int             indel = 0
    ) {
        b = bam_init1();
        bam_set1(b, 1, "r", 0, 0, 0, 0,
                 static_cast<size_t>(n_cigar), cigar,
                 0, 0, 0,
                 static_cast<size_t>(l_seq), seq, nullptr, 0);
        pl      = {};
        pl.b    = b;
        pl.qpos = qpos;
        pl.indel = indel;
    }
    ~IndelReadHelper() { bam_destroy1(b); }
};


TEST_CASE ("eval_support_snp") {
    VarHelper var("G", "A");

    SECTION ("read base matches alt") {
        ReadHelper rd("TACGT", 1);  // base at qpos 1 is 'A'
        REQUIRE (eval_support_snp(var.v, &rd.pl));
    }

    SECTION ("read base is ref, not alt") {
        ReadHelper rd("TGCGT", 1);  // base at qpos 1 is 'G' (ref)
        REQUIRE_FALSE (eval_support_snp(var.v, &rd.pl));
    }

    SECTION ("read base is neither ref nor alt") {
        ReadHelper rd("TCCGT", 1);  // base at qpos 1 is 'C'
        REQUIRE_FALSE (eval_support_snp(var.v, &rd.pl));
    }

    SECTION ("invalid qpos (deleted/refskipped base)") {
        ReadHelper rd("TACGT", -1);
        REQUIRE_FALSE (eval_support_snp(var.v, &rd.pl));
    }
}


TEST_CASE ("eval_support_mnp") {
    VarHelper var("GG", "AC");

    SECTION ("all bases match alt") {
        ReadHelper rd("ACGT", 0);  // bases at qpos 0,1 are 'A','C'
        REQUIRE (eval_support_mnp(var.v, &rd.pl));
    }

    SECTION ("first base matches, second does not") {
        ReadHelper rd("AGGT", 0);  // 'A' ok, 'G' != 'C'
        REQUIRE_FALSE (eval_support_mnp(var.v, &rd.pl));
    }

    SECTION ("no bases match alt") {
        ReadHelper rd("GGGT", 0);  // "GG" != "AC"
        REQUIRE_FALSE (eval_support_mnp(var.v, &rd.pl));
    }

    SECTION ("invalid qpos") {
        ReadHelper rd("ACGT", -1);
        REQUIRE_FALSE (eval_support_mnp(var.v, &rd.pl));
    }
}


TEST_CASE ("eval_support_del") {
    // deletion of 3 bases: get_var_indel_len = 1 - 4 = -3
    VarHelper var("ACGT", "A");

    SECTION ("indel length matches deletion") {
        ReadHelper rd("ACGT", 0, -3);
        REQUIRE (eval_support_del(var.v, &rd.pl));
    }

    SECTION ("wrong deletion length") {
        ReadHelper rd("ACGT", 0, -2);
        REQUIRE_FALSE (eval_support_del(var.v, &rd.pl));
    }

    SECTION ("insertion in pileup instead of deletion") {
        ReadHelper rd("ACGT", 0, 3);
        REQUIRE_FALSE (eval_support_del(var.v, &rd.pl));
    }

    SECTION ("no indel in pileup") {
        ReadHelper rd("ACGT", 0, 0);
        REQUIRE_FALSE (eval_support_del(var.v, &rd.pl));
    }

    SECTION ("invalid qpos") {
        ReadHelper rd("ACGT", -1, -3);
        REQUIRE_FALSE (eval_support_del(var.v, &rd.pl));
    }
}


TEST_CASE ("eval_support_ins") {
    // ref="G", alt="GAC" — insertion of "AC" (len 2) after anchor "G"
    VarHelper var("G", "GAC");

    SECTION ("insertion length and bases match") {
        // CIGAR: 1M 2I 3M — sequence: G(anchor) AC(inserted) ATG(match)
        // bam_plp_insertion requires p->indel > 0 to proceed
        uint32_t cigar[] = {
            bam_cigar_gen(1, BAM_CMATCH),
            bam_cigar_gen(2, BAM_CINS),
            bam_cigar_gen(3, BAM_CMATCH)
        };
        IndelReadHelper rd("GACATG", 6, cigar, 3, 0, 2);
        REQUIRE (eval_support_ins(var.v, &rd.pl));
    }

    SECTION ("insertion bases mismatch (correct length, wrong bases)") {
        // CIGAR: 1M 2I 3M — inserted "TT" != "AC"
        uint32_t cigar[] = {
            bam_cigar_gen(1, BAM_CMATCH),
            bam_cigar_gen(2, BAM_CINS),
            bam_cigar_gen(3, BAM_CMATCH)
        };
        IndelReadHelper rd("GTTATG", 6, cigar, 3, 0, 2);
        REQUIRE_FALSE (eval_support_ins(var.v, &rd.pl));
    }

    SECTION ("insertion length mismatch") {
        // CIGAR: 1M 1I 3M — only 1 inserted base, variant expects 2
        uint32_t cigar[] = {
            bam_cigar_gen(1, BAM_CMATCH),
            bam_cigar_gen(1, BAM_CINS),
            bam_cigar_gen(3, BAM_CMATCH)
        };
        IndelReadHelper rd("GATCG", 5, cigar, 3, 0, 1);
        REQUIRE_FALSE (eval_support_ins(var.v, &rd.pl));
    }

    SECTION ("no insertion in pileup (p->indel == 0)") {
        // bam_plp_insertion gates on p->indel > 0 and returns 0
        ReadHelper rd("GATCG", 0, 0);
        REQUIRE_FALSE (eval_support_ins(var.v, &rd.pl));
    }

    SECTION ("invalid qpos") {
        ReadHelper rd("GATCG", -1);
        REQUIRE_FALSE (eval_support_ins(var.v, &rd.pl));
    }
}
