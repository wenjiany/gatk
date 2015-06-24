/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.variant;

import htsjdk.variant.variantcontext.Allele;

/**
 * This class contains any constants (primarily FORMAT/INFO keys) in VCF files used by the GATK.
 * Note that VCF-standard constants are in VCFConstants, in htsjdk.  Keys in header lines should
 * have matching entries in GATKVCFHeaderLines
 */
public final class GATKVCFConstants {

    //INFO keys
    public static final String ALLELE_BALANCE_HET_KEY =             "ABHet";
    public static final String ALLELE_BALANCE_HOM_KEY =             "ABHom";
    public static final String ORIGINAL_AC_KEY =                    "AC_Orig"; //SelectVariants
    public static final String BEAGLE_AC_COMP_KEY =                 "ACH"; //BeagleOutputToVCF
    public static final String ORIGINAL_AF_KEY =                    "AF_Orig"; //SelectVariants
    public static final String BEAGLE_AF_COMP_KEY =                 "AFH"; //BeagleOutputToVCF
    public static final String ORIGINAL_AN_KEY =                    "AN_Orig"; //SelectVariants
    public static final String BEAGLE_AN_COMP_KEY =                 "ANH"; //BeagleOutputToVCF
    public static final String BASE_COUNTS_KEY =                    "BaseCounts";
    public static final String BASE_QUAL_RANK_SUM_KEY =             "BaseQRankSum";
    public static final String GENOTYPE_AND_VALIDATE_STATUS_KEY =   "callStatus";
    public static final String CLIPPING_RANK_SUM_KEY =              "ClippingRankSum";
    public static final String CULPRIT_KEY =                        "culprit";
    public static final String SPANNING_DELETIONS_KEY =             "Dels";
    public static final String ORIGINAL_DP_KEY =                    "DP_Orig"; //SelectVariants
    public static final String DOWNSAMPLED_KEY =                    "DS";
    public static final String FISHER_STRAND_KEY =                  "FS";
    public static final String GC_CONTENT_KEY =                     "GC";
    public static final String GQ_MEAN_KEY =                        "GQ_MEAN";
    public static final String GQ_STDEV_KEY =                       "GQ_STDDEV";
    public static final String HAPLOTYPE_SCORE_KEY =                "HaplotypeScore";
    public static final String HI_CONF_DENOVO_KEY =                 "hiConfDeNovo";
    public static final String HOMOPOLYMER_RUN_KEY =                "HRun";
    public static final String HARDY_WEINBERG_KEY =                 "HW";
    public static final String AVG_INTERVAL_DP_KEY =                "IDP"; //DiagnoseTargets
    public static final String INTERVAL_GC_CONTENT_KEY =            "IGC";
    public static final String INBREEDING_COEFFICIENT_KEY =         "InbreedingCoeff";
    public static final String LIKELIHOOD_RANK_SUM_KEY =            "LikelihoodRankSum";
    public static final String LO_CONF_DENOVO_KEY =                 "loConfDeNovo";
    public static final String LOW_MQ_KEY =                         "LowMQ";
    public static final String MLE_ALLELE_COUNT_KEY =               "MLEAC";
    public static final String MLE_ALLELE_FREQUENCY_KEY =           "MLEAF";
    public static final String MLE_PER_SAMPLE_ALLELE_COUNT_KEY =    "MLPSAC";
    public static final String MLE_PER_SAMPLE_ALLELE_FRACTION_KEY = "MLPSAF";
    public static final String MAP_QUAL_RANK_SUM_KEY =              "MQRankSum";
    public static final String MENDEL_VIOLATION_LR_KEY =            "MVLR";
    public static final String NOCALL_CHROM_KEY =                   "NCC";
    public static final String NUMBER_OF_DISCOVERED_ALLELES_KEY =   "NDA";
    public static final String NEGATIVE_LABEL_KEY =                 "NEGATIVE_TRAIN_SITE";
    public static final String NUM_GENOTYPES_CHANGED_KEY =          "NumGenotypesChanged"; //BeagleOutputToVCF
    public static final String NON_DIPLOID_RATIO_KEY =              "OND";
    public static final String ORIGINAL_ALT_ALLELE_INFO_KEY =       "OriginalAltAllele"; //BeagleOutputToVCF
    public static final String ORIGINAL_CONTIG_KEY =                "OriginalChr"; //LiftoverVariants
    public static final String ORIGINAL_START_KEY =                 "OriginalStart"; //LiftoverVariants
    public static final String N_BASE_COUNT_KEY =                   "PercentNBase";
    public static final String RBP_INCONSISTENT_KEY =               "PhasingInconsistent"; //ReadBackedPhasing
    public static final String GENOTYPE_PRIOR_KEY =                 "PG";
    public static final String POSITIVE_LABEL_KEY =                 "POSITIVE_TRAIN_SITE";
    public static final String QUAL_BY_DEPTH_KEY =                  "QD";
    public static final String BEAGLE_R2_KEY =                      "R2"; //BeagleOutputToVCF
    public static final String READ_POS_RANK_SUM_KEY =              "ReadPosRankSum";
    public static final String REFSAMPLE_DEPTH_KEY =                "REFDEPTH";
    public static final String REPEATS_PER_ALLELE_KEY =             "RPA";
    public static final String REPEAT_UNIT_KEY =                    "RU";
    public static final String SAMPLE_LIST_KEY =                    "Samples";
    public static final String STRAND_ODDS_RATIO_KEY =              "SOR";
    public static final String STR_PRESENT_KEY =                    "STR";
    public static final String TRANSMISSION_DISEQUILIBRIUM_KEY =    "TDT";
    public static final String VARIANT_TYPE_KEY =                   "VariantType";
    public static final String VQS_LOD_KEY =                        "VQSLOD";

<<<<<<< HEAD
    //Custom INFO Keys
    public static final String PolyX =                              "PolyX";
    public static final String DP4T =                               "DP4T";
    public static final String FLANKING =                           "FLANKING";
    public static final String MMQ =                                "MMQ";

=======
>>>>>>> 616e6962c6d5644ef80fa3258ee561beaf294c87
    //FORMAT keys
    public static final String ALLELE_BALANCE_KEY =                 "AB";
    public static final String PL_FOR_ALL_SNP_ALLELES_KEY =         "APL";
    public static final String RBP_HAPLOTYPE_KEY =                  "HP"; //ReadBackedPhasing
    public static final String AVG_INTERVAL_DP_BY_SAMPLE_KEY =      "IDP"; //DiagnoseTargets
    public static final String JOINT_LIKELIHOOD_TAG_NAME =          "JL"; //FamilyLikelihoodsUtils
    public static final String JOINT_POSTERIOR_TAG_NAME =           "JP"; //FamilyLikelihoodsUtils
    public static final String LOW_COVERAGE_LOCI =                  "LL"; //DiagnoseTargets
    public final static String MIN_DP_FORMAT_KEY =                  "MIN_DP";
    public static final String MAPPING_QUALITY_ZERO_BY_SAMPLE_KEY = "MQ0";
    public static final String ORIGINAL_GENOTYPE_KEY =              "OG"; //BeagleOutputToVCF
    public static final String HAPLOTYPE_CALLER_PHASING_GT_KEY =    "PGT";
    public static final String HAPLOTYPE_CALLER_PHASING_ID_KEY =    "PID";
    public static final String PHRED_SCALED_POSTERIORS_KEY =        "PP"; //FamilyLikelihoodsUtils / PosteriorLikelihoodsUtils
    public static final String REFERENCE_GENOTYPE_QUALITY =         "RGQ";
    public static final String STRAND_COUNT_BY_SAMPLE_KEY =         "SAC";
    public static final String STRAND_BIAS_BY_SAMPLE_KEY =          "SB";
    public final static String TRANSMISSION_PROBABILITY_KEY =       "TP"; //PhaseByTransmission
    public static final String ZERO_COVERAGE_LOCI =                 "ZL"; //DiagnoseTargets

    //FILTERS
    /* Note that many filters used throughout GATK (most notably in VariantRecalibration) are dynamic,
       their names (or descriptions) depend on some threshold.  Those filters are not included here
     */
    public static final String BEAGLE_MONO_FILTER_NAME =            "BGL_SET_TO_MONOMORPHIC";
    public static final String LOW_QUAL_FILTER_NAME =               "LowQual";

    // Symbolic alleles
    public final static String SYMBOLIC_ALLELE_DEFINITION_HEADER_TAG = "ALT";
    public final static String NON_REF_SYMBOLIC_ALLELE_NAME = "NON_REF";
    public final static Allele NON_REF_SYMBOLIC_ALLELE = Allele.create("<"+NON_REF_SYMBOLIC_ALLELE_NAME+">", false); // represents any possible non-ref allele at this site
    public final static String SPANNING_DELETION_SYMBOLIC_ALLELE_NAME = "*:DEL";
    public final static Allele SPANNING_DELETION_SYMBOLIC_ALLELE = Allele.create("<"+SPANNING_DELETION_SYMBOLIC_ALLELE_NAME+">", false); // represents any possible spanning deletion allele at this site
}
