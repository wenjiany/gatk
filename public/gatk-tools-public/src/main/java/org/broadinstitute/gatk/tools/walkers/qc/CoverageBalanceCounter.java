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

package org.broadinstitute.gatk.tools.walkers.qc;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import java.io.PrintStream;

/**
 *
 * Count strand biases by bases, based on GATKPaperGenotyper
 *
 * @author Wenjian Yang & Colton Smith
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_TOY, extraDocs = {CommandLineGATK.class} )
public class CoverageBalanceCounter extends LocusWalker<Double[], Double[]> implements TreeReducible<Double[]> {

    public static final double HUMAN_SNP_HETEROZYGOSITY = 1e-3;

    // the possible diploid genotype strings
    private static enum BASE {A, C, G, T, N};

    @Output
    private PrintStream out;

    /**
     * Return the homozygous base, return BASE.N if the site is not homozygous
     *
     */
    private BASE getHomozygousBase(final ReadBackedPileup pileup) {
        int[] counts = pileup.getBaseCounts();

        int maxIndex = 0;
        int sumCounts = 0;
        for (int i = 0; i < counts.length; i++) {
            if (counts[i] > counts[maxIndex]) {
                maxIndex = i;
            }
            sumCounts += counts[i];
        }

        // check if homozygous, i.e. total reads = max read
        if (counts[maxIndex] == sumCounts) {
            return BASE.values()[maxIndex];
        }
        return BASE.N;
    }

    /**
     * our map function, which takes the reads spanning this locus, any associated reference ordered data,
     * and the reference information.  We output the log ratio of reads by base as the result of this function
     *
     * @param tracker the reference ordered data tracker
     * @param ref     the reference information
     * @param context the locus context, which contains all of the read information
     * @return an array, {log ratio of reads, number of bases}
     */
    public Double[] map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (ref.getBase() == 'N' || ref.getBase() == 'n') return null; // we don't deal with the N ref base case

        // if (ref.getBase() == 'C' || ref.getBase() == 'G') return new Double[] {0.0, 0.0, 0.0, 0.0};;

        ReadBackedPileup pileup = context.getBasePileup().getBaseAndMappingFilteredPileup(20, 30);

        // ignore non-homozygous sites
        BASE currBase = getHomozygousBase(pileup);

        if (currBase == BASE.N) return null;

        int totalCoverage = pileup.depthOfCoverage();

        if (pileup.depthOfCoverage() >= 10 & pileup.depthOfCoverage() <= 100) {
            int negCoverage = pileup.getNegativeStrandPileup().depthOfCoverage();
            int posCoverage = totalCoverage - negCoverage;

            out.printf("%s\t%s\t%d\t%d\n", context.getLocation(), (char)ref.getBase(), posCoverage, negCoverage);

            switch (currBase) {
                case A:
                    return new Double[]{1.0, (double)totalCoverage, (double)posCoverage, 0.0, 0.0, 0.0, (double)negCoverage, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, Math.log10((posCoverage + 0.5) / (negCoverage + 0.5)), 0.0};
                case T:
                    return new Double[]{1.0, (double)totalCoverage, 0.0, (double)posCoverage, 0.0, 0.0, 0.0, (double)negCoverage, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, Math.log10((posCoverage + 0.5) / (negCoverage + 0.5)), 0.0};
                case C:
                    return new Double[]{1.0, (double)totalCoverage, 0.0, 0.0, (double)posCoverage, 0.0, 0.0, 0.0, (double)negCoverage, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, Math.log10((posCoverage + 0.5) / (negCoverage + 0.5))};
                case G:
                    return new Double[]{1.0, (double)totalCoverage, 0.0, 0.0, 0.0, (double)posCoverage, 0.0, 0.0, 0.0, (double)negCoverage, 0.0, 0.0, 0.0, 1.0, 0.0, Math.log10((posCoverage + 0.5) / (negCoverage + 0.5))};
                case N:
                    break;
            }
        }
        return null;
    }

    /**
        * Provide an initial value for reduce computations. In this case we simply return an empty list
        *
        * @return Initial value of reduce.
        */
    public Double[] reduceInit() {
        return new Double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    }

    /**
     * Outputs the number of bases and sum of log ratio of strandness.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Double[] reduce(Double[] value, Double[] sum) {
        if (value == null) {
            return sum;
        }
        for (int i = 0; i < value.length; i++) {
            sum[i] += value[i];
        }
        return sum;
    }

    /**
     * A composite, 'reduce of reduces' function.
     *
     * @param lhs 'left-most' portion of data in the composite reduce.
     * @param rhs 'right-most' portion of data in the composite reduce.
     * @return The composite reduce type.
     */
    public Double[] treeReduce(Double lhs[], Double rhs[]) {
        Double [] sum = new Double[4];
        for (int i = 0; i < sum.length; i++) {
            sum[i] = lhs[i] + rhs[i];
        }
        return sum;
    }

    /**
     * when we finish traversing, close the result list
     * @param result the final reduce result
     */
    public void onTraversalDone(Double[] result) {
        out.print("TotalCount=");
        for (Double d : result) {
            out.print(d + ",");
        }
        out.println();
    }
}


