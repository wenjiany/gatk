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

package org.broadinstitute.gatk.tools.walkers.annotator;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.io.File;

/**
 * Proportion of low quality reads
 *
 * <p>This annotation tells you what fraction of reads have a mapping quality of less than the given threshold of 10 (including 0). Note that certain tools may impose a different minimum mapping quality threshold. For example, HaplotypeCaller excludes reads with MAPQ<20.</p>
 *
 * <h3>Calculation</h3>
 * $$ LowMQ = \frac{# reads with MAPQ=0 + # reads with MAPQ<10}{total # reads} $$
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZero.php">MappingQualityZero</a></b> gives the count of reads with MAPQ=0 across all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZeroBySample.php">MappingQualityZeroBySample</a></b> gives the count of reads with MAPQ=0 for each individual sample.</li>
 * </ul>
 */
public class TrimEndCount extends InfoFieldAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int trim_n_bases = ((VariantAnnotator) walker).trim_n_bases;

        int n_ref_fwd = 0;
        int n_ref_rev = 0;
        int n_alt_fwd = 0;
        int n_alt_rev = 0;

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() )
		{
            if (!vc.isSNP()) continue;
            byte ref_base = vc.getAlleles().get(0).getBases()[0];
            byte alt_base = vc.getAlleles().get(1).getBases()[0];

//            SamReader ss = new SamReader(sample.getValue().getBasePileup().getReads().get(0).mFileSource.getReader());

            for ( PileupElement p : sample.getValue().getBasePileup() )
			{
                int x = sample.getValue().getBasePileup().getBaseFilteredPileup(30).getReads().size();
                int offset = p.getOffset();

                // skip reads with mapping quality less than 30 and base quality less than 20
                if (offset < trim_n_bases || offset > p.getRead().getReadLength() - trim_n_bases || p.getMappingQual() < 30 || p.getRead().getBaseQualities()[offset] < 20) {
                    continue;
                }

                if (p.getRead().getReadBases()[offset]==ref_base) {
                    if (!p.getRead().getReadNegativeStrandFlag()) {
                        n_ref_fwd++;
                    } else {
                        n_ref_rev++;
                    }
                } else {
                    if (!p.getRead().getReadNegativeStrandFlag()) {
                        n_alt_fwd++;
                    } else {
                        n_alt_rev++;
                    }
                }
            }
        }
        Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%d,%d,%d,%d", n_ref_fwd, n_ref_rev, n_alt_fwd, n_alt_rev));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.DP4T); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}
