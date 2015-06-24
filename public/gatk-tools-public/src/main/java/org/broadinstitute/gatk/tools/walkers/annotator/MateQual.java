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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.*;
import java.io.File;

/**
 * Average Mapping Quality of mate pairs by Alleles
 *
 * <p>This annotation tells you what the average mapping quality of mates. No mapping quality filter is used.</p>
 *
 * <h3>Calculation</h3>
 * $$ MMQ = {sum ref mate MQ/n ref }{sum alt mate MQ/n alt } $$
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZero.php">MappingQualityZero</a></b> gives the count of reads with MAPQ=0 across all samples.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_MappingQualityZeroBySample.php">MappingQualityZeroBySample</a></b> gives the count of reads with MAPQ=0 for each individual sample.</li>
 * </ul>
 */
public class MateQual extends InfoFieldAnnotation {

    // BAMFileReader this.reader = ((VariantAnnotator) walker).secondReader;

    SamReader secondReader;

    int sampleSize = 25;

    @Override
    public void initialize(AnnotatorCompatible walker, GenomeAnalysisEngine toolkit, Set<VCFHeaderLine> headerLines) {
        super.initialize(walker, toolkit, headerLines);
        SAMReaderID id = ((ArrayList<SAMReaderID>)toolkit.getReadsDataSource().getReaderIDs()).get(0);
        String bamfile = id.getSamFile().getAbsolutePath();
        secondReader =  SamReaderFactory.makeDefault().open(new File(bamfile));
    }

    @Override
    protected void finalize() throws Throwable {
        secondReader.close();
        super.finalize();
    }

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        int nRef = 0;
        int nAlt = 0;

        int refMateQual = 0;
        int altMateQual = 0;

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() )
		{
            if (!vc.isSNP()) continue;
            byte ref_base = vc.getAlleles().get(0).getBases()[0];
            byte alt_base = vc.getAlleles().get(1).getBases()[0];

            // too many reads to check..., maybe try sampling later

            if (sample.getValue().getBasePileup().depthOfCoverage() > 200) {
                continue;
            }

            ArrayList<String> dp4 = (ArrayList<String>) vc.getAttribute("DP4");

            double refSample = sampleSize/(1.0 + Integer.parseInt(dp4.get(0)) + Integer.parseInt(dp4.get(1)));
            double altSample = sampleSize/(1.0 + Integer.parseInt(dp4.get(2)) + Integer.parseInt(dp4.get(3)));

            for ( PileupElement p : sample.getValue().getBasePileup() )
			{
                // int x = sample.getValue().getBasePileup().getBaseFilteredPileup(30).getReads().size();

                if (p.getBase()==ref_base && (refSample > 0.9 || Math.random() < refSample)) {
                    nRef++;
                    refMateQual += secondReader.queryMate(p.getRead()).getMappingQuality();
                } else {
                    if (altSample > 0.9 || Math.random() < altSample) {
                        nAlt++;
                        altMateQual += secondReader.queryMate(p.getRead()).getMappingQuality();
                    }
                }
            }
        }

        if (nRef > 0) refMateQual = refMateQual / nRef;
        if (nAlt > 0) altMateQual = altMateQual / nAlt;

        Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%d,%d", refMateQual, altMateQual));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.MMQ); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}
