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

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;

import htsjdk.variant.vcf.VCFInfoHeaderLine;

import org.broadinstitute.gatk.utils.pileup.PileupElement;

import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Number of identical nucleutides as REF
 *
 * <p>Number of identical nucleutides as REF</p>
 *
 * <h3>Calculation</h3>
 * $$ PolyX = Number of identical nucleutides $$
 *
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b></b></li>
 *     <li><b></b></li>
 * </ul>
 */
public class PolyX extends InfoFieldAnnotation {

    private int countRepeat(int offSet, byte[] bases) {
        byte base = bases[offSet];

        int forward = 0;
        for(int i = offSet; i<bases.length;i++) {
            if (base == bases[i]) {
                forward++;
            } else {
                break;
            }
        }

        int backward = 0;
        for(int i = offSet; i>0; i--) {
            if (base == bases[i]) {
                backward++;
            } else {
                break;
            }
        }
        return forward + backward - 1;
    };

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        byte[] bases = ref.getBases();

        int currPos = (bases.length - 1)/2;

        int currPoly = countRepeat(currPos, bases);

        // Check flanking repeats;
        // Give less weight to flanking repeats, ie. adjacent repeats have to be two bases longer.

        if (bases[currPos-1] != bases[currPos]) {
            int prevPoly = countRepeat(currPos - 1, bases);
            if (currPoly < prevPoly + 2) currPoly = prevPoly;
        }

        if (bases[currPos+1] != bases[currPos]) {
            int nextPoly = countRepeat(currPos + 1, bases);
            if (currPoly < nextPoly + 2) currPoly = nextPoly;
        }

        Map<String, Object> map = new HashMap<>();
        map.put(getKeyNames().get(0), String.format("%d", currPoly));
        return map;
    }

    public List<String> getKeyNames() { return Arrays.asList(GATKVCFConstants.PolyX); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(GATKVCFHeaderLines.getInfoLine(getKeyNames().get(0))); }
}
