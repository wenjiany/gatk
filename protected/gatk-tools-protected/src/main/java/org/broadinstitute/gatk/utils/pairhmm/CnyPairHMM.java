/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.gatk.utils.pairhmm;

import org.broadinstitute.gatk.utils.haplotype.Haplotype;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import htsjdk.variant.variantcontext.Allele;

import java.io.File;
import java.lang.reflect.Field;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public final class CnyPairHMM extends PairHMM implements BatchPairHMM {
    private static class HmmInput {
        public byte[] readBases;
        public byte[] readQuals;
        public byte[] insertionGOP;
        public byte[] deletionGOP;
        public byte[] overallGCP;
        public List<Haplotype> haplotypes;
    }

    public static class ResultQueue {
        private int offset;
        private List<double[]> batchResults;

        public ResultQueue() {
            batchResults = new LinkedList<>();
            offset = 0;
        }

        public void push(double[] results) {
            batchResults.add(results);
        }

        public double pop() {
            double[] results = batchResults.get(0);
            double top = results[offset++];
            if (offset == results.length) {
                batchResults.remove(0);
                offset = 0;
            }
            return top;
        }
    }

    final static String libPath = "/opt/convey/personalities/32100.1.1.1.0";
    final static String libName = "gmvhdl_gatk_hmm";

    private static boolean loaded = false;
    private List<HmmInput> batchRequests = new LinkedList<>();
    private ResultQueue resultQueue = new ResultQueue();

    static public boolean isAvailable() {
        if (!loaded) {
            File library = new File(libPath + "/lib" + libName + ".so");
            return library.exists();
        }
        return true;
    }

    private native void initFpga();
    private native int dequeueRequirement(int reflen, int readlen);
    private native int enqueue(byte[] haplotypeBases,
                               byte[] readBases,
                               byte[] readQuals,
                               byte[] insertionGOP,
                               byte[] deletionGOP,
                               byte[] overallGCP,
                               int hapStartIndex,
                               boolean recacheReadValues);
    private native int flushQueue();
    private native int dequeue(double[] results);
    private native double softHmm(byte[] haplotypeBases,
                                  byte[] readBases,
                                  byte[] readQuals,
                                  byte[] insertionGOP,
                                  byte[] deletionGOP,
                                  byte[] overallGCP,
                                  int hapStartIndex,
                                  boolean recacheReadValues);

    public native void reportStats();

    public void initialize( final int READ_MAX_LENGTH, final int HAPLOTYPE_MAX_LENGTH ) {
        if (!loaded) {
            addLibraryPath(libPath);
            System.loadLibrary(libName);
            initFpga();
            loaded = true;
            System.out.println("FPGA HMM Initialized");
        }
    }

    public void batchAdd(final List<Haplotype> haplotypes,
                         final byte[] readBases,
                         final byte[] readQuals,
                         final byte[] insertionGOP,
                         final byte[] deletionGOP,
                         final byte[] overallGCP) {
        final int numHaplotypes = haplotypes.size();
        HmmInput test = new HmmInput();
        test.readBases = readBases;
        test.readQuals = readQuals;
        test.insertionGOP = insertionGOP;
        test.deletionGOP = deletionGOP;
        test.overallGCP = overallGCP;
        test.haplotypes = haplotypes;
        batchRequests.add(test);
        for (int jjj = 0; jjj < numHaplotypes; jjj++) {
            final boolean recacheReadValues = (jjj == 0);
            final Haplotype haplotype = haplotypes.get(jjj);
            enqueuePrepare(haplotype.getBases(), readBases);
            if (enqueue(haplotype.getBases(), readBases, readQuals, insertionGOP, deletionGOP, overallGCP, 0, recacheReadValues) == 0)
                throw new RuntimeException("FPGA queue overflow in batchAdd");
        }
    }

    public double[] batchGetResult() {
        double[] results;

        int n = flushQueue();
        if (n > 0) {
            results = new double[n];
            if (dequeue(results) != n)
                System.out.println("queue underflow in enqueuePrepare");
            resultQueue.push(results);
        }

        final HmmInput test = batchRequests.remove(0);
        final int numHaplotypes = test.haplotypes.size();
        results = new double[numHaplotypes];
        for (int jjj = 0; jjj < numHaplotypes; jjj++) {
            results[jjj] = resultQueue.pop();
            if (results[jjj]<-60.0) {
                final Haplotype haplotype = test.haplotypes.get(jjj);
                results[jjj]=softHmm(haplotype.getBases(), test.readBases, test.readQuals, test.insertionGOP, test.deletionGOP, test.overallGCP, 0, true);
            }
        }
        return results;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public PerReadAlleleLikelihoodMap computeLikelihoods(final List<GATKSAMRecord> reads, final Map<Allele, Haplotype> alleleHaplotypeMap, final Map<GATKSAMRecord, byte[]> GCPArrayMap){

        // initialize the pairHMM if necessary
        if (! initialized) {
            int readMaxLength = findMaxReadLength(reads);
            int haplotypeMaxLength = findMaxHaplotypeLength(alleleHaplotypeMap);
            initialize(readMaxLength, haplotypeMaxLength);
        }

        // Pass the read bases/quals, and the haplotypes as a list into the HMM
        performBatchAdditions(reads, alleleHaplotypeMap, GCPArrayMap);

        // Get the log10-likelihoods for each read/haplotype ant pack into the results map
        final PerReadAlleleLikelihoodMap likelihoodMap = new PerReadAlleleLikelihoodMap();
        collectLikelihoodResults(reads, alleleHaplotypeMap, likelihoodMap);

        return likelihoodMap;
    }

    private void collectLikelihoodResults(List<GATKSAMRecord> reads, Map<Allele, Haplotype> alleleHaplotypeMap, PerReadAlleleLikelihoodMap likelihoodMap) {
        for(final GATKSAMRecord read : reads){
            final double[] likelihoods = batchGetResult();
            int jjj = 0;
            for (Allele allele : alleleHaplotypeMap.keySet()){
                final double log10l = likelihoods[jjj];
                likelihoodMap.add(read, allele, log10l);
                jjj++;
            }
        }
    }

    private void performBatchAdditions(List<GATKSAMRecord> reads, Map<Allele, Haplotype> alleleHaplotypeMap, Map<GATKSAMRecord, byte[]> GCPArrayMap) {
        final List<Haplotype> haplotypeList = getHaplotypeList(alleleHaplotypeMap);
        for(final GATKSAMRecord read : reads){
            final byte[] readBases = read.getReadBases();
            final byte[] readQuals = read.getBaseQualities();
            final byte[] readInsQuals = read.getBaseInsertionQualities();
            final byte[] readDelQuals = read.getBaseDeletionQualities();
            final byte[] overallGCP = GCPArrayMap.get(read);

            batchAdd(haplotypeList, readBases, readQuals, readInsQuals, readDelQuals, overallGCP);
        }
    }


    protected double subComputeReadLikelihoodGivenHaplotypeLog10( final byte[] haplotypeBases,
                                                                  final byte[] readBases,
                                                                  final byte[] readQuals,
                                                                  final byte[] insertionGOP,
                                                                  final byte[] deletionGOP,
                                                                  final byte[] overallGCP,
                                                                  final int hapStartIndex,
                                                                  final boolean recacheReadValues,
                                                                  final int nextHapStartIndex) {
        return 0.0;
    }

    private List<Haplotype> getHaplotypeList(Map<Allele, Haplotype> alleleHaplotypeMap){
        final List<Haplotype> haplotypeList = new LinkedList<>();
        for (Allele a : alleleHaplotypeMap.keySet()){
             haplotypeList.add(alleleHaplotypeMap.get(a));
        }
        return haplotypeList;
    }

    private void enqueuePrepare(byte[] haplotypeBases, byte[] readBases) {
        double[] results = null;
        int n = dequeueRequirement(haplotypeBases.length, readBases.length);
        if (n>0) {
            results = new double[n];
            if (dequeue(results)!=n)
                System.out.println("queue underflow in enqueuePrepare");
        } else if (n<0) {
            n = flushQueue();
            if (n > 0) {
                results = new double[n];
                if (dequeue(results) != n)
                    System.out.println("queue underflow in enqueuePrepare");
            }
        }

        if (results != null)
            resultQueue.push(results);
    }

    public static void addLibraryPath(String pathToAdd) {
        try {
            final Field usrPathsField = ClassLoader.class.getDeclaredField("usr_paths");
            usrPathsField.setAccessible(true);

            //get array of paths
            final String[] paths = (String[])usrPathsField.get(null);

            //check if the path to add is already present
            for(String path : paths) {
                if(path.equals(pathToAdd)) {
                    return;
                }
            }

            //add the new path
            final String[] newPaths = Arrays.copyOf(paths, paths.length + 1);
            newPaths[newPaths.length-1] = pathToAdd;
            usrPathsField.set(null, newPaths);
        } catch (Exception ex) {
        }
    }
}