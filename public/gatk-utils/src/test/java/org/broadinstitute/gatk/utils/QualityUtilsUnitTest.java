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

package org.broadinstitute.gatk.utils;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/21/12
 */

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Basic unit test for QualityUtils class
 */
public class QualityUtilsUnitTest extends BaseTest {
    final private static double TOLERANCE = 1e-9;

    @BeforeClass
    public void init() {
    }

    @DataProvider(name = "QualTest")
    public Object[][] makeMyDataProvider() {
        final List<Object[]> tests = new ArrayList<>();

        for ( int qual = 0; qual < 255; qual++ ) {
            tests.add(new Object[]{(byte)(qual & 0xFF), Math.pow(10.0, ((double)qual)/-10.0)});
        }

        return tests.toArray(new Object[][]{});
    }

    /**
     * Example testng test using MyDataProvider
     */
    @Test(dataProvider = "QualTest")
    public void testMyData(final byte qual, final double errorRate) {
        final double trueRate = 1 - errorRate;

        final double actualErrorRate = QualityUtils.qualToErrorProb(qual);
        Assert.assertEquals(actualErrorRate, errorRate, TOLERANCE);
        final double actualTrueRate = QualityUtils.qualToProb(qual);
        Assert.assertEquals(actualTrueRate, trueRate, TOLERANCE);

        // log10 tests
        final double actualLog10ErrorRate = QualityUtils.qualToErrorProbLog10(qual);
        Assert.assertEquals(actualLog10ErrorRate, Math.log10(errorRate), TOLERANCE);
        final double actualLog10TrueRate = QualityUtils.qualToProbLog10(qual);
        Assert.assertEquals(actualLog10TrueRate, Math.log10(trueRate), TOLERANCE);

        // test that we can convert our error rates to quals, accounting for boundaries
        final int expectedQual = Math.max(Math.min(qual & 0xFF, QualityUtils.MAX_SAM_QUAL_SCORE), 1);
        final byte actualQual = QualityUtils.trueProbToQual(trueRate);
        Assert.assertEquals(actualQual, expectedQual & 0xFF);
        final byte actualQualFromErrorRate = QualityUtils.errorProbToQual(errorRate);
        Assert.assertEquals(actualQualFromErrorRate, expectedQual & 0xFF);

        for ( int maxQual = 10; maxQual < QualityUtils.MAX_SAM_QUAL_SCORE; maxQual++ ) {
            final byte maxAsByte = (byte)(maxQual & 0xFF);
            final byte expectedQual2 = (byte)(Math.max(Math.min(qual & 0xFF, maxQual), 1) & 0xFF);
            final byte actualQual2 = QualityUtils.trueProbToQual(trueRate, maxAsByte);
            Assert.assertEquals(actualQual2, expectedQual2, "Failed with max " + maxQual);
            final byte actualQualFromErrorRate2 = QualityUtils.errorProbToQual(errorRate, maxAsByte);
            Assert.assertEquals(actualQualFromErrorRate2, expectedQual2, "Failed with max " + maxQual);

            // test the integer routines
            final byte actualQualInt2 = QualityUtils.trueProbToQual(trueRate, maxQual);
            Assert.assertEquals(actualQualInt2, expectedQual2, "Failed with max " + maxQual);
            final byte actualQualFromErrorRateInt2 = QualityUtils.errorProbToQual(errorRate, maxQual);
            Assert.assertEquals(actualQualFromErrorRateInt2, expectedQual2, "Failed with max " + maxQual);
        }
    }

    @Test
    public void testTrueProbWithMinDouble() {
        final byte actual = QualityUtils.trueProbToQual(Double.MIN_VALUE);
        Assert.assertEquals(actual, 1, "Failed to convert true prob of min double to 1 qual");
    }

    @Test
    public void testTrueProbWithVerySmallValue() {
        final byte actual = QualityUtils.trueProbToQual(1.7857786272673852E-19);
        Assert.assertEquals(actual, 1, "Failed to convert true prob of very small value 1.7857786272673852E-19 to 1 qual");
    }

    @Test
    public void testQualCaches() {
        Assert.assertEquals(QualityUtils.qualToErrorProb((byte) 20), 0.01, 1e-6);
        Assert.assertEquals(QualityUtils.qualToErrorProbLog10((byte) 20), -2.0, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProb((byte) 20), 0.99, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProbLog10((byte) 20), -0.0043648054, 1e-6);

        Assert.assertEquals(QualityUtils.qualToErrorProb((byte) 30), 0.001, 1e-6);
        Assert.assertEquals(QualityUtils.qualToErrorProbLog10((byte) 30), -3.0, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProb((byte) 30), 0.999, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProbLog10((byte) 30), -0.000434511774, 1e-6);

        Assert.assertEquals(QualityUtils.qualToErrorProb((byte) 40), 0.0001, 1e-6);
        Assert.assertEquals(QualityUtils.qualToErrorProbLog10((byte) 40), -4.0, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProb((byte) 40), 0.9999, 1e-6);
        Assert.assertEquals(QualityUtils.qualToProbLog10((byte) 40), -4.34316198e-5, 1e-6);
    }

    @Test()
    public void testBoundingDefault() {
        for ( int qual = 0; qual < 1000; qual++ ) {
            final byte expected = (byte)Math.max(Math.min(qual, QualityUtils.MAX_SAM_QUAL_SCORE), 1);
            Assert.assertEquals(QualityUtils.boundQual(qual), expected);
        }
    }

    @Test()
    public void testBoundingWithMax() {
        for ( int max = 10; max < 255; max += 50 ) {
            for ( int qual = 0; qual < 1000; qual++ ) {
                final int expected = Math.max(Math.min(qual, max), 1);
                Assert.assertEquals(QualityUtils.boundQual(qual, (byte)(max & 0xFF)) & 0xFF, expected & 0xFF, "qual " + qual + " max " + max);
            }
        }
    }

    @DataProvider(name = "PhredScaleDoubleOps")
    public Object[][] makePhredDoubleTest() {
        final List<Object[]> tests = new ArrayList<>();

        tests.add(new Object[]{0.0, -10 * Math.log10(Double.MIN_VALUE)});
        tests.add(new Object[]{1.0, 0.0});
        for ( int pow = 1; pow < 20; pow++ ) {
            tests.add(new Object[]{Math.pow(10.0, -1.0 * pow), pow * 10});
            tests.add(new Object[]{Math.pow(10.0, -1.5 * pow), pow * 15});
        }

        return tests.toArray(new Object[][]{});
    }

    @Test()
    public void testQualToErrorProbDouble() {
        for ( double qual = 3.0; qual < 255.0; qual += 0.1 ) {
            final double expected = Math.pow(10.0, qual / -10.0);
            Assert.assertEquals(QualityUtils.qualToErrorProb(qual), expected, TOLERANCE, "failed qual->error prob for double qual " + qual);
        }
    }


    @Test(dataProvider = "PhredScaleDoubleOps")
    public void testPhredScaleDoubleOps(final double errorRate, final double expectedPhredScaled) {
        final double actualError = QualityUtils.phredScaleErrorRate(errorRate);
        Assert.assertEquals(actualError, expectedPhredScaled, TOLERANCE);
        final double trueRate = 1 - errorRate;
        final double actualTrue = QualityUtils.phredScaleCorrectRate(trueRate);
        if ( trueRate == 1.0 ) {
            Assert.assertEquals(actualTrue, QualityUtils.MIN_PHRED_SCALED_QUAL);
        } else {
            final double tol = errorRate < 1e-10 ? 10.0 : 1e-3;
            Assert.assertEquals(actualTrue, expectedPhredScaled, tol);
        }
    }
}