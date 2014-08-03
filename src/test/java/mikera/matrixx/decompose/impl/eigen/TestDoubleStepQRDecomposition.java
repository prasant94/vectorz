/*
 * Copyright (c) 2009-2013, Peter Abeles. All Rights Reserved.
 *
 * This file is part of Efficient Java Matrix Library (EJML).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package mikera.matrixx.decompose.impl.eigen;

import org.junit.Test;

import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrix;
import mikera.matrixx.Matrixx;
import mikera.matrixx.algo.Multiplications;
import mikera.matrixx.decompose.LUP;
import mikera.matrixx.decompose.impl.lu.AltLU;
import mikera.vectorz.Vector2;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;


/**
 * @author Peter Abeles
 */
public class TestDoubleStepQRDecomposition {

    boolean computeVectors;
    public static double EPS = Math.pow(2,-52);

    private DoubleStepQRDecomposition createDecomposition()
    {
        return new DoubleStepQRDecomposition(computeVectors);
    }
    
    @Test
    public void allTests() {
        computeVectors = true;

        checkRandom();
        checkKnownReal();
        checkKnownComplex();
        checkRandomSymmetric();
        checkExceptional();
        checkIdentity();
        checkAllZeros();
        checkLargeValue(false);
        checkLargeValue(true);
    }

    /**
     * Tests for when it just computes eigenvalues
     */
    @Test
    public void justEigenValues() {
        computeVectors = false;

        checkKnownReal_JustValue();
        checkKnownSymmetric_JustValue();
    }

    /**
     * Create a variety of different random matrices of different sizes and sees if they pass the standard
     * eigen decompositions tests.
     */
    public void checkRandom() {
        int sizes[] = new int[]{1,2,5,10,20,50,100,200};

        DoubleStepQRDecomposition alg = createDecomposition();

        for( int s = 2; s < sizes.length; s++ ) {
            int N = sizes[s];
//            System.out.println("N = "+N);

            for( int i = 0; i < 2; i++ ) {
                Matrix A = Matrix.createRandom(N,N);
                A.multiply(2);
                A.sub(1);
                
                assertNotNull(alg._decompose(A));

                performStandardTests(alg,A,-1);
            }
        }
    }

    /**
     * Compare results against a simple matrix with known results where all the eigenvalues
     * are real.  Octave was used to test the known values.
     */
    public void checkKnownReal() {
        Matrix A = Matrix.create(new double[][] {{0.907265, 0.832472, 0.255310}, 
                                                 {0.667810, 0.871323, 0.612657}, 
                                                 {0.025059, 0.126475, 0.427002}});

        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(A));
        performStandardTests(alg,A,-1);

        testForEigenpair(alg,1.686542,0,-0.739990,-0.667630,-0.081761);
        testForEigenpair(alg,0.079014,0,-0.658665,0.721163,-0.214673);
        testForEigenpair(alg,0.440034,0,-0.731422,0.211711,0.648229);
    }

    /**
     * Sees if it correctly computed the eigenvalues.  Does not check eigenvectors.
     */
    public void checkKnownReal_JustValue() {
        Matrix A = Matrix.create(new double[][] {{0.907265, 0.832472, 0.255310}, 
                                                 {0.667810, 0.871323, 0.612657}, 
                                                 {0.025059, 0.126475, 0.427002}});

        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(A));

        testForEigenvalue(alg,A,1.686542,0,1);
        testForEigenvalue(alg,A,0.079014,0,1);
        testForEigenvalue(alg,A,0.440034,0,1);
    }

    /**
     * Sees if it correctly computed the eigenvalues.  Does not check eigenvectors.
     */
    public void checkKnownSymmetric_JustValue() {
        Matrix A = Matrix.create(new double[][]
               {{0.98139,   0.78650,   0.78564},
                {0.78650,   1.03207,   0.29794},
                {0.78564,   0.29794,   0.91926}});
        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(A));

        testForEigenvalue(alg,A,0.00426,0,1);
        testForEigenvalue(alg,A,0.67856,0,1);
        testForEigenvalue(alg,A,2.24989,0,1);
    }

    /**
     * Compare results against a simple matrix with known results where some the eigenvalues
     * are real and some are complex.
     */
    public void checkKnownComplex() {
        Matrix A = Matrix.create(new double[][] {{-0.418284, 0.279875, 0.452912}, 
                                                 {-0.093748, -0.045179, 0.310949}, 
                                                 {0.250513, -0.304077, -0.031414}});

        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(A));
        performStandardTests(alg,A,-1);

        testForEigenpair(alg,-0.39996,0,0.87010,0.43425,-0.23314);
        testForEigenpair(alg,-0.04746,0.02391);
        testForEigenpair(alg,-0.04746,-0.02391);
    }

    /**
     * Check results against symmetric matrices that are randomly generated
     */
    public void checkRandomSymmetric() {
        for( int N = 1; N <= 15; N++ ) {
            for( int i = 0; i < 20; i++ ) {
//                DenseMatrix64F A = RandomMatrices.createSymmetric(N,-1,1,rand);
                AMatrix z = Matrixx.createRandomMatrix(3, 3);
                AMatrix A = z.innerProduct(z.getTranspose());

                DoubleStepQRDecomposition alg = createDecomposition();

                assertNotNull(alg._decompose(A));
            }
        }
    }

    /**
     * For some eigenvector algorithms this is a difficult matrix that requires a special
     * check for.  If it fails that check it will either loop forever or exit before converging.
     */
    public void checkExceptional() {
        Matrix A = Matrix.create(new double [][] 
                {{0, 0, 0, 0, 1}, 
                {1, 0, 0, 0, 0}, 
                {0, 1, 0, 0, 0}, 
                {0, 0, 1, 0, 0}, 
                {0, 0, 0, 1, 0}});

        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(A));
        
        performStandardTests(alg,A,1);
    }

    public void checkIdentity() {
        Matrix I = Matrix.createIdentity(4);

        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(I));

        performStandardTests(alg,I,4);

        testForEigenpair(alg,1,0,1,0,0,0);
        testForEigenpair(alg,1,0,0,1,0,0);
        testForEigenpair(alg,1,0,0,0,1,0);
        testForEigenpair(alg,1,0,0,0,0,1);
    }

    public void checkAllZeros() {
        Matrix A = Matrix.create(5,5);

        DoubleStepQRDecomposition alg = createDecomposition();

        assertNotNull(alg._decompose(A));

        performStandardTests(alg,A,5);
        testEigenvalues(alg,0);
    }

//    public void checkWithSomeRepeatedValuesSymm() {
//        DoubleStepQRDecomposition alg = createDecomposition();
//
//        checkSymmetricMatrix(alg,2,-3,-3,-3);
//        checkSymmetricMatrix(alg,2,-3,2,2);
//        checkSymmetricMatrix(alg,1,1,1,2);
//    }
//
//    public void checkWithSingularSymm() {
//
//        DoubleStepQRDecomposition alg = createDecomposition();
//
//        checkSymmetricMatrix(alg,1,0,1,2);
//    }
//
//    /**
//     * Creates a random symmetric matrix with the specified eigenvalues.  It then
//     * checks to see if it has the expected results.
//     */
//    private void checkSymmetricMatrix(DoubleStepQRDecomposition alg , double ...ev ) {
//        int numRepeated[] = new int[ ev.length ];
//
//        for( int i = 0; i < ev.length ; i++ ) {
//            int num = 0;
//
//            for (double anEv : ev) {
//                if (ev[i] == anEv)
//                    num++;
//            }
//            numRepeated[i] = num;
//        }
//
//        for( int i = 0; i < 200; i++ ) {
//            Matrix A = RandomMatrices.createEigenvaluesSymm(ev.length,rand,ev);
//
//            assertTrue(safeDecomposition(alg,A));
//
//            performStandardTests(alg,A,ev.length);
//
//            for( int j = 0; j < ev.length; j++ ) {
//                testForEigenvalue(alg,A,ev[j],0,numRepeated[j]);
//            }
//        }
//    }

    public void checkLargeValue( boolean symmetric) {

        DoubleStepQRDecomposition alg = createDecomposition();

        for( int i = 0; i < 20; i++ ) {
            Matrix z = Matrix.createRandom(3, 3);
            Matrix y = z.innerProduct(z.getTranspose());  // symmetric
            Matrix A = symmetric ?
                    y :
                    Matrix.createRandom(4,4);

//            CommonOps.scale(1e-200,A);
            A.scale(1e100);

            assertNotNull(alg._decompose(A));

            performStandardTests(alg,A,-1);
        }
    }

    /**
     * If the eigenvalues are all known, real, and the same this can be used to check them.
     */
    public void testEigenvalues( DoubleStepQRDecomposition alg , double expected ) {

        for( int i = 0; i < alg.getNumberOfEigenvalues(); i++ ) {
            Vector2 c = alg.getEigenvalue(i);

            assertTrue(c.y==0);

            assertEquals(expected,c.x,1e-8);
        }
    }

    /**
     * Preforms standard tests that can be performed on any decomposition without prior knowledge of
     * what the results should be.
     */
    public void performStandardTests( DoubleStepQRDecomposition alg , Matrix A , int numReal )
    {

        // basic sanity tests
        assertEquals(A.rowCount(),alg.getNumberOfEigenvalues());

        if( numReal >= 0 ) {
            for( int i = 0; i < A.rowCount(); i++ ) {
                Vector2 v = alg.getEigenvalue(i);

                assertFalse( Double.isNaN(v.x ));
                if( v.y==0 )
                    numReal--;
                else if( Math.abs(v.y) < 10*EPS)
                    numReal--;
            }

            // if there are more than the expected number of real eigenvalues this will
            // be negative
            assertEquals(0,numReal);
        }

//        checkCharacteristicEquation(alg,A);
        if( computeVectors ) {
            testPairsConsistent(alg,A);
            testVectorsLinearlyIndependent(alg);
        }
    }

    /**
     * Checks to see if an eigenvalue is complex then the eigenvector is null.  If it is real it
     * then checks to see if the equation A*v = lambda*v holds true.
     */
    public void testPairsConsistent( DoubleStepQRDecomposition alg , Matrix A )
    {
//        System.out.println("-------------------------------------------------------------------------");
        int N = alg.getNumberOfEigenvalues();

        for( int i = 0; i < N; i++ ) {
            Vector2 c = alg.getEigenvalue(i);
            AMatrix v = alg.getEigenVector(i);

            if( Double.isInfinite(c.x) || Double.isNaN(c.x) ||
                    Double.isInfinite(c.y) || Double.isNaN(c.y))
                fail("Uncountable eigenvalue");

            if( Math.abs(c.y) > 1e-20 ) {
                assertTrue(v==null);
            } else {
                assertTrue(v != null );
//                if( MatrixFeatures.hasUncountable(v)) {
//                    throw new RuntimeException("Egads");
//                }
                assertFalse(v.hasUncountable());

//                CommonOps.mult(A,v,tempA);
                AMatrix tempA = Multiplications.multiply(A, v);
//                CommonOps.scale(c.real,v,tempB);
                AMatrix tempB = v.multiplyCopy(c.x);

//                double max = NormOps.normPInf(A);
                double max = normPInf(A);
                if( max == 0 ) max = 1;
                
                double error = diffNormF(tempA,tempB)/max;
                
                if( error > 1e-12 ) {
//                    System.out.println("Original matrix:");
//                    System.out.println(A);
//                    A.print();
//                    System.out.println("Eigenvalue = "+c.x);
//                    Eigenpair p = EigenOps.computeEigenVector(A,c.real);
//                    p.vector.print();
//                    v.print();
//
//
//                    CommonOps.mult(A,p.vector,tempA);
//                    CommonOps.scale(c.real,p.vector,tempB);
//
//                    max = NormOps.normPInf(A);
//
//                    System.out.println("error before = "+error);
//                    error = SpecializedOps.diffNormF(tempA,tempB)/max;
//                    System.out.println("error after = "+error);
//                    A.print("%f");
//                    System.out.println();
                    fail("Error was too large");
                }

                assertTrue(error <= 1e-12);
            }
        }
    }
    
    private static double normPInf( AMatrix A ) {
        return A.absCopy().elementMax();
    }

    private double diffNormF(AMatrix tempA, AMatrix tempB)
    {
        AMatrix temp = tempA.copy();
        temp.sub(tempB);
        double total = temp.elementSquaredSum();
        temp.abs();
        double scale = temp.elementMax();
        return Math.abs(scale-0) > 1e-12 ? total/scale : 0;
    }

    private double inducedPInf(Matrix A)
    {
        double max = 0;

        int m = A.rowCount();
        int n = A.columnCount();

        for( int i = 0; i < m; i++ ) {
            double total = 0;
            for( int j = 0; j < n; j++ ) {
                total += Math.abs(A.get(i,j));
            }
            if( total > max ) {
                max = total;
            }
        }

        return max;
    }

    /**
     * See if eigenvalues cause the characteristic equation to have a value of zero
     */
    public void checkCharacteristicEquation( DoubleStepQRDecomposition alg ,
                                             Matrix A ) {
        int N = alg.getNumberOfEigenvalues();

        Matrix a = Matrix.create(A);

        for( int i = 0; i < N; i++ ) {
            Vector2 c = alg.getEigenvalue(i);

            if( Math.abs(c.y - 0) < 1e-8 ) {
                // test using the characteristic equation
                Matrix temp = Matrix.createIdentity(A.columnCount());
                temp.scale(c.x);
                temp.sub(a);
                double det = temp.determinant();

                // extremely crude test.  given perfect data this is probably considered a failure...  However,
                // its hard to tell what a good test value actually is.
                assertEquals(0, det, 0.1);
            }
        }
    }

    /**
     * Checks to see if all the real eigenvectors are linearly independent of each other.
     */
    public void testVectorsLinearlyIndependent( DoubleStepQRDecomposition alg ) {
        int N = alg.getNumberOfEigenvalues();

        // create a matrix out of the eigenvectors
        Matrix A = Matrix.create(N,N);

        int off = 0;
        for( int i = 0; i < N; i++ ) {
            AMatrix v = alg.getEigenVector(i);

            // it can only handle real eigenvectors
            if( v == null )
                off++;
            else {
                for( int j = 0; j < N; j++ ) {
                    A.set(i-off,j,v.get(j,0));
                }
            }
        }

        // see if there are any real eigenvectors
        if( N == off )
            return;

//        A.reshape(N-off,N, false);
        A = A.reshape(N-off, N);
        
        AltLU lu = new AltLU();
        lu._decompose(A);
        assertFalse(lu.isSingular());
//        assertTrue(MatrixFeatures.isRowsLinearIndependent(A));
    }

    /**
     * Sees if the pair of eigenvalue and eigenvector was found in the decomposition.
     */
    public void testForEigenpair( DoubleStepQRDecomposition alg , double valueReal ,
                                  double valueImg , double... vector )
    {
        int N = alg.getNumberOfEigenvalues();

        int numMatched = 0;
        for( int i = 0; i < N; i++ ) {
            Vector2 c = alg.getEigenvalue(i);
            
            if( Math.abs(c.x-valueReal) < 1e-4 && Math.abs(c.y-valueImg) < 1e-4) {

//                if( c.isReal() ) {
                if(Math.abs(c.y - 0) < 1e-8)
                    if( vector.length > 0 ) {
                        AMatrix v = alg.getEigenVector(i);
                        AMatrix e = Matrix.create(vector);
                        e = e.getTranspose();

                        double error = diffNormF(e,v);
//                        CommonOps.changeSign(e);
                        e.multiply(-1);
                        double error2 = diffNormF(e,v);


                        if(error < 1e-3 || error2 < 1e-3)
                            numMatched++;
                    } else {
                        numMatched++;
                } else if( Math.abs(c.y-0) > 1e-8 ) {
                    numMatched++;
                }
            }
        }

        assertEquals(1,numMatched);
    }

    public void testForEigenvalue( DoubleStepQRDecomposition alg ,
                                   Matrix A,
                                   double valueReal ,
                                   double valueImg , int numMatched ) {
        int N = alg.getNumberOfEigenvalues();

        int numFound = 0;
        for( int i = 0; i < N; i++ ) {
            Vector2 c = alg.getEigenvalue(i);

            if( Math.abs(c.x-valueReal) < 1e-4 && Math.abs(c.y-valueImg) < 1e-4) {
                numFound++;
            }
        }

        assertEquals(numMatched,numFound);
    }
}
