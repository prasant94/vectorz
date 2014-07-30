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

import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrix;
import mikera.matrixx.algo.Multiplications;
import mikera.matrixx.solve.impl.TriangularSolver;
import mikera.matrixx.solve.impl.lu.LUSolver;
import mikera.vectorz.Vector2;

/**
 * @author Peter Abeles
 */
public class DoubleStepQREigenvector {
    
    public static double EPS = Math.pow(2,-52);

    DoubleStepQREigen implicit;

    // Q matrix from double step QR
    Matrix Q;


    Matrix eigenvectors[];

    Matrix eigenvectorTemp;

    LUSolver solver;

    Vector2 origEigenvalues[];
    int N;

    int splits[];
    int numSplits;

    int x1,x2;

    int indexVal;
    boolean onscript;

    public boolean process( DoubleStepQREigen implicit , AMatrix A , AMatrix Q_h )
    {
        this.implicit = implicit;

        if( N != A.rowCount() ) {
            N = A.rowCount();
            Q = Matrix.create(N,N);
            splits = new int[N];
            origEigenvalues = new Vector2[N];
            eigenvectors = new Matrix[N];
            eigenvectorTemp = Matrix.create(N,1);

            solver = new LUSolver();
        } else {
//            UtilEjml.setnull(eigenvectors);
            eigenvectors = new Matrix[N];
        }
        System.arraycopy(implicit.eigenvalues,0,origEigenvalues,0,N);

        implicit.setup(A);
        implicit.setQ(Q);
        numSplits = 0;
        onscript = true;

//        System.out.println("Orig A");
//        A.print("%12.10f");

        if( !findQandR() )
            return false;

        return extractVectors(Q_h);
    }

    public boolean extractVectors( AMatrix Q_h ) {

//        UtilEjml.memset(eigenvectorTemp.data,0);
        eigenvectorTemp.set(0);
        // extract eigenvectors from the shur matrix
        // start at the top left corner of the matrix
        boolean triangular = true;
        
        for( int i = 0; i < N; i++ ) {

            Vector2 c = implicit.eigenvalues[N-i-1];

            if( triangular && !(c.y==0) )
                triangular = false;

            if( (c.y==0) && eigenvectors[N-i-1] == null) {
                solveEigenvectorDuplicateEigenvalue(c.x,i,triangular);
            }
        }
        
        // translate the eigenvectors into the frame of the original matrix
        if( Q_h != null ) {
            Matrix temp = Matrix.create(N,1);
            for( int i = 0; i < N; i++ ) {
                Matrix v = eigenvectors[i];

                if( v != null ) {
//                    CommonOps.mult(Q_h,v,temp);
                    temp = Multiplications.multiply(Q_h, v);
                    eigenvectors[i] = temp;
                    temp = v;
                }
            }
        }
        return true;
    }

    private void solveEigenvectorDuplicateEigenvalue( double real , int first , boolean isTriangle ) {

        double scale = Math.abs(real);
        if( scale == 0 ) scale = 1;

        eigenvectorTemp = eigenvectorTemp.reshape(N,1);
        eigenvectorTemp.fill(0);

        if( first > 0 ) {
            if( isTriangle ) {
                solveUsingTriangle(real, first , eigenvectorTemp);
            } else {
                solveWithLU(real, first , eigenvectorTemp);
            }
        }

        eigenvectorTemp = eigenvectorTemp.reshape(N,1);

        for( int i = first; i < N; i++ ) {
            Vector2 c = implicit.eigenvalues[N-i-1];

            if( (c.y==0) && Math.abs(c.x-real)/scale < 100.0*EPS ) {
                eigenvectorTemp.data[i] = 1;

                Matrix v = Matrix.create(N,1);
//                CommonOps.multTransA(Q,eigenvectorTemp,v);
                v = Multiplications.multiply(Q.getTranspose(), eigenvectorTemp);
                eigenvectors[N-i-1] = v;
//                NormOps.normalizeF(v);
                v.divide(Math.sqrt(v.elementSquaredSum()));

                eigenvectorTemp.data[i] = 0;
                
            }
        }
        
    }

    private void solveUsingTriangle(double real, int index, Matrix r ) {
        for( int i = 0; i < index; i++ ) {
            implicit.A.addAt(i,i,-real);
        }

//        SpecializedOps.subvector(implicit.A,0,index,index,false,0,r);
        AMatrix sub = implicit.A.subArray(new int[] {0, index}, new int[] {index, 1});
        r.set(sub.reshape(r.rowCount(), r.columnCount()));
//        CommonOps.changeSign(r);
        r.multiply(-1);

        TriangularSolver.solveU(implicit.A.data,r.data,implicit.A.rowCount(),0,index);

        for( int i = 0; i < index; i++ ) {
            implicit.A.addAt(i,i,real);
        }
    }

    private void solveWithLU(double real, int index, Matrix r ) {
//        Matrix A = Matrix.create(index,index);

//        CommonOps.extract(implicit.A,0,index,0,index,A,0,0);
        AMatrix A = implicit.A.subMatrix(0, index, 0, index);

        for( int i = 0; i < index; i++ ) {
            A.addAt(i,i,-real);
        }

        r = r.reshape(index,1);
        eigenvectorTemp = r;

//        SpecializedOps.subvector(implicit.A,0,index,index,false,0,r);
        AMatrix sub = implicit.A.subArray(new int[] {0, index}, new int[] {index, 1});
        r.set(sub.reshape(r.rowCount(), r.columnCount()));
//      CommonOps.changeSign(r);
        r.multiply(-1);

        // TODO this must be very inefficient
        if( solver.setA(A) == null)
            throw new RuntimeException("Solve failed");
        r.set(solver.solveEigenHelper(r, 1e-30));
    }

    public boolean findQandR() {
//        CommonOps.setIdentity(Q);
        Q.setElements(Matrix.createIdentity(Q.rowCount(), Q.columnCount()).data);

        x1 = 0;
        x2 = N-1;

        // use the already computed eigenvalues to recompute the Q and R matrices
        indexVal = 0;
        while( indexVal < N ) {
            if (!findNextEigenvalue()) {
                return false;
            }
        }

//        Q.print("%1.10f");
//
//        implicit.A.print("%1.10f");

        return true;
    }

    private boolean findNextEigenvalue() {
        boolean foundEigen = false;
        while( !foundEigen && implicit.steps < implicit.maxIterations ) {
//            implicit.A.print();
            implicit.incrementSteps();

            if( x2 < x1 ) {
                moveToNextSplit();
            } else if( x2-x1 == 0 ) {
                implicit.addEigenAt(x1);
                x2--;
                indexVal++;
                foundEigen = true;
            } else if( x2-x1 == 1 && !implicit.isReal2x2(x1,x2)) {
                implicit.addComputedEigen2x2(x1,x2);
                x2 -= 2;
                indexVal += 2;
                foundEigen = true;
            } else if( implicit.steps-implicit.lastExceptional > implicit.exceptionalThreshold ) {
//                implicit.A.print("%e");
                //System.err.println("If it needs to do an exceptional shift then something went very bad.");
//                return false;
                implicit.exceptionalShift(x1,x2);
                implicit.lastExceptional = implicit.steps;
            } else if( implicit.isZero(x2,x2-1)) {
                // check for convergence
                implicit.addEigenAt(x2);
                foundEigen = true;
                x2--;
                indexVal++;
            } else {
                checkSplitPerformImplicit();
            }
        }
        return foundEigen;
    }


    private void checkSplitPerformImplicit() {
        // check for splits
        for( int i = x2; i > x1; i-- ) {
            if( implicit.isZero(i,i-1)) {
                x1 = i;
                splits[numSplits++] = i-1;
                // reduce the scope of what it is looking at
                return;
            }
        }
        // first try using known eigenvalues in the same order they were originally found
        if( onscript) {
            if( implicit.steps > implicit.exceptionalThreshold/2  ) {
                onscript = false;
            } else {
                Vector2 a = origEigenvalues[indexVal];

                // if no splits are found perform an implicit step
                if( a.y==0 ) {
                    implicit.performImplicitSingleStep(x1,x2, a.x);
                } else if( x2 < N-2 ) {
                    implicit.performImplicitDoubleStep(x1,x2, a.x,a.y);
                } else {
                    onscript = false;
                }
            }
        } else {
            // that didn't work so try a modified order
            if( x2-x1 >= 1 && x2 < N-2 )
                implicit.implicitDoubleStep(x1,x2);
            else
                implicit.performImplicitSingleStep(x1,x2,implicit.A.get(x2,x2));
        }
    }


    private void moveToNextSplit() {
        if( numSplits <= 0 )
            throw new RuntimeException("bad");

        x2 = splits[--numSplits];

        if( numSplits > 0 ) {
            x1 = splits[numSplits-1]+1;
        } else {
            x1 = 0;
        }
    }

    public Matrix getQ() {
        return Q;
    }

    public DoubleStepQREigen getImplicit() {
        return implicit;
    }

    public Matrix[] getEigenvectors() {
        return eigenvectors;
    }

    public Vector2[] getEigenvalues() {
        return implicit.eigenvalues;
    }
}
