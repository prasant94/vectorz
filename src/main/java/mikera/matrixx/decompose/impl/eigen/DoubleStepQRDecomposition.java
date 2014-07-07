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
import mikera.matrixx.decompose.impl.hessenberg.HessenbergResult;
import mikera.matrixx.decompose.impl.hessenberg.HessenbergSimilarDecomposition;
import mikera.vectorz.AVector;
import mikera.vectorz.Vector2;

/**
 * <p>
 * Finds the eigenvalue decomposition of an arbitrary square matrix using the implicit double-step QR algorithm.
 * Watched is included in its name because it is designed to print out internal debugging information.  This
 * class is still underdevelopment and has yet to be optimized.
 * </p>
 *
 * <p>
 * Based off the description found in:<br>
 * David S. Watkins, "Fundamentals of Matrix Computations." Second Edition.
 * </p>
 *
 * @author Peter Abeles
 */
//TODO looks like there might be some pointless copying of arrays going on
public class DoubleStepQRDecomposition {

    HessenbergResult hessenbergResult;
    DoubleStepQREigenvalue algValue;
    DoubleStepQREigenvector algVector;

    AMatrix H;

    // should it compute eigenvectors or just eigenvalues
    boolean computeVectors;

    public DoubleStepQRDecomposition( boolean computeVectors ) {
        algValue = new DoubleStepQREigenvalue();
        algVector = new DoubleStepQREigenvector();

        this.computeVectors = computeVectors;
    }

    public boolean decompose(AMatrix A) {

        hessenbergResult = HessenbergSimilarDecomposition.decompose(A);
        if( hessenbergResult == null )
            return false;

        H = hessenbergResult.getH();

        algValue.getImplicitQR().createR = false;
//        algValue.getImplicitQR().setChecks(true,true,true);

        if( !algValue.process(H) )
            return false;
        
//        for( int i = 0; i < A.numRows; i++ ) {
//            System.out.println(algValue.getEigenvalues()[i]);
//        }

        algValue.getImplicitQR().createR = true;

        if( computeVectors )
            return algVector.process(algValue.getImplicitQR(), H, hessenbergResult.getQ());
        else
            return true;
    }

    public int getNumberOfEigenvalues() {
        return algValue.getEigenvalues().length;
    }

    public Vector2 getEigenvalue(int index) {
        return algValue.getEigenvalues()[index];
    }

    public AVector getEigenVector(int index) {
        return algVector.getEigenvectors()[index];
    }
}
