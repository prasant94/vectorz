package mikera.matrixx.decompose;

import mikera.matrixx.AMatrix;
import mikera.matrixx.decompose.impl.eigen.DoubleStepQRDecomposition;

public class Eigen {
    public static IEigenResult decompose(AMatrix A, boolean computeVectors) {
        return DoubleStepQRDecomposition.decompose(A, computeVectors);
    }
    public static IEigenResult decompose(AMatrix A) {
        return decompose(A, true);
    }
}
