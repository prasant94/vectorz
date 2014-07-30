package mikera.matrixx.decompose.impl.eigen;

import mikera.matrixx.Matrix;
import mikera.matrixx.decompose.IEigenResult;
import mikera.vectorz.AVector;
import mikera.vectorz.Vector;
import mikera.vectorz.Vector2;

public class EigenResult implements IEigenResult {
    
    private final AVector[] eigenVectors;
    private final Vector2[] eigenValues;
    
    public EigenResult(Vector2[] eigenValues, Matrix[] eigenVectors) {
        this.eigenValues = eigenValues;
        Matrix[] vectors = eigenVectors;
        this.eigenVectors  = new AVector[vectors.length];
        for(int i=0; i<eigenVectors.length; i++) {
            this.eigenVectors[i] = vectors[i] == null ? null : Vector.create(vectors[i].asDoubleArray());
        }
    }

    public EigenResult(Vector2[] eigenValues) {
        this.eigenValues = eigenValues;
        this.eigenVectors  = null;
    }
    
    public Vector2[] getEigenvalues() {
        return eigenValues;
    }

    public AVector[] getEigenVectors() {
        if (eigenVectors == null)
            throw new UnsupportedOperationException("EigenVectors not computed");
        return eigenVectors;
    }
}
