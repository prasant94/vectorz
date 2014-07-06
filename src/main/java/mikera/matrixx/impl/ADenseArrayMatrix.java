package mikera.matrixx.impl;

import mikera.arrayz.INDArray;
import mikera.arrayz.impl.IDenseArray;
import mikera.matrixx.AMatrix;
import mikera.vectorz.AVector;
import mikera.vectorz.Vector;
import mikera.vectorz.Vectorz;
import mikera.vectorz.impl.ADenseArrayVector;
import mikera.vectorz.impl.AStridedVector;
import mikera.vectorz.util.DoubleArrays;
import mikera.vectorz.util.ErrorMessages;

/**
 * Abstract base class for matrices wrapping a dense (rows*cols) subset of a double[] array
 * @author Mike
 *
 */
public abstract class ADenseArrayMatrix extends AStridedMatrix implements IFastRows, IDenseArray {
	private static final long serialVersionUID = -2144964424833585026L;

	protected ADenseArrayMatrix(double[] data, int rows, int cols) {
		super(data, rows, cols);
	}

	@Override
	public abstract int getArrayOffset();
	
	@Override
	public boolean isPackedArray() {
		return (getArrayOffset()==0) && (data.length ==(rows*cols));
	}
	
	@Override
	public boolean isZero() {
		return DoubleArrays.isZero(data, getArrayOffset(), rows*cols);
	}
	
	@Override
	public boolean isUpperTriangular() {
		// triangular test, taking into account cache layout to access via rows
		int rc=rowCount();
		int cc=columnCount();
		int offset=getArrayOffset();
		for (int i=1; i<rc; i++) {
			if (!DoubleArrays.isZero(data, offset+i*cc, Math.min(cc, i))) return false;
		}
		return true;
	}
	
	@Override
	public boolean isLowerTriangular() {
		// triangular test, taking into account cache layout to access via rows
		int offset=getArrayOffset();
		int cc=columnCount();
		int testRows=Math.min(cc, rowCount());
		for (int i=0; i<testRows; i++) {
			if (!DoubleArrays.isZero(data, offset+i+1, cc-i-1)) return false;
			offset+=cc;
		}
		return true;
	}
	
	@Override
	public int rowStride() {
		return cols;
	}
	
	@Override
	public int columnStride() {
		return 1;
	}
	
	@Override
	public double unsafeGet(int i, int j) {
		return data[index(i,j)];
	}
	
	@Override
	public void set(AVector v) {
		int rc=rowCount();
		if (v.length()!=cols) throw new IllegalArgumentException(ErrorMessages.incompatibleShapes(this, v));
		double[] data=getArray();
		int offset=getArrayOffset();
		for (int i=0; i<rc; i++) {
			v.getElements(data, offset+i*cols);
		}
	}
	
	@Override
	public void set(AMatrix m) {
		if (!isSameShape(m)) throw new IllegalArgumentException(ErrorMessages.incompatibleShapes(this, m));
		double[] data=getArray();
		int offset=getArrayOffset();
		m.getElements(data, offset);
	}
	
	@Override
	public ADenseArrayVector getRowView(int i) {
		return Vectorz.wrap(data, getArrayOffset()+i*cols, cols);
	}
	
	@Override
	public AStridedVector getColumnView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i, rows, cols);
	}
	
	@Override
	public double elementSum() {
		return DoubleArrays.elementSum(data,getArrayOffset(), rows*cols);
	}
	
	@Override
	public double elementSquaredSum() {
		return DoubleArrays.elementSquaredSum(data,getArrayOffset(), rows*cols);
	}
	
	@Override
	public double elementMax(){
		return DoubleArrays.elementMax(data,getArrayOffset(), rows*cols);
	}
	
	@Override
	public double elementMin(){
		return DoubleArrays.elementMin(data,getArrayOffset(), rows*cols);
	}
	
	@Override
	public void copyRowTo(int row, double[] dest, int destOffset) {
		System.arraycopy(data, getArrayOffset()+row*cols, dest, destOffset, cols);
	}
	
	@Override
	public void unsafeSet(int i, int j,double value) {
		data[index(i,j)]=value;
	}
	
	protected int index(int row, int col) {
		return getArrayOffset()+(row*cols)+col;
	}
	
	@Override
	public void transform(AVector source, AVector dest) {
		if ((source instanceof Vector )&&(dest instanceof Vector)) {
			transform ((Vector)source, (Vector)dest);
			return;
		}
		if(rows!=dest.length()) throw new IllegalArgumentException(ErrorMessages.wrongDestLength(dest));
		if(cols!=source.length()) throw new IllegalArgumentException(ErrorMessages.wrongSourceLength(source));
		double[] data=getArray();
		int offset=getArrayOffset();
		for (int i=0; i<rows; i++) {
			dest.unsafeSet(i,source.dotProduct(data, offset+ i*cols));
		}
	}
	
	@Override
	public void add(AVector v) {
		int rc=rowCount();
		int cc=columnCount();
		if(cc!=v.length()) throw new IllegalArgumentException(ErrorMessages.mismatch(this, v));
		double[] data=getArray();
		int offset=getArrayOffset();

		for (int i=0; i<rc; i++) {
			v.addToArray(data, offset+i*cc);
		}		
	}
	
	@Override
	public void add(AMatrix a) {
		if (!isSameShape(a)) {
			throw new IllegalArgumentException(ErrorMessages.incompatibleShapes(this, a));
		}
		a.addToArray(getArray(), getArrayOffset());
	}
	
	@Override
	public void addToArray(double[] data, int offset) {
		DoubleArrays.add(getArray(), getArrayOffset(), data, offset, rows*cols);
	}
	
	@Override
	public boolean equals(AMatrix a) {
		if (!isSameShape(a)) return false;
		return a.equalsArray(getArray(), getArrayOffset());
	}
	
	@Override
	public boolean equals(INDArray a) {
		if (!isSameShape(a)) return false;
		return a.equalsArray(getArray(), getArrayOffset());
	}
	
	@Override
	public boolean equalsArray(double[] data, int offset) {
		return DoubleArrays.equals(getArray(), getArrayOffset(), data, offset, rows*cols);
	}
	
	@Override
	public boolean equals(ADenseArrayMatrix m) {
		if (!isSameShape(m)) return false;
		
		return DoubleArrays.equals(getArray(), getArrayOffset(), m.getArray(), m.getArrayOffset(), rows*cols);
	}

}
