package mikera.matrixx.impl;

import java.util.Iterator;

import mikera.arrayz.impl.IStridedArray;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrixx;
import mikera.vectorz.AVector;
import mikera.vectorz.Op;
import mikera.vectorz.Vectorz;
import mikera.vectorz.impl.AStridedVector;
import mikera.vectorz.util.ErrorMessages;

/**
 * Abstract base class for arbitrary strided matrices
 * 
 * @author Mike
 */
public abstract class AStridedMatrix extends AArrayMatrix implements IStridedArray {
	private static final long serialVersionUID = -8908577438753599161L;

	protected AStridedMatrix(double[] data, int rows, int cols) {
		super(data, rows, cols);
	}

	public abstract int getArrayOffset();

	public abstract int rowStride();
	
	public abstract int columnStride();	
	
	@Override
	public AStridedMatrix subMatrix(int rowStart, int rowCount, int colStart, int colCount) {
		if ((rowStart<0)||(rowStart>=this.rows)||(colStart<0)||(colStart>=this.cols)) throw new IndexOutOfBoundsException(ErrorMessages.position(rowStart,colStart));
		if ((rowStart+rowCount>this.rows)||(colStart+colCount>this.cols)) throw new IndexOutOfBoundsException(ErrorMessages.position(rowStart+rowCount,colStart+colCount));
		int rowStride=rowStride();
		int colStride=columnStride();
		int offset=getArrayOffset();
		return StridedMatrix.wrap(data, rowCount, colCount, offset+rowStart*rowStride+colStart*colStride, rowStride, colStride);
	}
	
	@Override
	public AStridedVector getRowView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i*rowStride(), cols, columnStride());
	}
	
	@Override
	public double diagonalProduct() {
		int n=Math.min(rowCount(), columnCount());
		int offset=getArrayOffset();
		int st=rowStride()+columnStride();
		double[] data=getArray();
		double result=1.0;
		for (int i=0; i<n; i++) {
			result*=data[offset];
			offset+=st;
		}
		return result;
	}
	
	@Override
	public double trace() {
		int n=Math.min(rowCount(), columnCount());
		int offset=getArrayOffset();
		int st=rowStride()+columnStride();
		double[] data=getArray();
		double result=0.0;
		for (int i=0; i<n; i++) {
			result+=data[offset];
			offset+=st;
		}
		return result;
	}
	
	@Override
	public AStridedVector getColumnView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i*columnStride(), rows, rowStride());
	}
	
	@Override
	public AStridedVector getBand(int i) {
		int cs=columnStride();
		int rs=rowStride();
		if ((i>cols)||(i<-rows)) throw new IndexOutOfBoundsException(ErrorMessages.invalidBand(this, i));
		return Vectorz.wrapStrided(data, getArrayOffset()+bandStartColumn(i)*cs+bandStartRow(i)*rs, bandLength(i), rs+cs);
	}
	
	@Override
	public void add(AVector v) {
		checkColumnCount(v.length());
		int offset=getArrayOffset();
		int colStride=columnStride();
		int rowStride=rowStride();
		for (int i=0; i<rows; i++) {
			v.addToArray(data, offset+i*rowStride, colStride);
		}
	}
	
	@Override
	public void addToArray(double[] dest, int destOffset) {
		int offset=getArrayOffset();
		int colStride=columnStride();
		int rowStride=rowStride();
		for (int i=0; i<rows; i++) {
			int ro=offset+i*rowStride;
			for (int j=0; j<cols; j++) {
				dest[destOffset++]+=data[ro+j*colStride];
			}
		}
	}
	
	@Override
	public void applyOp(Op op) {
		int offset=getArrayOffset();
		int colStride=columnStride();
		int rowStride=rowStride();
		for (int i=0; i<rows; i++) {
			int ro=offset+i*rowStride;
			for (int j=0; j<cols; j++) {
				int ix=ro+j*colStride;
				data[ix]=op.apply(data[ix]);
			}
		}
	}
	
	@Override
	public void add(AMatrix m) {
		checkSameShape(m);
		int offset=getArrayOffset();
		int colStride=columnStride();
		int rowStride=rowStride();
		for (int i=0; i<rows; i++) {
			m.getRow(i).addToArray(data, offset+i*rowStride, colStride);
		}
	}
	
	@Override
	public abstract void copyRowTo(int row, double[] dest, int destOffset);
	
	@Override
	public abstract void copyColumnTo(int col, double[] dest, int destOffset);
	
	@Override
	public int[] getStrides() {
		return new int[] {rowStride(), columnStride()};
	}
	
	@Override
	public int getStride(int dimension) {
		switch (dimension) {
			case 0: return rowStride();
			case 1: return columnStride();
			default: throw new IllegalArgumentException(ErrorMessages.invalidDimension(this, dimension));
		}
	}
	
	@Override
	public Iterator<Double> elementIterator() {
		return new StridedElementIterator(this);
	}
	
	
	@Override
	public AMatrix getTranspose() {
		return getTransposeView();
	}
	
	@Override
	public AMatrix getTransposeView() {
		return Matrixx.wrapStrided(getArray(),columnCount(),rowCount(),getArrayOffset(),columnStride(),rowStride());
	}
	
	@Override
	public boolean isPackedArray() {
		return (getArrayOffset()==0)&&(columnStride()==1)&&(rowStride()==columnCount())&&(getArray().length==elementCount());
	}
	
	@Override
	public double[] asDoubleArray() {
		if (isPackedArray()) return getArray();
		return null;
	}
}
