package mikera.matrixx.impl;

import mikera.arrayz.impl.IStridedArray;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrixx;
import mikera.vectorz.Vectorz;
import mikera.vectorz.impl.AStridedVector;
import mikera.vectorz.util.ErrorMessages;

public abstract class AStridedMatrix extends AArrayMatrix implements IStridedArray {
	private static final long serialVersionUID = -8908577438753599161L;

	protected AStridedMatrix(double[] data, int rows, int cols) {
		super(data, rows, cols);
	}

	public abstract int getArrayOffset();

	public abstract int rowStride();
	
	public abstract int columnStride();	
	
	@Override
	public AStridedMatrix subMatrix(int rowStart, int rows, int colStart, int cols) {
		if ((rowStart<0)||(rowStart>=rows)||(colStart<0)||(colStart>=cols)) throw new IndexOutOfBoundsException(ErrorMessages.position(rowStart,colStart));
		if ((rowStart+rows>this.rows)||(colStart+cols>this.cols)) throw new IndexOutOfBoundsException(ErrorMessages.position(rowStart+rows,colStart+cols));
		if ((rows<1)||(cols<1)) throw new IllegalArgumentException(ErrorMessages.illegalSize(rows,cols));
		return StridedMatrix.wrap(data, rows, cols, 
				getArrayOffset()+rowStart*rowStride()+colStart*columnStride(), this.rowStride(), this.columnStride());
	}
	
	@Override
	public AStridedVector getRowView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i*rowStride(), cols, columnStride());
	}
	
	@Override
	public AStridedVector getColumnView(int i) {
		return Vectorz.wrapStrided(data, getArrayOffset()+i*columnStride(), rows, rowStride());
	}
	
	@Override
	public void copyRowTo(int row, double[] dest, int destOffset) {
		int cc=columnCount();
		for (int i=0; i<cc; i++) {
			dest[i+destOffset]=unsafeGet(row,i);
		}
	}
	
	@Override
	public void copyColumnTo(int col, double[] dest, int destOffset) {
		int rc=rowCount();
		for (int i=0; i<rc; i++) {
			dest[i+destOffset]=unsafeGet(i,col);
		}
	}
	
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
