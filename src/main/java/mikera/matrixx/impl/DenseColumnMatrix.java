package mikera.matrixx.impl;

import java.util.Arrays;

import mikera.arrayz.impl.IDense;
import mikera.matrixx.AMatrix;
import mikera.matrixx.Matrix;
import mikera.vectorz.AVector;
import mikera.vectorz.Op;
import mikera.vectorz.impl.ArraySubVector;
import mikera.vectorz.util.DoubleArrays;

/**
 * A densely packed matrix organised in column-major format.
 * 
 * Transposes to/from a regular dense Matrix
 * 
 * @author Mike
 *
 */
public class DenseColumnMatrix extends AStridedMatrix implements IFastColumns, IDense {
	private static final long serialVersionUID = 5459617932072332096L;

	private DenseColumnMatrix(int rowCount, int columnCount, double[] data) {
		super(data, rowCount, columnCount);
	}
	
	private DenseColumnMatrix(int rowCount, int columnCount) {
		this(rowCount, columnCount, Matrix.createStorage(rowCount, columnCount));
	}
	
	public static DenseColumnMatrix wrap(int rows, int cols, double[] data) {
		return new DenseColumnMatrix(rows, cols, data);
	}
	
	@Override
	public int getArrayOffset() {
		return 0;
	}

	@Override
	public int rowStride() {
		return 1;
	}

	@Override
	public int columnStride() {
		return rows;
	}
	
	@Override
	public ArraySubVector getColumnView(int j) {
		return ArraySubVector.wrap(data, j*rows, rows);
	}

	@Override
	public void copyRowTo(int i, double[] dest, int destOffset) {
		for (int j=0; j<cols; j++) {
			dest[destOffset+j]=data[i+j*rows];
		}
	}

	@Override
	public void copyColumnTo(int j, double[] dest, int destOffset) {
		System.arraycopy(data, j*rows, dest, destOffset, rows);
	}
	
	@Override
	public void setRow(int i, AVector row) {
		int cc = checkColumnCount(row.length());
		for (int j = 0; j < cc; i++) {
			data[index(i, j)] = row.unsafeGet(i);
		}
	}

	@Override
	public void setColumn(int j, AVector col) {
		int rc = checkRowCount(col.length());
		col.getElements(data, j * rc);
	}
	
	@Override
	public void addMultiple(AMatrix m, double factor) {
		checkRowCount(m.rowCount());
		int cc=checkColumnCount(m.columnCount());
		
		for (int i=0; i<cc; i++) {
			getColumnView(i).addMultiple(m.getColumn(i), factor);
		}
	}

	@Override
	protected int index(int i, int j) {
		return i+j*rows;
	}
	
	@Override
	public double get(int i, int j) {
		if ((i < 0) || (i >= rows)) {
			// we only need to check i is in range: out of range j will trigger exception anyway
			throw new IndexOutOfBoundsException("Row: "+i);
		}
		return data[(j * rows) + i];
	}
	
	@Override
	public void set(int i, int j, double value) {
		if ((i < 0) || (i >= rows)) {
			// we only need to check i is in range: out of range j will trigger exception anyway
			throw new IndexOutOfBoundsException("Row: "+i);
		}
		data[(j * rows) + i] = value;
	}

	@Override
	public void unsafeSet(int i, int j, double value) {
		data[(j * rows) + i] = value;
	}

	@Override
	public double unsafeGet(int i, int j) {
		return data[(j * rows) + i];
	}
	
	@Override
	public void addAt(int i, int j, double d) {
		data[(j * rows) + i] += d;
	}

	@Override
	public boolean isFullyMutable() {
		return true;
	}
	
	@Override
	public boolean isPackedArray() {
		return (cols<=1);
	}
	
	@Override
	public boolean isBoolean() {
		return DoubleArrays.isBoolean(data, 0, data.length);
	}

	@Override
	public boolean isZero() {
		return DoubleArrays.isZero(data, 0, data.length);
	}
	
	@Override
	public double elementSum() {
		return DoubleArrays.elementSum(data);
	}

	@Override
	public double elementSquaredSum() {
		return DoubleArrays.elementSquaredSum(data);
	}

	@Override
	public double elementMax() {
		return DoubleArrays.elementMax(data);
	}

	@Override
	public double elementMin() {
		return DoubleArrays.elementMin(data);
	}

	@Override
	public void abs() {
		DoubleArrays.abs(data);
	}

	@Override
	public void signum() {
		DoubleArrays.signum(data);
	}

	@Override
	public void square() {
		DoubleArrays.square(data);
	}

	@Override
	public void exp() {
		DoubleArrays.exp(data);
	}

	@Override
	public void log() {
		DoubleArrays.log(data);
	}
	
	@Override
	public void applyOp(Op op) {
		op.applyTo(data);
	}

	@Override
	public long nonZeroCount() {
		return DoubleArrays.nonZeroCount(data);
	}
	
	@Override
	public void add(double d) {
		DoubleArrays.add(data, d);
	}

	@Override
	public void multiply(double factor) {
		DoubleArrays.multiply(data, factor);
	}
	
	@Override
	public void set(double value) {
		Arrays.fill(data, value);
	}
	
	@Override
	public void reciprocal() {
		DoubleArrays.reciprocal(data, 0, data.length);
	}

	@Override
	public void clamp(double min, double max) {
		DoubleArrays.clamp(data, 0, data.length, min, max);
	}

	
	@Override
	public Matrix getTranspose() {
		return getTransposeView();
	}
	
	@Override
	public Matrix getTransposeView() {
		return Matrix.wrap(cols, rows, data);
	}

	@Override
	public DenseColumnMatrix exactClone() {
		return new DenseColumnMatrix(rows,cols,data.clone());
	}
	
	@Override
	public DenseColumnMatrix dense() {
		return this;
	}
	
	@Override
	public DenseColumnMatrix copy() {
		return exactClone();
	}
	
	@Override
	public DenseColumnMatrix clone() {
		return exactClone();
	}
	
	@Override
	public Matrix toMatrixTranspose() {
		return Matrix.wrap(cols, rows, data);
	}

}
