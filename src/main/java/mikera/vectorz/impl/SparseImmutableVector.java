package mikera.vectorz.impl;

import mikera.indexz.Index;
import mikera.matrixx.AMatrix;
import mikera.matrixx.impl.AVectorMatrix;
import mikera.vectorz.AVector;
import mikera.vectorz.Op;
import mikera.vectorz.Vector;
import mikera.vectorz.util.DoubleArrays;
import mikera.vectorz.util.ErrorMessages;
import mikera.vectorz.util.VectorzException;

/**
 * Indexed sparse immutable vector
 * 
 * Efficient for sparse vectors. Maintains a indexed array of strictly non-zero elements. 
 * 
 * @author Mike
 *
 */
public class SparseImmutableVector extends ASparseIndexedVector {
	private static final long serialVersionUID = 750093598603613879L;

	private final Index index;
	private final double[] data;
	private final int dataLength;
	
	private SparseImmutableVector(int length, Index index) {
		this(length,index,new double[index.length()]);
	}
	
	private SparseImmutableVector(int length, Index index, double[] data) {
		super(length);
		this.index=index;
		this.data=data;
		dataLength=data.length;
	}
	
	private SparseImmutableVector(int length, Index index, AVector data) {
		this(length,index,data.toDoubleArray());
	}
	
	/**
	 * Creates a SparseImmutableVector with the specified index and data values.
	 * 
	 * WARNING: Performs no checking - Index must be distinct and sorted, and data must be non-zero.
	 */
	public static SparseImmutableVector wrap(int length, Index index, double[] data) {
		assert(index.length()==data.length);
		assert(index.isDistinctSorted());
		return new SparseImmutableVector(length, index,data);
	}
	
	@Override
	double[] internalData() {
		return data;
	}

	@Override
	Index internalIndex() {
		return index;
	}
	
	/**
	 * Creates a SparseImmutableVector using the given sorted Index to identify the indexes of non-zero values,
	 * and a double[] array to specify all the non-zero element values
	 */
	public static AVector create(int length, Index index, double[] data) {
		int dataLength=data.length;
		if (!index.isDistinctSorted()) {
			throw new IllegalArgumentException("Index must be sorted and distinct");
		}
		if (!(index.length()==dataLength)) {
			throw new IllegalArgumentException("Length of index: mismatch woth data");			
		}
		if (dataLength==0) return ZeroVector.create(length);
		if (dataLength==length) return ImmutableVector.create(data);
		return new SparseImmutableVector(length, index.clone(),DoubleArrays.copyOf(data));
	}
	
	/**
	 * Creates a SparseImmutableVector using the given sorted Index to identify the indexes of non-zero values,
	 * and a dense vector to specify all the non-zero element values
	 */
	public static AVector create(int length, Index index, AVector data) {
		int dataLength=data.length();
		if (!index.isDistinctSorted()) {
			throw new IllegalArgumentException("Index must be sorted and distinct");
		}
		if (!(index.length()==dataLength)) {
			throw new IllegalArgumentException("Length of index: mismatch woth data");			
		}
		if (dataLength==0) return ZeroVector.create(length);
		if (dataLength==length) return ImmutableVector.create(data);
		return wrap(length, index.clone(), data.toDoubleArray());
	}
	
	/** 
	 * Creates a SparseIndexedVector from the given vector, ignoring the zeros in the source.
	 * 
	 */
	public static AVector create(AVector source) {
		int length = source.length();
		if (length==0) return Vector0.INSTANCE;
		int dataLength=(int) source.nonZeroCount();
		if (dataLength==length) return ImmutableVector.create(source);
		if (dataLength==0) return ZeroVector.create(length);
		
		int[] indexes=new int[dataLength];
		double[] vals=new double[dataLength];
		int pos=0;
		for (int i=0; i<length; i++) {
			double v=source.unsafeGet(i);
			if (v!=0.0) {
				indexes[pos]=i;
				vals[pos]=v;
				pos++;
			}
		}
		return wrap(length,Index.wrap(indexes),vals); 
	}
	
	/** Creates a SparseIndexedVector from a row of an existing matrix */
	public static AVector createFromRow(AMatrix m, int row) {
		if (m instanceof AVectorMatrix) return create(m.getRow(row));
		return create(m.getRow(row));
	}
	
	@Override
	public int nonSparseElementCount() {
		return dataLength;
	}
	
	@Override
	public boolean isZero() {
		// never zero, since we maintain invariant of always having one non-zero value
		return false;
	}
	
	@Override
	public double maxAbsElement() {
		double result=data[0];
		for (int i=1; i<dataLength; i++) {
			double d=Math.abs(data[i]);
			if (d>result) result=d; 
		}
		return result;
	}
	
	@Override
	public int maxElementIndex(){
		double result=data[0];
		int di=0;
		for (int i=1; i<dataLength; i++) {
			double d=data[i];
			if (d>result) {
				result=d; 
				di=i;
			}
		}
		if (result<0.0) { // need to find a sparse element
			int ind=sparseElementIndex();
			if (ind>0) return ind;
		}
		return index.get(di);
	}
	
 
	@Override
	public int maxAbsElementIndex(){
		double result=Math.abs(data[0]);
		int di=0;
		for (int i=1; i<dataLength; i++) {
			double d=Math.abs(data[i]);
			if (d>result) {
				result=d; 
				di=i;
			}
		}
		return index.get(di);
	}
	
	@Override
	public int minElementIndex(){
		double result=data[0];
		int di=0;
		for (int i=1; i<dataLength; i++) {
			double d=data[i];
			if (d<result) {
				result=d; 
				di=i;
			}
		}
		if (result>0.0) { // need to find a sparse element
			int ind=sparseElementIndex();
			if (ind>=0) return ind;
		}
		return index.get(di);
	}
	
	/**
	 * Return this index of a sparse zero element, or -1 if not sparse
	 * @return
	 */
	private int sparseElementIndex() {
		for (int i=0; i<length; i++) {
			if (!index.contains(i)) return i;
		}
		throw new VectorzException(ErrorMessages.impossible());
	}

	
	@Override
	public void negate() {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public void applyOp(Op op) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public void abs() {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}

	@Override
	public double get(int i) {
		checkIndex(i);
		return unsafeGet(i);
	}
	
	@Override
	public double unsafeGet(int i) {
		int ip=index.indexPosition(i);
		if (ip<0) return 0.0;
		return data[ip];
	}
	
	@Override
	public boolean isFullyMutable() {
		return false;
	}
	
	@Override
	public boolean isMutable() {
		return false;
	}
	
	@Override
	public long nonZeroCount() {
		return dataLength;
	}
	
	@Override
	public int[] nonZeroIndices() {
		return index.data.clone();
	}

	@Override
	public void add(ASparseVector v) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public void set(AVector v) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}

	@Override
	public void set(int i, double value) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));

	}
	
	@Override
	public void addAt(int i, double value) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}

	@Override
	public Vector nonSparseValues() {
		return Vector.wrap(data);
	}
	
	@Override
	protected double[] nonZeroValues() {
		return DoubleArrays.copyOf(data);
	}
	
	@Override
	public Index nonSparseIndexes() {
		return index;
	}

	@Override
	public boolean includesIndex(int i) {
		return index.indexPosition(i)>=0;
	}
	
	@Override
	public Vector dense() {
		Vector v=Vector.createLength(length);
		addToArray(v.data,0);
		return v;
	}
	
	@Override
	public SparseIndexedVector mutable() {
		return SparseIndexedVector.create(length, index, data);
	}
	
	@Override
	public SparseIndexedVector clone() {
		return SparseIndexedVector.create(length, index, data);
	}
	
	@Override
	public SparseIndexedVector sparseClone() {
		return SparseIndexedVector.create(length, index, data);
	}
	
	@Override
	public SparseImmutableVector exactClone() {
		return new SparseImmutableVector(length,index.clone(),data.clone());
	}
	
	@Override
	public void validate() {
		if (data.length==0) throw new VectorzException("SparseImmutableVector must have some non-zero values");
		if (index.length()!=data.length) throw new VectorzException("Inconsistent data and index!");
		if (!index.isDistinctSorted()) throw new VectorzException("Invalid index: "+index);
		for (int i=0; i<data.length; i++) {
			if (data[i]==0) throw new VectorzException("Should be no zero values in data array!");
		}
		super.validate();
	}


}
