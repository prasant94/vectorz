package mikera.vectorz.impl;

import java.io.ObjectStreamException;
import java.util.Iterator;

import mikera.arrayz.INDArray;
import mikera.indexz.Index;
import mikera.matrixx.AMatrix;
import mikera.randomz.Hash;
import mikera.vectorz.AVector;
import mikera.vectorz.Scalar;
import mikera.vectorz.Vector;
import mikera.vectorz.Vectorz;
import mikera.vectorz.util.DoubleArrays;
import mikera.vectorz.util.ErrorMessages;
import mikera.vectorz.util.IntArrays;

/**
 * Specialised immutable vector containing nothing but zeros.
 * 
 * Must have length 1 or more: use Vector0 instead for immutable length 0 vectors.
 * 
 * @author Mike
 */
public final class ZeroVector extends ASparseVector {
	private static final long serialVersionUID = -7928191943246067239L;
	
	private static final int ZERO_VECTOR_CACHE_SIZE=30;
	private static final ZeroVector[] ZERO_VECTORS = new ZeroVector[ZERO_VECTOR_CACHE_SIZE];
	private static final Double ZERO_DOUBLE=0.0;
	
	private static ZeroVector last=new ZeroVector(ZERO_VECTOR_CACHE_SIZE);
	
	static {
		for (int i=1; i<ZERO_VECTOR_CACHE_SIZE; i++) {
			ZERO_VECTORS[i]=new ZeroVector(i);
		}
	}
	
	private ZeroVector(int dimensions) {
		super(dimensions);
	}
	
	public static ZeroVector create(int dimensions) {
		if (dimensions<=0) throw new IllegalArgumentException("Can't create length "+dimensions+" ZeroVector. Use Vector0 instead");
		return new ZeroVector(dimensions);
	}
	
	public static ZeroVector createNew(int dimensions) {
		if (dimensions<=0) throw new IllegalArgumentException("Can't create length "+dimensions+" ZeroVector. Use Vector0 instead");
		return new ZeroVector(dimensions);
	}
	
	public static ZeroVector createCached(int dimensions) {
		if (dimensions<=0) throw new IllegalArgumentException("Can't create length "+dimensions+" ZeroVector. Use Vector0 instead");
		ZeroVector zv=tryCreate(dimensions);
		if (zv!=null) return zv;
		zv= new ZeroVector(dimensions);
		last=zv;
		return zv;
	}
	
	/**
	 * Creates a ZeroVector with the same number of elements as the given array.
	 * @param arraySize
	 * @return
	 */
	public static AVector create(INDArray array) {
		int n=Vectorz.safeLongToInt(array.elementCount());
		return ZeroVector.create(n);
	}

	
	private static ZeroVector tryCreate(int dimensions) {
		if (dimensions<ZERO_VECTOR_CACHE_SIZE) {
			return ZERO_VECTORS[dimensions];
		}
		if (dimensions==last.length) return last;
		return null;
	}

	@Override
	public double dotProduct(AVector v) {
		checkSameLength(v);
		return 0.0;
	}
	
	@Override
	public double dotProduct(double[] data, int offset) {
		return 0.0;
	}
	
	@Override
	public AVector innerProduct(AMatrix m) {
		checkLength(m.rowCount());
		return ZeroVector.create(m.columnCount());
	}
	
	@Override
	public ImmutableScalar innerProduct(AVector a) {
		checkSameLength(a);
		return ImmutableScalar.ZERO;
	}
	
	@Override
	public Scalar innerProduct(Vector v) {
		checkSameLength(v);
		return Scalar.create(0.0);
	}
	
	@Override
	public ZeroVector innerProduct(double a) {
		return this;
	}

	@Override
	public double get(int i) {
		checkIndex(i);
		return 0.0;
	}

	@Override
	public void set(int i, double value) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public double unsafeGet(int i) {
		return 0.0;
	}

	@Override
	public void unsafeSet(int i, double value) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public void add(ASparseVector v) {
		throw new UnsupportedOperationException(ErrorMessages.immutable(this));
	}
	
	@Override
	public AVector addCopy(AVector a) {
		checkSameLength(a);
		return a.copy();
	}
	
	@Override
	public AVector subCopy(AVector a) {
		checkSameLength(a);
		return a.negateCopy();
	}
	
	@Override
	public double magnitudeSquared() {
		return 0.0;
	}
	
	@Override
	public double magnitude() {
		return 0.0;
	}
	
	@Override
	public double elementSum() {
		return 0.0;
	}
	
	@Override
	public double elementProduct() {
		return 0.0;
	}
	
	@Override
	public double elementMax(){
		return 0.0;
	}
	
	@Override
	public double elementMin(){
		return 0.0;
	}
	
	@Override
	public int maxElementIndex(){
		return 0;
	}
	
	@Override
	public double maxAbsElement(){
		return 0.0;
	}
	
	@Override
	public int maxAbsElementIndex(){
		return 0;
	}
	
	@Override
	public int minElementIndex(){
		return 0;
	}
	
	@Override
	public long nonZeroCount() {
		return 0;
	}

	@Override
	public boolean isZero() {
		return true;
	}
	
	@Override
	public boolean isRangeZero(int start, int length) {
		return true;
	}
	
	@Override
	public boolean isBoolean() {
		return true;
	}
	
	@Override
	public boolean isMutable() {
		return false;
	}
	
	@Override
	public boolean isFullyMutable() {
		return false;
	}
	
	@Override
	public boolean isUnitLengthVector() {
		return false;
	}
	
	@Override
	public void addToArray(int offset, double[] array, int arrayOffset, int length) {
		// do nothing!
	}
	
	@Override
	public void copyTo(int offset, double[] dest, int destOffset, int length, int stride) {
		for (int i=0; i<length; i++) {
			dest[destOffset+i*stride]=0.0;
		}
	}
	
	@Override
	public void addToArray(double[] dest, int offset, int stride) {
		// do nothing!
	}
	
	@Override
	public void addMultipleToArray(double factor,int offset, double[] array, int arrayOffset, int length) {
		// do nothing!
	}
	
	@Override
	public final ImmutableScalar slice(int i) {
		checkIndex(i);
		return ImmutableScalar.ZERO;
	}
	
	@Override
	public Iterator<Double> iterator() {
		return new RepeatedElementIterator(length,ZERO_DOUBLE);
	}

	@Override
	public double density() {
		return 0.0;
	}
	
	@Override
	public AVector subVector(int offset, int length) {
		int len=checkRange(offset,length);
		if (length==0) return Vector0.INSTANCE;
		if (length==len) return this;
		return ZeroVector.create(length);
	}

	public ZeroVector join(ZeroVector a) {
		return ZeroVector.create(length+a.length);
	}
	
	@Override
	public AVector tryEfficientJoin(AVector a) {
		if (a instanceof ZeroVector) {
			return join((ZeroVector)a);
		} else if (a instanceof AxisVector) {
			AxisVector av=(AxisVector)a;
			return AxisVector.create(av.getAxis()+length, av.length()+length);
		} else if (a instanceof SingleElementVector) {
			SingleElementVector sev=(SingleElementVector)a;
			return SingleElementVector.create(sev.value, length+sev.index, sev.length+length);
		}
		return null;
	}
	
	@Override
	public AVector reorder(int[] order) {
		int n=order.length;
		if (n==length) return this;
		return createNew(n);
	}	
	
	@Override
	public AVector reorder(int dim, int[] order) {
		if (dim!=0) throw new IndexOutOfBoundsException(ErrorMessages.invalidDimension(this, dim));
		return reorder(order);
	}	
	
	/**
	 * readResolve method to ensure we always use the singleton
	 */
	private Object readResolve() throws ObjectStreamException {
		ZeroVector zv=tryCreate(length);
		if (zv!=null) return zv;
		return this;
	}

	@Override
	public int nonSparseElementCount() {
		return 0;
	}

	@Override
	public AVector nonSparseValues() {
		return Vector0.INSTANCE;
	}

	@Override
	public Index nonSparseIndexes() {
		return Index.EMPTY;
	}
	
	@Override
	public int[] nonZeroIndices() {
		return IntArrays.EMPTY_INT_ARRAY;
	}
	
	@Override
	public AVector squareCopy() {
		return this;
	}

	@Override
	public boolean includesIndex(int i) {
		return false;
	}
	
	@Override
	public int hashCode() {
		return Hash.zeroVectorHash(length);
	}
	
	@Override
	public ZeroVector exactClone() {
		return new ZeroVector(length);
	}
	
	@Override
	public AVector sparseClone() {
		return Vectorz.createSparseMutable(length);
	}
	
	@Override
	public double[] toDoubleArray() {
		return new double[length];
	}
	
	@Override
	public AVector selectClone(int... inds) {
		return Vectorz.newVector(inds.length);
	}
	
	@Override
	public AVector selectView(int... inds) {
		return Vectorz.createZeroVector(inds.length);
	}
	
	@Override
	public boolean equals(AVector v) {
		if (!isSameShape(v)) return false;
		return v.isZero();
	}
	
	@Override
	public boolean equalsArray(double[] data, int offset) {
		return DoubleArrays.isZero(data, offset, length);
	}
	
	@Override
	public boolean elementsEqual(double value) {
		return value==0.0;
	}

	@Override
	public boolean hasUncountable() {
		return false;
	}

}
