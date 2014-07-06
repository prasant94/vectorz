package mikera.vectorz.impl;

import java.nio.DoubleBuffer;
import java.util.Arrays;

import mikera.arrayz.INDArray;
import mikera.arrayz.impl.IDenseArray;
import mikera.vectorz.AVector;
import mikera.vectorz.Op;
import mikera.vectorz.Vector;
import mikera.vectorz.Vectorz;
import mikera.vectorz.util.DoubleArrays;
import mikera.vectorz.util.IntArrays;
import mikera.vectorz.util.VectorzException;

/**
 * Base class for mutable dense vectors backed by a double[] array with a fixed stride of 1
 * 
 * The double array can be directly accessed for performance purposes
 * 
 * @author Mike
 */
@SuppressWarnings("serial")
public abstract class ADenseArrayVector extends AStridedVector implements IDenseArray {
	
	protected ADenseArrayVector(int length, double[] data) {
		super(length,data);
	}

	/**
	 * ADenseArrayVector has a fixed stride of 1, which enables efficient operations
	 * on arrays
	 */
	@Override
	public final int getStride() {
		return 1;
	}

	/**
	 * Returns a vector referencing a sub-vector of the current vector
	 * 
	 * @param offset
	 * @param length
	 * @return
	 */
	public AVector subVector(int offset, int length) {
		int len = checkRange(offset,length);
		if (length == 0) return Vector0.INSTANCE;
		if (length == len) return this;
		return new ArraySubVector(this, offset, length);
	}

	@Override
	public ArrayIndexScalar slice(int position) {
		checkIndex(position);
		return new ArrayIndexScalar(getArray(), getArrayOffset() + position);
	}
	
	@Override
	public AVector selectView(int... inds) {
		inds=inds.clone();
		IntArrays.add(inds,getArrayOffset());
		return IndexedArrayVector.wrap(getArray(), inds);
	}

	public boolean isPackedArray() {
		return (getArrayOffset() == 0) && (length() == getArray().length);
	}

	@Override
	public boolean isView() {
		// ArrayVector is usually a view
		return true;
	}

	@Override
	public boolean isZero() {
		return DoubleArrays.isZero(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public boolean isRangeZero(int start, int length) {
		return DoubleArrays.isZero(getArray(), getArrayOffset()+start, length);
	}

	@Override
	public void toDoubleBuffer(DoubleBuffer dest) {
		dest.put(getArray(), getArrayOffset(), length());
	}

	@Override
	public void getElements(double[] data, int offset) {
		System.arraycopy(getArray(), getArrayOffset(), data, offset, length());
	}
	
	@Override
	public ADenseArrayVector dense() {
		// we are already dense!
		return this;
	}

	@Override
	public void fillRange(int offset, int length, double value) {
		assert ((offset >= 0) && (length >= 0) && ((offset + length) <= length()));
		double[] arr = getArray();
		int off = getArrayOffset();
		Arrays.fill(arr, off + offset, off + offset + length, value);
	}

	@Override
	public void set(AVector a) {
		checkSameLength(a);
		a.getElements(getArray(), getArrayOffset());
	}

	@Override
	public void set(AVector a, int offset) {
		assert (offset >= 0);
		assert (offset + length() <= a.length());
		a.copyTo(offset, this, 0, length());
	}

	@Override
	public void setRange(int offset, double[] data, int dataOffset, int length) {
		System.arraycopy(data, dataOffset, getArray(), getArrayOffset()
				+ offset, length);
	}

	@Override
	public void setElements(double[] values, int offset, int length) {
		assert (length == this.length());
		System.arraycopy(values, offset, getArray(), getArrayOffset(), length);
	}

	@Override
	public abstract double get(int i);

	@Override
	public abstract void set(int i, double value);

	@Override
	public abstract double unsafeGet(int i);

	@Override
	public abstract void unsafeSet(int i, double value);

	@Override
	public void add(AVector src) {
		src.addToArray(0, getArray(), getArrayOffset(), length());
	}

	@Override
	public void add(ADenseArrayVector v) {
		add(v, 0);
	}

	@Override
	public void add(AVector src, int srcOffset) {
		src.addToArray(srcOffset, getArray(), getArrayOffset(), length);
	}

	@Override
	public void add(int offset, AVector src) {
		src.addToArray(0, getArray(), getArrayOffset() + offset, length);
	}

	public void add(int offset, ADenseArrayVector src) {
		int length = src.length();
		DoubleArrays.add(src.getArray(), src.getArrayOffset(), getArray(),
				offset + getArrayOffset(), length);
	}

	public void add(int offset, ADenseArrayVector src, int srcOffset, int length) {
		DoubleArrays.add(src.getArray(), src.getArrayOffset() + srcOffset,
				getArray(), offset + getArrayOffset(), length);
	}

	@Override
	public void addMultiple(AVector v, double factor) {
		int length = checkSameLength(v);
		v.addMultipleToArray(factor, 0, getArray(), getArrayOffset(), length);
	}

	@Override
	public void scaleAdd(double factor, double constant) {
		DoubleArrays.scaleAdd(getArray(), getArrayOffset(), length(), factor,
				constant);
	}

	@Override
	public void add(double constant) {
		DoubleArrays.add(getArray(), getArrayOffset(), length(), constant);
	}

	@Override
	public void addProduct(AVector a, int aOffset, AVector b, int bOffset,
			double factor) {
		int length = length();
		double[] array = getArray();
		int offset = getArrayOffset();
		a.addProductToArray(factor, aOffset, b, bOffset, array, offset, length);
	}

	@Override
	public void addToArray(int offset, double[] destData, int destOffset, int length) {
		double[] data = getArray();
		int dataOffset = getArrayOffset() + offset;
		DoubleArrays.add(data, dataOffset, destData, destOffset, length);
	}
	
	@Override
	public void addToArray(double[] dest, int offset, int stride) {
		double[] data = getArray();
		int dataOffset = getArrayOffset();

		for (int i = 0; i < length; i++) {
			dest[offset + i*stride] += data[dataOffset + i];
		}
	}
	
	@Override
	public void addProduct(AVector a, AVector b) {
		int len = length();
		assert (len == a.length());
		assert (len == b.length());
		double[] array = getArray();
		int offset = getArrayOffset();
		if (b instanceof ADenseArrayVector) {
			a.addProductToArray(1.0, 0, (ADenseArrayVector)b, 0, array, offset, len);
		} else {
			a.addProductToArray(1.0, 0, b, 0, array, offset, len);			
		}
	}	

	@Override
	public void addProduct(AVector a, AVector b, double factor) {
		if (factor==0) return;
		int len = length();
		assert (len == a.length());
		assert (len == b.length());
		double[] array = getArray();
		int offset = getArrayOffset();
		if (b instanceof ADenseArrayVector) {
			a.addProductToArray(factor, 0, (ADenseArrayVector)b, 0, array, offset, len);
		} else {
			a.addProductToArray(factor, 0, b, 0, array, offset, len);			
		}
	}

	@Override
	public void addMultipleToArray(double factor, int offset, double[] array,
			int arrayOffset, int length) {
		if (factor==0) return;
		double[] data = getArray();
		int dataOffset = getArrayOffset() + offset;

		for (int i = 0; i < length; i++) {
			array[i + arrayOffset] += factor * data[i + dataOffset];
		}
	}

	@Override
	public void addProductToArray(double factor, int offset, AVector other,
			int otherOffset, double[] array, int arrayOffset, int length) {
		if (factor==0) return;
		if (other instanceof ADenseArrayVector) {
			addProductToArray(factor, offset, (ADenseArrayVector) other,
					otherOffset, array, arrayOffset, length);
			return;
		}
		assert (offset >= 0);
		assert (offset + length <= length());
		double[] thisArray = getArray();
		offset += getArrayOffset();
		for (int i = 0; i < length; i++) {
			array[i + arrayOffset] += factor * thisArray[i + offset]
					* other.unsafeGet(i + otherOffset);
		}
	}

	@Override
	public void addProductToArray(double factor, int offset,
			ADenseArrayVector other, int otherOffset, double[] array,
			int arrayOffset, int length) {
		if (factor==0) return;
		assert (offset >= 0);
		assert (offset + length <= length());
		double[] otherArray = other.getArray();
		otherOffset += other.getArrayOffset();
		double[] thisArray = getArray();
		offset += getArrayOffset();
		for (int i = 0; i < length; i++) {
			array[i + arrayOffset] += factor * thisArray[i + offset]
					* otherArray[i + otherOffset];
		}
	}

	public void add(ADenseArrayVector src, int srcOffset) {
		src.checkRange(srcOffset,length);
		double[] vdata = src.getArray();
		double[] data = getArray();
		int offset = getArrayOffset();
		int voffset = src.getArrayOffset() + srcOffset;
		for (int i = 0; i < length; i++) {
			data[offset + i] += vdata[voffset + i];
		}
	}
	
	@Override
	public void add(double[] data, int offset) {
		DoubleArrays.add(data, offset, getArray(), getArrayOffset(), length);
	}

	@Override
	public void addAt(int i, double v) {
		assert ((i >= 0) && (i < length()));
		double[] data = getArray();
		int offset = getArrayOffset();
		data[i + offset] += v;
	}

	@Override
	public double dotProduct(double[] data, int offset) {
		return DoubleArrays.dotProduct(getArray(), getArrayOffset(), data, offset, length());
	}

	@Override
	public double dotProduct(AVector v) {
		int length = checkSameLength(v);
		
		if (v instanceof ADenseArrayVector) {
			ADenseArrayVector vv = (ADenseArrayVector) v;
			return DoubleArrays.dotProduct(getArray(), getArrayOffset(),
					vv.getArray(), vv.getArrayOffset(), length);
		} else {
			return v.dotProduct(this.getArray(), this.getArrayOffset());
		}
	}

	@Override
	public void abs() {
		DoubleArrays.abs(getArray(), getArrayOffset(), length());
	}

	@Override
	public void log() {
		DoubleArrays.log(getArray(), getArrayOffset(), length());
	}

	@Override
	public void exp() {
		DoubleArrays.exp(getArray(), getArrayOffset(), length());
	}

	@Override
	public void applyOp(Op op) {
		op.applyTo(getArray(), getArrayOffset(), length());
	}

	@Override
	public double elementSum() {
		return DoubleArrays.elementSum(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public double elementProduct() {
		return DoubleArrays.elementProduct(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public double elementMax(){
		return DoubleArrays.elementMax(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public double elementMin(){
		return DoubleArrays.elementMin(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public double maxAbsElement(){
		return DoubleArrays.elementMaxAbs(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public int minElementIndex(){
		return DoubleArrays.elementMinIndex(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public int maxElementIndex(){
		return DoubleArrays.elementMaxIndex(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public int maxAbsElementIndex(){
		return DoubleArrays.elementMaxAbsIndex(getArray(), getArrayOffset(), length());
	}

	@Override
	public long nonZeroCount() {
		return DoubleArrays
				.nonZeroCount(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public int[] nonZeroIndices() {
		return DoubleArrays.nonZeroIndices(getArray(), getArrayOffset(), length());
	}

	@Override
	public void square() {
		DoubleArrays.square(getArray(), getArrayOffset(), length());
	}
	
	@Override
	public Vector squareCopy() {
		double[] ds=toDoubleArray();
		DoubleArrays.square(ds);
		return Vector.wrap(ds);
	}

	@Override
	public void sqrt() {
		DoubleArrays.sqrt(getArray(), getArrayOffset(), length());
	}

	/**
	 * Sets each component of the vector to its sign value (-1, 0 or 1)
	 */
	public void signum() {
		DoubleArrays.signum(getArray(), getArrayOffset(), length());
	}

	@Override
	public void multiply(AVector v) {
		v.multiplyTo(getArray(), getArrayOffset());
	}

	@Override
	public void multiply(double[] src, int srcOffset) {
		int len = length();
		double[] cdata = getArray();
		int coffset = getArrayOffset();
		for (int i = 0; i < len; i++) {
			cdata[i + coffset] *= src[i + srcOffset];
		}
	}

	@Override
	public void multiplyTo(double[] dest, int destOffset) {
		DoubleArrays.arraymultiply(getArray(), getArrayOffset(), dest, destOffset,
				length());
	}

	@Override
	public void divide(AVector v) {
		v.divideTo(getArray(), getArrayOffset());
	}

	@Override
	public void divide(double[] data, int offset) {
		DoubleArrays.arraydivide(data, offset, getArray(), getArrayOffset(), length());
	}

	@Override
	public void divideTo(double[] data, int offset) {
		DoubleArrays.arraydivide(getArray(), getArrayOffset(), data, offset, length());
	}

	@Override
	public void copyTo(int start, AVector dest, int destOffset, int length) {
		if (dest instanceof ADenseArrayVector) {
			copyTo(start, (ADenseArrayVector) dest, destOffset, length);
			return;
		}
		double[] src = getArray();
		int off = getArrayOffset();
		for (int i = 0; i < length; i++) {
			dest.unsafeSet(destOffset + i, src[off + start + i]);
		}
	}

	public void copyTo(int offset, ADenseArrayVector dest, int destOffset, int length) {
		double[] src = getArray();
		int off = getArrayOffset();
		double[] dst = dest.getArray();
		System.arraycopy(src, off + offset, dst, dest.getArrayOffset()
				+ destOffset, length);
	}
	
	@Override
	public void copyTo(int offset, double[] dest, int destOffset, int length, int stride) {
		for (int i=0; i<length; i++) {
			dest[destOffset+i*stride]=data[i+offset];
		}
	}

	@Override
	public void copyTo(int offset, double[] dest, int destOffset, int length) {
		double[] src = getArray();
		int off = getArrayOffset();
		System.arraycopy(src, off + offset, dest, destOffset, length);
	}

	public void addMultiple(ADenseArrayVector v, double factor) {
		int length = checkSameLength(v);
		v.addMultipleToArray(factor, 0, getArray(), getArrayOffset(), length);
	}
	
	@Override
	public void addMultiple(int offset, AVector src, int srcOffset, int length, double factor) {
		checkRange(offset,length);
		src.checkRange(srcOffset,length);
		if (factor==0.0) return;
		int tOffset=offset+this.getArrayOffset();
		src.addMultipleToArray(factor, srcOffset, this.getArray(), tOffset, length);
	}

	@Override
	public double magnitudeSquared() {
		return DoubleArrays.elementSquaredSum(data, getArrayOffset(), length);
	}

	@Override
	public double magnitude() {
		return Math.sqrt(magnitudeSquared());
	}

	@Override
	public void fill(double value) {
		int offset = getArrayOffset();
		Arrays.fill(getArray(), offset, offset + length, value);
	}

	@Override
	public void pow(double exponent) {
		int len = length();
		double[] data = getArray();
		int offset = getArrayOffset();
		DoubleArrays.pow(data, offset, len, exponent);
	}

	@Override
	public void reciprocal() {
		DoubleArrays.reciprocal(getArray(), getArrayOffset(), length());
	}

	@Override
	public void clamp(double min, double max) {
		DoubleArrays.clamp(getArray(), getArrayOffset(), length(), min, max);
	}

	@Override
	public void multiply(double factor) {
		DoubleArrays.multiply(getArray(), getArrayOffset(), length(), factor);
	}

	@Override
	public AVector tryEfficientJoin(AVector v) {
		if (v instanceof ADenseArrayVector) return join((ADenseArrayVector) v);
		if (v instanceof JoinedArrayVector) return join((JoinedArrayVector) v);
		return null;
	}

	public AVector join(ADenseArrayVector v) {
		if ((v.getArray() == getArray())
				&& ((getArrayOffset() + length) == v.getArrayOffset())) { 
			return Vectorz.wrap(getArray(), getArrayOffset(), length + v.length()); 
		}
		return JoinedArrayVector.joinVectors(this, v);
	}

	public JoinedArrayVector join(JoinedArrayVector v) {
		return JoinedArrayVector.wrap(this).join(v);
	}
	
	@Override
	public boolean equals(INDArray v) {
		if (v.dimensionality()!=1) return false;
		int len=length();
		if (len != v.getShape(0)) return false;
		
		return v.equalsArray(getArray(), getArrayOffset());
	}

	@Override
	public boolean equals(AVector v) {
		if (v==this) return true;
		int len = length();
		if (v.length() != len) return false;
		return v.equalsArray(getArray(), getArrayOffset());
	}

	@Override
	public boolean equalsArray(double[] data, int offset) {
		return DoubleArrays.equals(data, offset, getArray(), getArrayOffset(), length);
	}

	@Override
	public boolean equalsArray(double[] data) {
		if (length() != data.length) return false;
		return equalsArray(data,0);
	}
	
	@Override
	public boolean elementsEqual(double value) {
		int length=length();
		double[] data = getArray();
		int offset = getArrayOffset();
		return DoubleArrays.elementsEqual(data,offset,length,value);
	}

	@Override
	public void validate() {
		int length = length();
		double[] data = getArray();
		int offset = getArrayOffset();
		if ((offset < 0) || (offset + length > data.length)) {
			throw new VectorzException("ArrayVector out of bounds");
		}
		super.validate();
	}
}
