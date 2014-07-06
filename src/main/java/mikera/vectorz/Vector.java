package mikera.vectorz;

import java.nio.DoubleBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import mikera.arrayz.INDArray;
import mikera.indexz.AIndex;
import mikera.indexz.Index;
import mikera.randomz.Hash;
import mikera.vectorz.impl.ADenseArrayVector;
import mikera.vectorz.util.DoubleArrays;
import mikera.vectorz.util.ErrorMessages;
import mikera.vectorz.util.VectorzException;


/**
 * General purpose vector of arbitrary length, backed by an densely packed double[] array.
 * 
 * This is the most efficient type for general purpose dense 1D vectors.
 * 
 * @author Mike
 *
 */
public final class Vector extends ADenseArrayVector {
	private static final long serialVersionUID = 6283741614665875877L;

	public static final Vector EMPTY = wrap(DoubleArrays.EMPTY);

	private Vector(double... values) {
		super(values.length,values);
	}
	
	private Vector(Object... values) {
		super(values.length,new double[values.length]);
		for (int i=0; i<length; i++) {
			data[i]=Tools.toDouble(values[i]);
		}
	}

	private Vector(int length) {
		this(new double[length]);
	}

	/**
	 * Copy constructor from an arbitrary vector
	 * 
	 * @param source
	 */
	public Vector(AVector source) {
		this(source.toDoubleArray());
	}
	
	/**
	 * Wraps a double array into a Vector, does *no defensive copy* so use with caution
	 * @param source
	 * @return
	 */
	public static Vector wrap(double[] source) {
		return new Vector(source);
	}
	
	/**
	 * Creates a new Vector from the given double[] data. Takes a defensive copy.
	 */
	public static Vector create(double[] data) {
		if (data.length==0) return EMPTY;
		return wrap(data.clone());
	}
	
	public static Vector create(ArrayList<Double> al) {
		int n=al.size();
		Vector v=Vector.createLength(n);
		for (int i=0; i<n; i++) {
			v.unsafeSet(i,al.get(i));
		}
		return v;
	}
	
	public static Vector create(List<Double> al) {
		int n=al.size();
		Vector v=Vector.createLength(n);
		for (int i=0; i<n; i++) {
			v.unsafeSet(i,al.get(i));
		}
		return v;
	}
	
	/**
	 * Creates a new Vector from the given double[] data. Takes a defensive copy.
	 */
	public static Vector create(double[] data, int start, int length) {
		return wrap(DoubleArrays.copyOf(data,start,length));
	}
	
	/**
	 * Creates a new vector using the elements in the specified vector.
	 * Truncates or zero-pads the data as required to fill the new vector
	 * @param data
	 * @return
	 */
	public static Vector createFromVector(AVector source, int length) {
		Vector v=Vector.createLength(length);
		int n=Math.min(length, source.length());
		source.copyTo(0, v.data, 0, n);
		return v;
	}
	
	/**
	 * Create a vector with specific component values. 
	 * Creates a copy of the value array
	 * @param values
	 * @return
	 */
	public static Vector of(double... values) {
		return create(values);
	}
	
	/**
	 * Create an empty (zero-filled) vector of a specified length
	 * @param length
	 * @return
	 */
	public static Vector createLength(int length) {
		if (length<1) {
		  if (length<0) throw new IllegalArgumentException(ErrorMessages.illegalSize(length));
		  return EMPTY;
		}
		return new Vector(length);
	}

	public static Vector create(AVector a) {
		return new Vector(a.toDoubleArray());
	}
	
	public static Vector create(INDArray a) {
		return new Vector(a.toDoubleArray());
	}
	
	public static Vector create(AIndex a) {
		int n=a.length();
		Vector v=createLength(n);
		for (int i=0; i<n; i++) {
			v.unsafeSet(i,a.get(i));
		}
		return v;
	}

	@Override
	public double get(int i) {
		return data[i];
	}
	
	@Override
	public double unsafeGet(int i) {
		return data[i];
	}

	@Override
	public void set(int i, double value) {
		data[i]=value;
	}
	
	@Override
	public void unsafeSet(int i, double value) {
		data[i]=value;
	}
	
	@Override
	public void setRange(int offset, double[] data, int dataOffset, int length) {
		System.arraycopy(data, dataOffset, this.data, offset, length);
	}
	
	@Override
	public void set(AVector a) {
		if (a instanceof Vector) {
			if (a==this) return;
			Vector v=(Vector)a;
			System.arraycopy(v.data, 0, data, 0, data.length);
		} else {
			super.set(a);
		}
	}
	
	@Override
	public void getElements(double[] dest, int offset) {
		System.arraycopy(data, 0, dest, offset, data.length);
	}

	@Override
	public int getArrayOffset() {
		return 0;
	}
	
	@Override
	public void applyOp(Op op) {
		op.applyTo(data, 0, data.length);
	}
	
	@Override
	public void fill(double value) {
		Arrays.fill(data, value);
	}
	
	@Override
	public void clamp(double min, double max) {
		DoubleArrays.clamp(data,min,max);
	}
	
	@Override
	public void square() {
		DoubleArrays.square(data);	
	}
	
	@Override
	public void tanh() {
		DoubleArrays.tanh(data);		
	}
	
	@Override
	public void logistic() {
		DoubleArrays.logistic(data);		
	}
	
	@Override
	public double elementSum() {
		return DoubleArrays.elementSum(data);
	}
	
	@Override
	public double elementMax(){
		return DoubleArrays.elementMax(data);
	}
	
	@Override
	public double elementMin(){
		return DoubleArrays.elementMin(data);
	}
	
	@Override
	public long nonZeroCount() {
		return DoubleArrays.nonZeroCount(data);
	}	
	
	@Override
	public void signum() {
		DoubleArrays.signum(data);
	}
	
	@Override
	public void abs() {
		DoubleArrays.abs(data);
	}
	
	@Override
	public void add(ADenseArrayVector src, int srcOffset) {
		int length=length();
		assert(srcOffset>=0);
		assert(srcOffset+length<=src.length());
		double[] vdata=src.getArray();
		int voffset=src.getArrayOffset()+srcOffset;
		for (int i = 0; i < length; i++) {
			data[i] += vdata[voffset + i];
		}
	}
	
	@Override
	public void addMultiple(ADenseArrayVector v, double factor) {
		int length=checkSameLength(v);
		v.addMultipleToArray(factor,0,data, 0,length);
	}
	
	@Override
	public void add(AVector v) {
		checkSameLength(v);
		v.addToArray(data, 0);
	}
	
	@Override
	public void add(Vector v) {
		checkSameLength(v);
		add(v.data,0);
	}
	
	@Override
	public void add(double[] srcData, int srcOffset) {
		int length=length();
		DoubleArrays.add(srcData, srcOffset, data, 0, length);
	}
	
	@Override
	public void scaleAdd(double factor, double constant) {
		int length=length();
		for (int i=0; i<length; i++) {
			data[i]=(factor*data[i])+constant;
		}
	}

	@Override
	public void add(double constant) {
		DoubleArrays.add(data, 0, data.length, constant);
	}
	
	@Override
	public void addProduct(AVector a, AVector b) {
		if((a instanceof Vector)&&(b instanceof Vector)) {
			addProduct((Vector)a,(Vector)b);
			return;
		}
		super.addProduct(a,b);
	}
	
	public void addProduct(Vector a, Vector b) {
		int length=checkSameLength(a,b);
		for (int i = 0; i < length; i++) {
			data[i]+=(a.data[i]*b.data[i]);
		}
	}
	
	public void addProduct(Vector a, Vector b, double factor) {
		int length=checkSameLength(a,b);
		for (int i = 0; i < length; i++) {
			data[i]+=(a.data[i]*b.data[i])*factor;
		}
	}
	
	@Override
	public void addAt(int i, double v) {
		data[i]+=v;
	}
	
	@Override
	public void sub(AVector v) {
		if (v instanceof ADenseArrayVector) {sub(((ADenseArrayVector)v)); return;}
		int length=checkSameLength(v);
		
		for (int i = 0; i < length; i++) {
			data[i] -= v.unsafeGet(i);
		}
	}
	
	@Override
	public double dotProduct(AVector v, Index ix) {
		if ((v instanceof Vector)) return dotProduct((Vector)v,ix);
		int vl=v.length();
		assert(v.length()==ix.length());
		double result=0.0;
		int[] idata=ix.getData();
		for (int i=0; i<vl; i++) {
			result+=data[idata[i]]*v.unsafeGet(i);
		}
		return result;
	}
	
	public double dotProduct(Vector v, Index ix) {
		int vl=v.length();
		assert(v.length()==ix.length());
		double result=0.0;
		int[] idata=ix.getData();
		for (int i=0; i<vl; i++) {
			result+=data[idata[i]]*v.data[i];
		}
		return result;
	}
	
	@Override
	public double dotProduct(double[] data, int offset) {
		int len=length();
		double result=0.0;
		for (int i=0; i<len; i++) {
			result+=this.data[i]*data[offset+i];
		}
		return result;
	}
	
	@Override
	public double dotProduct(AVector v) {
		checkSameLength(v);
		return v.dotProduct(data,0);
	}
	
	@Override
	public Vector innerProduct(double a) {
		int n=length;
		double[] result=new double[n];
		for (int i=0; i<length; i++) {
			result[i]=data[i]*a;
		}
		return wrap(result);
	}
	
	@Override public Scalar innerProduct(AVector v) {
		return Scalar.create(dotProduct(v));
	}
	
	@Override public Scalar innerProduct(Vector v) {
		return Scalar.create(dotProduct(v));
	}
	
	public double dotProduct(Vector v) {
		int len=checkSameLength(v);

		double result=0.0;
		for (int i=0; i<len; i++) {
			result+=data[i]*v.data[i];
		}
		return result;
	}
	
	public double distanceSquared(Vector v) {
		int len=length();
		double total=0.0;
		for (int i=0; i<len; i++) {
			double d=data[i]-v.data[i];
			total+=d*d;
		}
		return total;
	}
	
	public double distance(Vector v) {
		return Math.sqrt(distanceSquared(v));
	}
	
	public double distance(AVector v) {
		if (v instanceof Vector) {
			return distance((Vector)v);
		}
		return super.distance(v);
	}
	
	public void sub(ADenseArrayVector v) {
		sub(v,0);
	}
	
	public void sub(ADenseArrayVector src,int srcOffset) {
		int length=length();
		assert(length==src.length());
		double[] srcData=src.getArray();
		int voffset=src.getArrayOffset()+srcOffset;
		for (int i = 0; i < length; i++) {
			data[i] -= srcData[voffset + i];
		}
	}
	
	@Override
	public void addMultiple(AVector v, double factor) {
		if (v instanceof ADenseArrayVector) {addMultiple(((ADenseArrayVector)v),factor); return;}
		v.addMultipleToArray(factor, 0, data, 0, length());
	}
	
	@Override
	public void addWeighted(AVector v, double factor) {
		if (v instanceof ADenseArrayVector) {addWeighted(((ADenseArrayVector)v),factor); return;}
		int length=checkSameLength(v);

		for (int i = 0; i < length; i++) {
			data[i] = (data[i]*(1.0-factor)) + (v.unsafeGet(i)*factor);
		}
	}
	
	public void addWeighted(ADenseArrayVector v, double factor) {
		int length=length();
		assert(length==v.length());
		double[] arr=v.getArray();
		int offset=v.getArrayOffset();
		for (int i = 0; i < length; i++) {
			data[i] = (data[i]*(1.0-factor)) + (arr[i+offset]*factor);
		}
	}
	
	@Override
	public void addMultiple(Vector source, Index index, double factor) {
		int len=source.length();
		if (index.length()!=len) throw new VectorzException(ErrorMessages.incompatibleShapes(index, source));
		for (int i=0; i<len; i++) {
			int j=index.data[i];
			this.data[j]+=source.data[i]*factor;
		}
	}
	
	@Override
	public void addMultiple(Index destToSource, Vector source, double factor) {
		int len=this.length();
		if (destToSource.length()!=len) throw new VectorzException("Index length must match this vector");
		for (int i=0; i<len; i++) {
			int j=destToSource.data[i];
			this.data[i]+=source.data[j]*factor;
		}
	}
	
	@Override
	public void multiply(double factor) {
		DoubleArrays.multiply(data, factor);
	}
	
	@Override
	public void multiply(AVector v) {
		if (v instanceof Vector) {multiply(((Vector)v)); return;}
		checkSameLength(v);
		v.multiplyTo(data, 0);	
	}
	
	public void multiply(Vector v) {
		checkSameLength(v);
		DoubleArrays.multiply(data, v.data);
	}
	
	@Override
	public void divide(AVector v) {
		if (v instanceof Vector) {divide(((Vector)v)); return;}
		checkSameLength(v);
		v.divideTo(data, 0);	
	}
	
	public void divide(Vector v) {
		int len=checkSameLength(v);
		for (int i = 0; i < len; i++) {
			data[i]=(data[i]/v.data[i]);
		}	
	}

	@Override
	public boolean isView() {
		return false;
	}
	
	@Override
	public Vector clone() {
		return Vector.wrap(DoubleArrays.copyOf(data));
	}
	
	@Override
	public int hashCode() {
		int hashCode = 1;
		int len=length();
		for (int i = 0; i < len; i++) {
			hashCode = 31 * hashCode + (Hash.hashCode(data[i]));
		}
		return hashCode;
	}
	
	@Override
	public Vector ensureMutable() {
		return this;
	}
	
	@Override 
	public Vector exactClone() {
		return clone();
	}
	
	@Override
	public boolean isPackedArray() {
		return true;
	}
	
	@Override
	public void toDoubleBuffer(DoubleBuffer dest) {
		dest.put(data);
	}
	
	@Override
	public Vector toNormal() {
		Vector v= Vector.create(this);
		v.normalise();
		return v;
	}
	
	@Override
	public Vector toVector() {
		return this;
	}
	
	@Override
	public double[] asDoubleArray() {
		return data;
	}
	
	@Override
	public Vector dense() {
		return this;
	}
	
	@Override
	public boolean equals(AVector v) {
		if (v.length()!=length) return false;
		return v.equalsArray(data,0);
	}
	
	@Override
	public boolean equals(ADenseArrayVector v) {
		if (length!=v.length()) return false;
		return v.equalsArray(data, 0);
	}
	
	@Override
	public boolean equalsArray(double[] arr, int offset) {
		return DoubleArrays.equals(data, 0, arr, offset, length);
	}
	
	@Override
	public boolean equalsArray(double[] arr) {
		return DoubleArrays.equals(data, arr, length);
	}

	@Override
	protected int index(int i) {
		return i;
	}

}
