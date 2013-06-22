package mikera.vectorz;

import mikera.vectorz.ops.AFunctionOp;
import mikera.vectorz.ops.Clamp;
import mikera.vectorz.ops.Exp;
import mikera.vectorz.ops.Identity;
import mikera.vectorz.ops.Linear;
import mikera.vectorz.ops.Logistic;
import mikera.vectorz.ops.NormalRBF;
import mikera.vectorz.ops.Power;
import mikera.vectorz.ops.Quadratic;
import mikera.vectorz.ops.ScaledLogistic;
import mikera.vectorz.ops.SoftPlus;
import mikera.vectorz.ops.Sqrt;
import mikera.vectorz.ops.StochasticBinary;
import mikera.vectorz.ops.Tanh;

/**
 * Static unary operator instances and constructor functions.
 * 
 * @author Mike
 *
 */
public final class Ops {
	public static final Op STOCHASTIC_BINARY=StochasticBinary.INSTANCE;
	public static final Op IDENTITY=Identity.INSTANCE;
	public static final Op LINEAR=Identity.INSTANCE;
	public static final Op LOGISTIC=Logistic.INSTANCE;
	public static final Op SCALED_LOGISTIC=ScaledLogistic.INSTANCE;
	public static final Op RECTIFIER=new Clamp(0.0,Double.MAX_VALUE);
	public static final Op STOCHASTIC_LOGISTIC=Op.compose(STOCHASTIC_BINARY,Logistic.INSTANCE);
	public static final Op TANH=Tanh.INSTANCE;
	public static final Op SOFTPLUS=SoftPlus.INSTANCE;
	public static final Op NEGATE=Linear.create(-1.0, 0.0);
	public static final Op SQUARE = Quadratic.create(1.0, 0.0, 0.0);
	public static final Op SQRT = Sqrt.INSTANCE;
	public static final Op CBRT = Power.create(1.0/3.0);
	public static final Op RBF_NORMAL = NormalRBF.INSTANCE;
	public static final Op TO_DEGREES = Linear.create(180.0/Math.PI, 0.0);
	public static final Op TO_RADIANS = Linear.create(Math.PI/180.0, 0.0);
	public static final Op EXP = new Exp();
	
	public static final Op RECIPROCAL = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return 1.0/x;
		}
		
		@Override
		public double derivative(double x) {
			return -1.0/(x*x);
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return -y*y;
		}
		
		@Override public double averageValue() {return 1.0;}
		@Override public boolean hasInverse() {return true;}
		@Override public Op getInverse() {return this;}
		@Override public boolean hasDerivative() {return true;}
	};
	
	public static final Op SIN = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return Math.sin(x);
		}
		
		@Override
		public double derivative(double x) {
			return Math.cos(x);
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return Math.asin(y);
		}
		
		@Override public boolean hasDerivative() {return true;}
		@Override public double minValue() {return -1.0;}
		@Override public double maxValue() {return 1.0;}
		@Override public Op getDerivativeOp() {return COS;}
	};

	public static final Op COS = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return Math.cos(x);
		}
		
		@Override
		public double derivative(double x) {
			return -Math.sin(x);
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return Math.acos(y);
		}
		
		@Override public boolean hasDerivative() {return true;}
		@Override public double minValue() {return -1.0;}
		@Override public double maxValue() {return 1.0;}
		@Override public Op getDerivativeOp() {return Ops.negate(SIN);}
	};
	
	public static final Op TAN = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return Math.tan(x);
		}
		
		@Override
		public double derivative(double x) {
			double sec=1.0/Math.cos(x);
			return sec*sec;
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return derivative(Math.atan(y));
		}
	};
	
	public static final Op ACOS = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return Math.acos(x);
		}
		
		@Override
		public double derivative(double x) {
			return -1.0/Math.sqrt(1.0-x*x);
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return derivative(Math.cos(y));
		}
		
		@Override public boolean hasDerivative() {return true;}
		@Override public boolean hasInverse() {return true;}
		@Override public Op getInverse() {return COS;}
		@Override public double minValue() {return -Math.PI;}
		@Override public double maxValue() {return Math.PI;}
		@Override public double minDomain() {return -1.0;}
		@Override public double maxDomain() {return 1.0;}
	};
	
	public static final Op ASIN = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return Math.asin(x);
		}
		
		@Override
		public double derivative(double x) {
			return 1.0/Math.sqrt(1.0-x*x);
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return derivative(Math.sin(y));
		}
		
		@Override public boolean hasDerivative() {return true;}
		@Override public boolean hasInverse() {return true;}
		@Override public Op getInverse() {return SIN;}
		@Override public double minValue() {return -Math.PI;}
		@Override public double maxValue() {return Math.PI;}
		@Override public double minDomain() {return -1.0;}
		@Override public double maxDomain() {return 1.0;}
	};
	
	public static final Op ATAN = new AFunctionOp() {
		@Override
		public double apply(double x) {
			return Math.atan(x);
		}
		
		@Override
		public double derivative(double x) {
			return 1.0/(1.0+x*x);
		}
		
		@Override
		public double derivativeForOutput(double y) {
			return derivative(Math.tan(y));
		}
		
		@Override public boolean hasDerivative() {return true;}
		@Override public boolean hasInverse() {return true;}
		@Override public Op getInverse() {return TAN;}
		@Override public double minValue() {return -Math.PI;}
		@Override public double maxValue() {return Math.PI;}
	};
	
	public static Op negate(Op op) {
		return NEGATE.compose(op);
	}

	public static final Op compose(Op a, Op b) {
		return a.compose(b);	
	}

	public static Op product(Op a, Op b) {
		return a.product(b);
	}
	
	public static Op sum(Op a, Op b) {
		return a.sum(b);
	}
	
	public static Op divide(Op a, Op b) {
		return a.divide(b);
	}
}
