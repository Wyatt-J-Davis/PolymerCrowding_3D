package org.opensourcephysics.sip.CPM;

/**
 * This class represents the Ellipsoid-Sphere overlap characteristic polynomial, 
 * which root gives one of the coordinates of the point on the ellipsoid, that is closest to the sphere.
 * The other coordinates can be found with this given coordinate, using the constraint of the point being on the ellipsoid (1D freedom).
 * <br>
 * @author Wei Kang Lim
 *
 */
public class ESOverlapPolynomial extends org.opensourcephysics.numerics.Polynomial {
	private double l1,l2,l3; // eigenvalues
	private double x0,y0,z0; // coordinate of sphere
	public ESOverlapPolynomial(double [] coefficients){
		super(coefficients);
	}

	/**
	 * Instantiates the function with the required variables.
	 * @param ellips double[] An array containing the ellipsoid radii squared.<br>
	 * 				 [0,1,2] =  The x,y,z ellipsoid radii.
	 * @param sphere double[] An array containing the sphere coordinates. <br>
	 * 				 [0,1,2] = The sphere coordinates, relative to the center of ellipsoid. <br>
	 */
	public ESOverlapPolynomial(double [] ellips, double [] sphere){
		this(new double [7]); // Creates polynomial with all coefficients equal to zero.
		update(ellips, sphere); // Updates the polynomial coefficients.
	}
	
	/**
	 * Updates the polynomial with new parameters.
	 * @param ellips double[] An array containing the ellipsoid radii squared.<br>
	 * 				 [0,1,2] =  The x,y,z ellipsoid radii.
	 * @param sphere double[] An array containing the sphere coordinates. <br>
	 * 				 [0,1,2] = The sphere coordinates, relative to the center of ellipsoid. <br>
	 */
	public void update(double [] ellips, double [] sphere){
		l1 = ellips[0]; l2=ellips[1]; l3=ellips[2];
		x0 = sphere[0]; y0 = sphere[1]; z0=sphere[2];
		//updateCoefficients();

		double a = l2/(l1*l1) * y0*y0;
		double b = l2/l1 - 1;
		double c = l3/(l1*l1) * z0*z0;
		double d = l3/l1 - 1;
		double e = x0*x0/l1;

		coefficients[6]=b*b*d*d*e;
		coefficients[5]=2*b*d*e*(b+d)*x0;
		coefficients[4]=(b*b*c+a*d*d-b*b*d*d+b*b*e+4*b*d*e+d*d*e)*x0*x0;
		coefficients[3]=2*(b*c+a*d-b*b*d-b*d*d+b*e+d*e)*Math.pow(x0,3);
		coefficients[2]=(a-b*b+c-4*b*d-d*d+e)*Math.pow(x0,4);
		coefficients[1]=-2*(b+d)*Math.pow(x0,5);
		coefficients[0]=-Math.pow(x0,6);
	}
	
	/**
	 * Updates the polynomial coefficients.
	 */
	private void updateCoefficients(){
		//double x0 = this.x0;
		//double y0 = this.y0;
		//double z0 = this.z0;
		double a = l2/(l1*l1) * y0*y0;
		double b = l2/l1 - 1;
		double c = l3/(l1*l1) * z0*z0;
		double d = l3/l1 - 1;
		double e = x0*x0/l1;
		//double e = this.x0*this.x0/l1;

		coefficients[6]=b*b*d*d*e;
		coefficients[5]=2*b*d*e*(b+d)*x0;
		coefficients[4]=(b*b*c+a*d*d-b*b*d*d+b*b*e+4*b*d*e+d*d*e)*x0*x0;
		coefficients[3]=2*(b*c+a*d-b*b*d-b*d*d+b*e+d*e)*Math.pow(x0,3);
		coefficients[2]=(a-b*b+c-4*b*d-d*d+e)*Math.pow(x0,4);
		coefficients[1]=-2*(b+d)*Math.pow(x0,5);
		coefficients[0]=-Math.pow(x0,6);
/*
		coefficients[6]=Math.pow(b,2)*Math.pow(d,2)*e;
		coefficients[5]=2*Math.pow(b,2)*d*e*x0+2*b*Math.pow(d,2)*e*x0;
		coefficients[4]=Math.pow(b,2)*c*Math.pow(x0,2)+a*Math.pow(d,2)*Math.pow(x0,2)-Math.pow(b,2)*Math.pow(d,2)*Math.pow(x0,2)+Math.pow(b,2)*e*Math.pow(x0,2)+4*b*d*e*Math.pow(x0,2)+Math.pow(d,2)*e*Math.pow(x0,2);
		coefficients[3]=2*b*c*Math.pow(x0,3)+2*a*d*Math.pow(x0,3)-2*Math.pow(b,2)*d*Math.pow(x0,3)-2*b*Math.pow(d,2)*Math.pow(x0,3)+2*b*e*Math.pow(x0,3)+2*d*e*Math.pow(x0,3);
		coefficients[2]=a*Math.pow(x0,4)-Math.pow(b,2)*Math.pow(x0,4)+c*Math.pow(x0,4)-4*b*d*Math.pow(x0,4)-Math.pow(d,2)*Math.pow(x0,4)+e*Math.pow(x0,4);
		coefficients[1]=-2*b*Math.pow(x0,5)-2*d*Math.pow(x0,5);
		coefficients[0]=-Math.pow(x0,6);
*/
	}
}
