package org.opensourcephysics.sip.CPM;

import org.opensourcephysics.numerics.Matrix3DTransformation;
import org.opensourcephysics.numerics.PBC;
import org.opensourcephysics.numerics.VectorMath;
import java.util.Arrays;
//import java.lang.*;

// NOTE: Solvent quality is set in toRadius method!

public class Polymer extends Particle{
	private static boolean exact;
	private static double tolerance;
	private static double shapeTolerance;
	private static double rotTolerance;
	private static double default_eX;
	private static double default_eY;
	private static double default_eZ;
	private static double q;
	private double eX;
	private double eY;
    	private double eZ;
    	private final double originAxis[] = {0,1,0};
    	//private final double originAxis[] = {0,0,1}; // Wrong!
    	private double oldAxis[] = {0,0,1};
    	private double newAxis[] = {0,1,0}; // transformation axis
    	public double [] closestPoint = new double[3];
    	public double [] overlapSphere = new double[3];
    	private ESOverlapPolynomial overlapPolynomial;
    	private UpdatableMatrix3DTransformation transform;
    	private boolean rotationIsDirty = false; // whether transform requires updating.
	public static long  rootFailCount = 0;
    	public boolean debug= false;



	public Polymer(double x_, double y_, double z_, double tolerance_, double eX, double eY, double eZ, double q_){
		setX(x_);
		setY(y_);
		setZ(z_);
		seteX(eX);
		seteY(eY);
		seteZ(eZ);
		tolerance = tolerance_;
		
		double [][] identity = {
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1}
		};
		
		transform = new UpdatableMatrix3DTransformation(identity);
	}
	
	public Polymer(double x_, double y_, double z_){
		setX(x_);
		setY(y_);
		setZ(z_);
		seteX(default_eX);
		seteY(default_eY);
		seteZ(default_eZ);
		
		double [][] identity = {
				{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1}
		};		
		
		transform = new UpdatableMatrix3DTransformation(identity);
	}
	
	/**
	 * Checks if a polymer overlaps a given particle.
	 * @param Particle A particle object.
	 * @return boolean true if an overlap occurs, false otherwise.
	 */
	public boolean overlap(Particle particle) {
		if(particle instanceof Polymer){
			return false;
		} else if(this.equals(particle)){
			return false;
		} else{
			Nano nano = (Nano) particle;
			
			// Check if polymer is spherical. 
			if(this.geteX() == this.geteY() && this.geteX() == this.geteZ()){
				return this.squaredSeparation(nano) < Math.pow(this.getrX() + nano.getrX(), 2); // Calculate the pbc distance directly.
			} else { // Polymer is ellipsoidal				

				if(exact){
					debug = false;
					// Initial filter, independent of orientation of ellipsoid.
					// If nanosphere center is outside sphere of radius equaling nano_radius plus longest ellipsoid radius, reject immediately.
					double maxEllipsRadius = Math.max(Math.max(this.getrX(), this.getrY()), this.getrZ());
					if(this.squaredSeparation(nano) > Math.pow(maxEllipsRadius + nano.getrX(), 2)){ 
						return false;
					}
					
					// If center of ellipsoid is inside nanosphere, accept immediately.
					if(this.squaredSeparation(nano) <= Math.pow(nano.getrX(), 2)){
						return true;
					}
					
					// Coordinates of nanosphere center from center of ellipsoid
					double [] sphereCoord = {PBC.separation(nano.getX()-this.getX(), Particle.getLx()), PBC.separation(nano.getY()-this.getY(), Particle.getLy()), PBC.separation(nano.getZ()-this.getZ(), Particle.getLz())};

					// Orientation of ellipsoid must be taken into account
					Matrix3DTransformation rotationTransformation =  this.getRotationTransformation();
					//Matrix3DTransformation rotationTransformation =  this.getRotationTransformationInternal();
					if(isRotated()){ // Rotations needed
						sphereCoord = rotationTransformation.inverse(sphereCoord);
						//rotationTransformation.inverse(sphereCoord);
					}
					
					// No overlap if nanosphere is obviously outside ellipsoid coating
					if(Math.pow(sphereCoord[0]/(this.getrX()+nano.getrX()),2) > 1){
						return false;
					}
					if(Math.pow(sphereCoord[1]/(this.getrY()+nano.getrX()),2) > 1){
						return false;
					}
					if(Math.pow(sphereCoord[2]/(this.getrZ()+nano.getrX()),2) > 1){
						return false;
					}

					// Overlap if nanosphere is inside "fattened ellipsoid"
					if(Math.pow(sphereCoord[0]/(this.getrX()+nano.getrX()),2) + 
					   Math.pow(sphereCoord[1]/(this.getrY()+nano.getrY()),2) + 
					   Math.pow(sphereCoord[2]/(this.getrZ()+nano.getrZ()),2) < 1){
						return true;
					}
					
					// Exact solution 
					// Nanosphere lies within ellipsoid coating of thickness nano_radius
					double [] ellipsEigen = {Math.pow(this.getrX(),2), Math.pow(this.getrY(),2), Math.pow(this.getrZ(),2)}; // ellipsoid eigenvalues
					// Update overlap polynomial with new ellipsoid shape and sphere coordinates
					if(overlapPolynomial == null){ // If object not yet instantiated, instantiate a new object
						overlapPolynomial = new ESOverlapPolynomial(ellipsEigen, sphereCoord);
					} else { // Object instantiated, reuse it instead of creating new references for performance
						overlapPolynomial.update(ellipsEigen, sphereCoord);
					}
					
					double yRatio = ellipsEigen[1] / ellipsEigen[0]; //  lambda2/lambda1
					double zRatio = ellipsEigen[2] / ellipsEigen[0]; //  lambda3/lambda1
					double [] roots = overlapPolynomial.rootsReal(); // get real roots

					double distSq;
					double distSqMin = 1.e+10;
					double xMin = 1.e+10;
					double yMin = 1.e+10;
					double zMin = 1.e+10;
					int rootCount = 0;

					for(double x : roots){
						rootCount ++;

						double y = yRatio*x*sphereCoord[1]/(sphereCoord[0] + (yRatio-1)*x);
						double z = zRatio*x*sphereCoord[2]/(sphereCoord[0] + (zRatio-1)*x);

						distSq = Math.pow(x-sphereCoord[0],2)+Math.pow(y-sphereCoord[1],2)+Math.pow(z-sphereCoord[2],2);

						// Filter roots outside sphere of radius equaling nano_radius plus longest ellipsoid radius
						if(distSq > Math.pow(maxEllipsRadius + nano.getrX(), 2)){ // AD
							//System.out.println(" distSq = " + distSq);
							continue;
						}

// Find root corresponding to closest point:
						if(distSq < distSqMin) {
							distSqMin = distSq;
							xMin = x;
							yMin = y;
							zMin = z;
						}

					} // end of loop over roots

					double [] closest = {xMin,yMin,zMin};

					boolean exactOverlap = distSqMin < Math.pow(nano.getrX(),2);
					
					// Check whether closest root is on the ellipsoid:
					//double error = xMin*xMin/ellipsEigen[0] + yMin*yMin/ellipsEigen[1] + zMin*zMin/ellipsEigen[2] - 1;
					//if(exactOverlap && Math.abs(error) > 1.0e-06){
						//System.out.println(rootCount + " roots  " + " error = " + error + " distSqMin = " + distSqMin);
						//return false;
					//}

					if(isRotated()){ // Rotations needed
						closest = rotationTransformation.direct(closest);
					}
						
					if(this.geteX() == this.geteY() && this.geteX() == this.geteZ() && exactOverlap && this.squaredSeparation(nano) >= Math.pow(this.getrX() + nano.getrX(), 2)){
						System.out.println("Inconsistency in overlap.");
					}
						
					overlapSphere[0] = sphereCoord[0];
					overlapSphere[1] = sphereCoord[1];
					overlapSphere[2] = sphereCoord[2];
					closestPoint[0] = closest[0] + this.getX();
					closestPoint[1] = closest[1] + this.getY();
					closestPoint[2] = closest[2] + this.getZ();

					if(exactOverlap){
						debug = true;
						return true;
					}

					return false; // no exact overlap

				} else{ // inexact overlap algorithm

					// Initial filter, independent of orientation of ellipsoid:
					// If nanosphere center is outside sphere of radius equaling nano_radius plus longest ellipsoid radius, reject immediately.
					double maxEllipsRadius = Math.max(Math.max(this.getrX(), this.getrY()), this.getrZ());
					if(this.squaredSeparation(nano) > Math.pow(maxEllipsRadius + nano.getrX(), 2)){ 
						return false;
					}
					
					// If center of ellipsoid is inside nanosphere, accept immediately.
					if(this.squaredSeparation(nano) <= Math.pow(nano.getrX(), 2)){
						return true;
					}

					 //Check whether nanosphere is within the fattened ellipsoid
					double [] nanoCoord = {PBC.separation(nano.getX()-this.getX(), Particle.getLx()), PBC.separation(nano.getY()-this.getY(), Particle.getLy()), PBC.separation(nano.getZ()-this.getZ(), Particle.getLz())};
					if(isRotated()){
						Matrix3DTransformation transformation = this.getRotationTransformation();
						//Matrix3DTransformation transformation = this.getRotationTransformationInternal();
						transformation.inverse(nanoCoord);
					}
					
					if(Math.pow(nanoCoord[0]/(this.getrX()+nano.getrX()),2) + 
					   Math.pow(nanoCoord[1]/(this.getrY()+nano.getrY()), 2) + 
					   Math.pow(nanoCoord[2]/(this.getrZ()+nano.getrZ()),2) < 1){
						return true;
					} 
					
					return false;
				}
			} 
		}
	}
	
	/** Performs a shape change
	 * 
	 */
	public void shapeChange(){
		//Tolerances should have ratios 10:3:1 and shape tolerance should be around 0.001

// start of Alan's edits:
                double oldEX = this.geteX();
                double oldEY = this.geteY();
                double oldEZ = this.geteZ();
		double [] tolFac = {10, 3, 1};
//Wyatt
		//double [] E = {oldEX, oldEY, oldEZ};
		double [] E = new double[3];
                double [] Eold = new double[3];
		E[0] = this.geteX(); 
		E[1] = this.geteY(); 
		E[2] = this.geteZ();
		Eold=E;
		Arrays.sort(E);
//Wyatt
		
		
		for(int i = 0; i < 3; i++) {
			//if(E[i] == Math.max(oldEX, oldEY, oldEZ)){
			if(Eold[i] == E[2]){//Wyatt
				tolFac[i] = 10;
                	//}else if(E[i] == Math.min(oldEX, oldEY, oldEZ)){
			}else if(Eold[i] == E[0]){//Wyatt
				tolFac[i] = 1;
                	}else{
                        	tolFac[i] = 3;
			}

                }
 
		double newEX = oldEX + tolFac[0] * shapeTolerance * 2. * (Math.random() - 0.5);
		double newEY = oldEY + tolFac[1] * shapeTolerance * 2. * (Math.random() - 0.5);
		double newEZ = oldEZ + tolFac[2] * shapeTolerance * 2. * (Math.random() - 0.5);

		// Reject trial changes for negative eigenvalue immediately
		if(newEX < 0 || newEY < 0 || newEZ < 0){ // AD: newEY --> newEX
			return;
		}
// end of Alan's edits:

		this.seteX(newEX);
		this.seteY(newEY);
		this.seteZ(newEZ);
	}
	
	/**
	 * Performs a rotation on the polymer.
	 */
	public void rotate(){
		this.setOldAxis(this.getNewAxis());
		double [] a = this.getOldAxis();
		double [] v = new double[3];
		
		// generate randomly oriented vector v 
		for(int i = 0 ; i < v.length; i++){
			v[i] = Math.random() - 0.5;
		}

		// normalize new (randomly oriented) vector v 
		VectorMath.normalize(v);
				
		// addition of the old and new vector 
		// Note: rotTolerance, which replaces rotMagnitude, should be << 1 (e.g. 0.1)
		for(int i = 0; i < v.length; i++){
			a[i] = a[i] + rotTolerance*v[i];
		}

		// normalize result
		VectorMath.normalize(a);
		this.setNewAxis(a);
		
		rotationIsDirty = true;
	}

	/**
	 * Moves the polymer with the tolerance specified by Polymer.tolerance.
	 */
	public void move(){
		super.move(getTolerance());
	}
	
	public static void setTolerance(double tolerance) {
		Polymer.tolerance = tolerance;
	}
	
	public static double getDefault_eX() {
		return default_eX;
	}

	public static void setDefault_eX(double defaulteX_) {
		Polymer.default_eX = defaulteX_;
	}
	public static void setDefault_eY(double defaulteY_) {
		Polymer.default_eY = defaulteY_;
	}
	public static void setDefault_eZ(double defaulteZ_) {
		Polymer.default_eZ = defaulteZ_;
	}
	public void seteX(double eX_){
		eX = eX_;
		setrX(toRadius(eX));
	}
	
	public void seteY(double eY_){
		eY = eY_;
		setrY(toRadius(eY));
	}
	
	public void seteZ(double eZ_){
		eZ = eZ_;
		setrZ(toRadius(eZ));
	}
	public void setOldAxis(double [] axis){
		if(axis.length == 3){
			this.oldAxis = axis;
		}
		rotationIsDirty = true;
	}

	public static void setQ(double q_){
		Polymer.q = q_;
	}
	public void setNewAxis(double newAxis[]) {
		this.newAxis = newAxis;
		rotationIsDirty = true;
	}


  // -------------------------------------
  // Getter methods
  // -------------------------------------
	private Matrix3DTransformation getRotationTransformationInternal(){
		if(rotationIsDirty){
			transform.allign(originAxis,newAxis);
			rotationIsDirty=false;
		}
		
		return transform;
	}
	
	public Matrix3DTransformation getRotationTransformation(){
		if(rotationIsDirty){
			transform.allign(originAxis,newAxis);
			rotationIsDirty=false;
		}
		
		return transform.clone();
	}

	
  public double[] getNewAxis() {
	  	double [] returnAxis = new double[newAxis.length];
	  	for(int i =0; i < newAxis.length; i++){
	  		returnAxis[i] = newAxis[i];
	  	}
		return returnAxis;
	}


	public double geteX(){
		return eX;
	}
	
	public double geteY(){
		return eY;
	}
	
	public double geteZ(){
		return eZ;
	}
	public double[] getOldAxis(){
		double a[] = new double[oldAxis.length];		
		for(int i = 0; i < oldAxis.length; i ++){
			a[i] = oldAxis[i];
		}
		return a;
	}
	public static double getQ(){
		return Polymer.q;
	}
	public static double getDefault_eZ() {
		return default_eZ; // AD: eY --> eZ
	}
	public static double getTolerance() {
		return tolerance;
	}
	public static double getDefault_eY() {
		return default_eY;
	}
	
	/**
	 * Converts radius eigenvalue to radius.
	 * @param ei Eigenvalue
	 * @param v enumeration representing the radius axis
	 * @return The radius
	 */
	public double toRadius(double ei) {
		//return (q/2.) * Math.sqrt(18. * ei); // theta solvent
		return 1.9634 * q * Math.sqrt(ei); // good solvent
	}
	
	public boolean updateRotation(){
		for(int i = 0; i < oldAxis.length; i++){
			if(oldAxis[i] != newAxis[i]) return true;
		}
		return false;
	}
	
	public boolean isRotated(){
		for(int i = 0; i < originAxis.length; i++){
			if(newAxis[i] != originAxis[i]) return true;
		}
		return false; 
	}
	
	public static double getRotTolerance() {
		return rotTolerance;
	}

	public static void setRotTolerance(double rotTolerance) {
		Polymer.rotTolerance = rotTolerance;
	}
	
	public static double getShapeTolerance() {
		return shapeTolerance;
	}

	public static void setShapeTolerance(double shapeTolerance) {
		Polymer.shapeTolerance = shapeTolerance;
	}
	
	public static boolean getExact(){
		return exact;
	}
	
	public static void setExact(boolean exact){
		Polymer.exact = exact;
	}

	public String toString(){
		return "oldAxis: " + oldAxis[0] + ", " + oldAxis[1] + ", " + oldAxis[2] + "\n" +
			   "newAxis: " + newAxis[0] + ", " + newAxis[1] + ", " + newAxis[2];
	}
	
	public String toPovray(){
		double [] oldAxis = this.getOldAxis();
		double [] newAxis = this.getNewAxis();
		String pString = String.format("sphere{\n"
				+ "\t<%.3f, %.3f, %.3f>, 1\n"
				+ "\tscale <%.3f, %.3f, %.3f>\n"
				+ "\ttexture {\n"
				+ "\t		pigment { color Red }\n"
				+ "\t}\n"
				+ "Reorient_Trans(<%.3f, %.3f, %.3f>, <%.3f, %.3f, %.3f>)"
				+ "}", 
				this.getX(), this.getY(), this.getZ(),
				this.getrX(), this.getrY(), this.getrZ(),
				oldAxis[0], oldAxis[1], oldAxis[2],
				newAxis[0], newAxis[1], newAxis[2]);
		return pString;
	}
	
	public String polynomial(){
		return this.overlapPolynomial.toString();
	}
}
