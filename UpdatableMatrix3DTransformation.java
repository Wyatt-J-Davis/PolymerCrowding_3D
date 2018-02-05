package org.opensourcephysics.sip.CPM;

import org.opensourcephysics.numerics.VectorMath;

public class UpdatableMatrix3DTransformation extends org.opensourcephysics.numerics.Matrix3DTransformation implements Cloneable{

	public UpdatableMatrix3DTransformation(double[][] matrix) {
		super(matrix);
		inverseMatrix = new double[3][3];
        inverseMatrix[0][0] = inverseMatrix[1][1] = inverseMatrix[2][2] = 1;
	}
	
	public UpdatableMatrix3DTransformation(double[][] matrix, double[][] inverse){
		super(matrix);
		inverseMatrix = new double[3][3];

	  if(inverse==null) { // identity matrix
        inverseMatrix[0][0] = inverseMatrix[1][1] = inverseMatrix[2][2] = 1;
        return;
      }
      for(int i = 0; i<inverseMatrix.length; i++) { // loop over the rows
        System.arraycopy(inverse[i], 0, inverseMatrix[i], 0, inverseMatrix[i].length);
      }
	}
	
	/**
	 * Updates this object to perform a theta angled-rotation around axis.  
	 * @param theta angle of rotation
	 * @param axis axis to rotate around
	 */
	private void rotate(double theta, double axis []){
	    double x = axis[0], y = axis[1], z = axis[2];
	    double norm = x*x+y*y+z*z;
	    if(norm!=1) { // this usually doesn't happen because of roundoff but is worth a try
	      norm = 1/Math.sqrt(norm);
	      x *= norm;
	      y *= norm;
	      z *= norm;
	    }
	    double c = Math.cos(theta), s = Math.sin(theta);
	    double t = 1-c;
	    // matrix elements not listed are zero
	    matrix[0][0] = t*x*x+c;
	    matrix[0][1] = t*x*y-s*z;
	    matrix[0][2] = t*x*z+s*y;
	    matrix[1][0] = t*x*y+s*z;
	    matrix[1][1] = t*y*y+c;
	    matrix[1][2] = t*y*z-s*x;
	    matrix[2][0] = t*x*z-s*y;
	    matrix[2][1] = t*y*z+s*x;
	    matrix[2][2] = t*z*z+c;
	    // inverse matrix is null but we know what it is so create it
	    double[][] inv = new double[3][3];
	    inverseMatrix = inv;
	    inv[0][0] = matrix[0][0];
	    inv[1][0] = matrix[0][1];
	    inv[2][0] = matrix[0][2];
	    inv[0][1] = matrix[1][0];
	    inv[1][1] = matrix[1][1];
	    inv[2][1] = matrix[1][2];
	    inv[0][2] = matrix[2][0];
	    inv[1][2] = matrix[2][1];
	    inv[2][2] = matrix[2][2];
	}
	
	/*
	 * Updates the rotation matrix internally to perform an allignment transformation of v1 along v2.
	 */
	public void allign(double[] v1, double[] v2){
		if(v1.length != v2.length || v1.length != 3){
			return;
		}
		
		boolean equal = true;
		for(int i=0;i<v1.length;i++){
			if(v1[i]!=v2[i]) equal=false;
		}
		
		if(equal){
			double[][] identity = {{1,0,0},{0,1,0},{0,0,1}};
			matrix = identity;
			return;
		}
		
	    v1 = VectorMath.normalize(v1.clone());
	    v2 = VectorMath.normalize(v2.clone());
	    double theta = Math.acos(VectorMath.dot(v1, v2));
	    double[] axis = VectorMath.cross3D(v1, v2);
	    rotate(theta,axis);
	}
	
	public UpdatableMatrix3DTransformation clone(){
		UpdatableMatrix3DTransformation clone = new UpdatableMatrix3DTransformation(this.matrix, this.inverseMatrix);
		return clone;
	}

}
