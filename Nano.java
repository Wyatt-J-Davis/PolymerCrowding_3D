package org.opensourcephysics.sip.CPM;

public class Nano extends Particle{
	private static double tolerance;
	private static double default_r; // AD: overridden by setDefault_r !
	private boolean intersect;

// AD: Note that first constructor is not used!
	public Nano(double x_, double y_, double z_, double tolerance_, double d){
		setX(x_);
		setY(y_);
		setZ(z_);
		setR(d);
		setTolerance(tolerance_);
	}
	
	public Nano(double x_, double y_, double z_){
		setX(x_);
		setY(y_);
		setZ(z_);
		setR(Nano.default_r);
	}
	
	@Override
	public boolean overlap(Particle particle) {
		if(particle instanceof Nano){
			return squaredSeparation(particle) < Math.pow(2.*getrX(), 2); // overlap if center-to-center distance is less than particle diameter
		} else if(this.equals(particle)){
			return false;
		} else{
			Polymer poly = (Polymer) particle;
			return poly.overlap(this);
		}
	}
	
	public void setR(double r_){
		setrX(r_);
		setrY(r_);
		setrZ(r_);
	}
	
	public double getR(){
		return getrX();
	}

	public void move(){
		move(getTolerance());
	}
	
	public boolean getIntersect(){
		return intersect;
	}
	
	public void setIntersect(boolean intersect_){
		intersect = intersect_;
	}

	public static double getTolerance() {
		return tolerance;
	}

	public static void setTolerance(double tolerance) {
		Nano.tolerance = tolerance;
	}

	public static double getDefault_r() {
		return default_r;
	}

	public static void setDefault_r(double defaultD) {
		Nano.default_r = defaultD;
	}
	
	
	public String toPovray(){
		String pString = String.format("sphere{\n"
				+ "\t<%.3f, %.3f, %.3f>, 1\n"
				+ "\tscale <%.3f, %.3f, %.3f>\n"
				+ "\ttexture {\n"
				+ "\t		pigment { color Black }\n"
				+ "\t}\n"
				+ "}", 
				this.getX(), this.getY(), this.getZ(),
				this.getrX(), this.getrY(), this.getrZ());
		return pString;
	}
}
