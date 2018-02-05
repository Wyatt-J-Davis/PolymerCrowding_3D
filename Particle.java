package org.opensourcephysics.sip.CPM;
import java.util.ArrayList;

import org.opensourcephysics.numerics.PBC;


public abstract class Particle {
	private double x;
	private double y;
	private double z;
	private double rX;
	private double rY;
	private double rZ;
	private static double Lx;
	private static double Ly;
	private static double Lz;
	public ArrayList<Particle> intersectPairs = new ArrayList<Particle>();

	
	public double squaredSeparation(Particle particle){
		double dx = PBC.separation(this.x-particle.x, Lx);
		double dy = PBC.separation(this.y-particle.y, Ly);
		double dz = PBC.separation(this.z-particle.z, Lz);
		return dx*dx+dy*dy+dz*dz;
	}
	
	public abstract boolean overlap(Particle particle);
	
	public void move(double tolerance){
		setX(PBC.position(getX() + tolerance*2.*(Math.random()-0.5), getLx()));
		setY(PBC.position(getY() + tolerance*2.*(Math.random()-0.5), getLy()));
		setZ(PBC.position(getZ() + tolerance*2.*(Math.random()-0.5), getLz()));
	}
	
	public static void setBoundaries(double Lx_, double Ly_, double Lz_){
		setLx(Lx_);
		setLy(Ly_);
		setLz(Lz_);
	}


	
	public void setrX(double rX_){
		if(rX_ > 0) rX = rX_;
	}
	
	public void setrY(double rY_){
		if(rY_ > 0) rY = rY_;
	}
	
	public void setrZ(double rZ_){
		if(rZ_ > 0 ) rZ = rZ_;
	}
	
	public double getrX(){
		return rX;
	}
	
	public double getrY(){
		return rY;
	}
	
	public double getrZ(){
		return rZ;
	}
	
	public void setX(double x_){
		x = x_;
	}
	
	public void setY(double y_){
		y = y_;
	}
	
	public void setZ(double z_){
		z = z_;
	}
	
	public double getX(){
		return x;
	}
	
	public double getY(){
		return y;
	}
	
	public double getZ(){
		return z;
	}

	public static double getLx() {
		return Lx;
	}

	public static void setLx(double lx) {
		Lx = lx;
	}

	public static double getLy() {
		return Ly;
	}

	public static void setLy(double ly) {
		Ly = ly;
	}

	public static double getLz() {
		return Lz;
	}

	public static void setLz(double lz) {
		Lz = lz;
	}
}
