package org.opensourcephysics.sip.CPM;

import org.opensourcephysics.numerics.PBC;

public class RDF {
	private int [] nr;
	public double [] gr;
	private Nano [] nanos;
	private double side;
	private double deltaR = 0.005;
	private int nshell;
	private double maxR;
	private int count;
	
	public RDF(Nano [] n, double side){
		this.nanos = n;
		this.side = side;
		maxR = side/2;
		nshell = (int) Math.ceil(maxR / deltaR) + 1;
		nr = new int [nshell];
		gr = new double [nshell];
	}
	
	public void update(){
		count++;
		for(Nano nano: nanos){
			for(Nano nano2: nanos){
				if(!nano.equals(nano2)){
					double rx = PBC.separation(Math.abs(nano2.getX() - nano.getX()), side);
					double ry = PBC.separation(Math.abs(nano2.getY() - nano.getY()), side);
					double rz = PBC.separation(Math.abs(nano2.getZ() - nano.getZ()), side);
		            double r2=rx*rx+ry*ry+rz*rz;
		            
		            if(r2 <= maxR*maxR){
		            	double r = Math.sqrt(r2);
		            	if(r < nanos[0].getrX()){
		            		System.out.println("{" + nano2.getX() + "," + nano2.getY() + "," + nano2.getZ() + "}");
		            		System.out.println("{" + nano.getX() + "," + nano.getY() + "," + nano.getZ() + "}");
		            		System.out.println("Separation is : " + r);
		            	}
		            	int n = (int) Math.ceil(r/deltaR); // nth bin of histogram
		            	nr[n]++;
		            }
				}
			}
		}
	}
	
	public double [] calcDistribution(){
		double rho = nanos.length / Math.pow(side,3);
		double normC = 4*Math.PI*rho*Math.pow(deltaR, 3)*nanos.length*count;
		
		for(int n = 1; n < gr.length; n++){
			double r = deltaR*n;
			if(r <= maxR){
				gr[n] = nr[n] / (normC*n*n);
			}
		}
		
		double [] grClone = new double[gr.length];
		System.arraycopy(gr, 0, grClone, 0, gr.length);
		return grClone;
	}
	
	public String distributionData(){
		StringBuffer data = new StringBuffer();
		gr = calcDistribution();
		for(int i =0; i < gr.length; i++){
			data.append(i*deltaR + " " + gr[i] + "\n");
		}
		String finalData = data.toString();
		return finalData;
	}
	
	public String nrData(){
		StringBuffer data = new StringBuffer();
		for(int i =0; i < nr.length; i++){
			data.append(i*deltaR + " " + nr[i] + "\n");
		}
		String finalData = data.toString();
		return finalData;
	}

}
