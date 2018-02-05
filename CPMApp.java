fsteppackage org.opensourcephysics.sip.CPM;

import java.awt.Color;
import java.text.DecimalFormat;

import org.opensourcephysics.controls.AbstractSimulation;
import org.opensourcephysics.controls.SimulationControl;
import org.opensourcephysics.display3d.simple3d.ElementEllipsoid;
import org.opensourcephysics.display3d.simple3d.ElementSphere;
import org.opensourcephysics.frames.Display3DFrame;
import org.opensourcephysics.frames.PlotFrame;
import java.util.Arrays;

/**
 * NanoPolyMixApp is a simulation framework for a binary mixture of
 * colloids(nanoparticles) and polymers model.
 * 
 * @author Wei Kang Lim, Alan Denton
 * @version 1.0 11-10-2013
 * 
 */
public class CPMApp extends AbstractSimulation {
	public enum WriteModes{WRITE_NONE,WRITE_SHAPES,WRITE_ROTATIONS,WRITE_RADIAL,WRITE_ALL;};
	CPM np = new CPM();
	Display3DFrame display3d = new Display3DFrame("3D Frame");
	PlotFrame plotframe = new PlotFrame("Monte Carlo Steps",
			"Number of Intersections", "Number of Intersections");
	double totalIntersections = 0;
	long snapshotIntervals = 1;
	double [] zAxis = {0,0,1};
	double [] xAxis = {1,0,0};
	RDF rdf;
	double polar;
	double azimuth;
	double volFraction;
	int dataPoints;
	int maxDataPoints;
	int runs;
	int currentRun =  1;
	int i =0;
	long timeStarted = 0;
	long timeElapsed = 0;
	final long START_MCS = 50000;
	long totalMCS = 0;
	String minuteInfo = "";
	WriteModes writeMode;
	ElementSphere nanoSphere[];
	ElementEllipsoid polySphere[];
	boolean added = false;
	boolean penetrationEnergyToggle;
	DataFile [] dataFiles;

	public CPMApp(){
		if(display3d != null) display3d.dispose();
		display3d = new Display3DFrame("3D Frame");

	}
	/**
	 * Initializes the model.
	 */
	public void initialize() {
		i++;
		totalIntersections = 0;
		dataPoints = 0;
		np.nP = control.getInt("N Polymers");
		np.nN = control.getInt("N Nano");
		np.q = control.getDouble("Polymer colloid ratio");
		//np.lc = control.getDouble("Lattice constant");
		np.volFraction = control.getDouble("Nanoparticle Volume Fraction");
		np.lc = Math.pow(4 * Math.PI/24/np.volFraction,1./3);
		//np.nano_r = 0.;
		//np.nano_r = control.getDouble("Nanoparticle radius");
		//np.epsilon = 4.40724/Math.pow(np.q,1.29932); // 
		np.epsilon = control.getDouble("Penetration free energy (epsilon):");
		//np.C = control.getDouble("Penetration free parameter:");
		np.nano_r = 0.5;
		np.init_eX = control.getDouble("x");
		np.init_eY = control.getDouble("y");
		np.init_eZ = control.getDouble("z");
		np.pomfRun = false;
		String configuration = control.getString("Initial configuration");
		np.tolerance = control.getDouble("Tolerance");
		np.rotTolerance = control.getDouble("Rotation tolerance");
		np.shapeTolerance = control.getDouble("Shape tolerance");
		np.trialMovesPerMcs = control.getInt("Trial moves per MCS");
		np.energyProfile = control.getBoolean("Energy profile");
		snapshotIntervals = control.getInt("Snapshot interval");
		maxDataPoints = control.getInt("Number of datapoints");
		runs = control.getInt("Number of runs");
		penetrationEnergyToggle =control.getBoolean("Penetration energy");
		Polymer.setExact(control.getBoolean("Exact overlap"));
		totalMCS = START_MCS + snapshotIntervals * maxDataPoints ;
		Polymer.rootFailCount = 0;

		
		switch(control.getInt("Write Mode")){
		case 0: writeMode = WriteModes.WRITE_NONE; break;
		case 1: writeMode = WriteModes.WRITE_SHAPES; break;
		case 2: writeMode = WriteModes.WRITE_ROTATIONS; break;
		case 3: writeMode = WriteModes.WRITE_RADIAL; break;
		case 4: writeMode = WriteModes.WRITE_ALL; break;
		}
		
		np.initialize(configuration, penetrationEnergyToggle);
		
		// add simple3d.Element particles to the arrays 
		if (i == 2) { // particles only allowed to be added once, this is to prevent clone particles when initialize is called twice by the simulation
			nanoSphere = new ElementSphere[np.nN];
			polySphere = new ElementEllipsoid[np.nP];

			for (int i = 0; i < np.nN; i++) {
				nanoSphere[i] = new ElementSphere();
				display3d.addElement(nanoSphere[i]);
			}

			for (int i = 0; i < np.nP; i++) {
				polySphere[i] = new ElementEllipsoid();
				display3d.addElement(polySphere[i]);
			}
			added = true; 
			
			// Initialize visualization elements for nano particles
			for (int i = 0; i < np.nN; i++) {
				nanoSphere[i].setSizeXYZ(2 * np.nanos[i].getrX(),
						2 * np.nanos[i].getrY(), 2 * np.nanos[i].getrZ());
				nanoSphere[i].getStyle().setFillColor(new Color(92, 146, 237, 100));  // light blue with half transparency
				nanoSphere[i].getStyle().setLineColor(new Color(92, 146, 237, 100));
				nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
						np.nanos[i].getZ());
			}

			// Initialize visualization elements for polymer particles
			for (int i = 0; i < np.nP; i++) {
				polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
						2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
				polySphere[i].getStyle().setFillColor(Color.RED);
				polySphere[i].getStyle().setLineColor(Color.RED);
				polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
						np.polymers[i].getZ());
	            if(np.polymers[i].updateRotation()){ 
	                    polySphere[i].setTransformation(np.polymers[i].getRotationTransformation());
	            }
			}

			plotframe.append(0, np.mcs, np.totalIntersectCount);
			display3d.setPreferredMinMax(0, np.Lx, 0, np.Ly, 0, np.Lz);
			display3d.setSquareAspect(true);
		}		
	}
	
	/**
	 * Does a simulation step.
	 */
	public void doStep() {
		// Initialize files for writing output data during the first step
		// This code is being ran here 
		if(np.mcs == 0 && writeMode != WriteModes.WRITE_NONE){
			DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
			DecimalFormat threeDecimal = new DecimalFormat("#0.###");
			DecimalFormat fourDecimal = new DecimalFormat("#0.####");
			String phi_n;
			if(np.volFraction < 0.001){
				phi_n = fourDecimal.format(np.volFraction);
			} else {
				phi_n = threeDecimal.format(np.volFraction);
			}

			String configurations = "# Number of Polymers: " + np.nP +
					"\n# Number of Nanoparticles: "+np.nN +
					"\n# Move Tolerance: "+threeDecimal.format(np.tolerance)+
					"\n# Shape Change Tolerance: "+threeDecimal.format(np.shapeTolerance)+
					"\n# Nanoparticle Volume Fraction: "+phi_n + 
					"\n# Polymer Colloid Ratio: "+threeDecimal.format(np.q)+
					"\n# Lattice Constant: " +threeDecimal.format(np.lc)+
					"\n# Exact overlap: " + Polymer.getExact() + 			
					"\n# Lattice Length: " +threeDecimal.format(np.Lx)+					
					"\n# Rotation Tolerance: "+threeDecimal.format(np.rotTolerance)+
					"\n# Trial Moves Per Mcs: "+np.trialMovesPerMcs+
					"\n# Snapshot Interval: "+largeDecimal.format(this.snapshotIntervals)+
					"\n# Number of Data Points: " + maxDataPoints +
					"\n# Run Number: " + currentRun + 
					"\n# Penetration Energy On: " + this.penetrationEnergyToggle
					;
			switch(writeMode){
			case WRITE_SHAPES:
				dataFiles = new DataFile[3];
				dataFiles[0] = new DataFile("eX", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				dataFiles[1] = new DataFile("eY", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				dataFiles[2] = new DataFile("eZ", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				break;
			case WRITE_ROTATIONS:
				dataFiles = new DataFile[2];
				dataFiles[0] = new DataFile("polar", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				dataFiles[1] = new DataFile("azimuth", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				break;
			case WRITE_RADIAL:
				dataFiles = new DataFile[1]; 
				dataFiles[0] = new DataFile("radial", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				rdf = new RDF(np.nanos, np.Lx);
				break;
			case WRITE_ALL:
				dataFiles = new DataFile[4];
				dataFiles[0] = new DataFile("eX", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				dataFiles[1] = new DataFile("eY", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				dataFiles[2] = new DataFile("eZ", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				dataFiles[3] = new DataFile("radial", configurations, DataFile.FileIdentifier.SIZE_AND_FRACTION, true);
				rdf = new RDF(np.nanos, np.Lx);
				break;
			default:
				break;
			}
			
			// Write out the initialization data
			for(DataFile df : dataFiles){
				df.write();
			}
			
		}
		
		if(np.mcs == 0){
			timeStarted = System.nanoTime();
		}

		// logical step in the CPM class
		np.step();
		
		// update intersect count and mcs step
		plotframe.append(0, np.mcs, np.totalIntersectCount);
		totalIntersections += np.totalIntersectCount;
		display3d.setMessage("Number of mcs steps: " + np.mcs);

		// Visualization updates
		if(control.getBoolean("Visualization on")){
			// nanoparticle visualization updates
			for (int i = 0; i < np.nN; i++) {
//				if(np.nanos[i].intersectPairs.isEmpty()){
//					nanoSphere[i].getStyle().setFillColor(new Color(92, 146, 237, 100));  // light blue with half transparency
//					nanoSphere[i].getStyle().setLineColor(new Color(92, 146, 237, 100));
//				}else{
//					nanoSphere[i].getStyle().setFillColor(new Color(0, 255, 0, 100));  // light blue with half transparency
//					nanoSphere[i].getStyle().setLineColor(new Color(0, 255, 0, 100));
//				}
				nanoSphere[i].setXYZ(np.nanos[i].getX(), np.nanos[i].getY(),
						np.nanos[i].getZ());
			}
			// polymer visualization updates
			for (int i = 0; i < np.nP; i++) {
				polySphere[i].setXYZ(np.polymers[i].getX(), np.polymers[i].getY(),
						np.polymers[i].getZ());
				polySphere[i].setSizeXYZ(2 * np.polymers[i].getrX(),
						2 * np.polymers[i].getrY(), 2 * np.polymers[i].getrZ());
				if(np.polymers[i].updateRotation()){ // this is done to save computation
					polySphere[i].setTransformation(np.polymers[i].getRotationTransformation());
				}
			}
		}
		
		// writing of data
		if (writeMode != WriteModes.WRITE_NONE && np.mcs % snapshotIntervals == 0) {
			if(dataPoints >= maxDataPoints){
				if(currentRun < runs){
					if(writeMode != WriteModes.WRITE_NONE){
						if(writeMode == WriteModes.WRITE_RADIAL){ // Radial Distribution Function 
							dataFiles[0].record(rdf.distributionData());
							System.out.println(rdf.nrData());
							dataFiles[0].write();
						}
						
						if(writeMode == WriteModes.WRITE_ALL){ // Radial Distribution Function 
							dataFiles[3].record(rdf.distributionData());
							System.out.println(rdf.nrData());
							dataFiles[3].write();
						}
						
						dataFiles[0].record("Root fails: " + Polymer.rootFailCount);
						dataFiles[0].record("Root fail (% of total possible comparions): " + (Polymer.rootFailCount / np.mcs) / (np.nN*np.nP) * 100);
						dataFiles[0].record("Root fail (per mcs): " + (Polymer.rootFailCount / np.mcs));
						
						for(DataFile df : dataFiles){
							int elapsedMinutes = (int) Math.floor(timeElapsed/(1000*60));
							int elapsedSeconds = (int) Math.round(timeElapsed/1000) % 60;
							String formatTimeElapsed = (elapsedMinutes == 0) ? elapsedSeconds + "s ": elapsedMinutes + "m " + elapsedSeconds + "s"; 
							df.record("# Total simulation time: " + formatTimeElapsed); 
							df.close();
						}
					}
					this.initialize();
					currentRun++;
					return;
				} else {
					control.setAdjustableValue("Save", true);
					this.stopAnimation();
					return;
				}
			} else{
				switch(writeMode){
					case WRITE_SHAPES:
						if(np.mcs > START_MCS){ // hardcoded 
							for(Polymer poly: np.polymers){
//Sort them 						
//MARK:
							double[] E = new double[3];
							E[0] = poly.geteX(); 
							E[1] = poly.geteY(); 
							E[2] = poly.geteZ();
							Arrays.sort(E);
							dataFiles[0].record(String.valueOf(E[2]));
							dataFiles[1].record(String.valueOf(E[1]));
							dataFiles[2].record(String.valueOf(E[0]));
							}
							dataPoints++;
						}
						break;
					case WRITE_ROTATIONS:
					for (Polymer poly : np.polymers) {
						double[] ellipseAxis = poly.getNewAxis();
						polar = ellipseAxis[2];
						azimuth = Math.atan(ellipseAxis[1]/ ellipseAxis[0]);
						dataFiles[0].record(String.valueOf(polar));
						dataFiles[1].record(String.valueOf(azimuth));
					} 
					dataPoints++;
					break;
					case WRITE_RADIAL:
						if(np.mcs > START_MCS){
							rdf.update();
							dataPoints++;
						}
						break;
					case WRITE_ALL:
						if(np.mcs > START_MCS){ // hardcoded 
							for(Polymer poly: np.polymers){
								 	double[] E = new double[3];
									E[0] = poly.geteX();
									E[1] = poly.geteY();
									E[2] = poly.geteZ();
									Arrays.sort(E); 
//Sort 
									dataFiles[0].record(String.valueOf(E[2]));
									dataFiles[1].record(String.valueOf(E[1]));
									dataFiles[2].record(String.valueOf(E[0]));
							}
							rdf.update();
							dataPoints++;
						}
					break;
					default:break;
				}
			}
		}
		
		// write data onto harddisk for every 100 data values
		if(writeMode != WriteModes.WRITE_NONE && np.mcs % (100*snapshotIntervals) == 0){
			for(DataFile df : dataFiles){
				df.write();
			}
		}
		
		// Simulation info
		timeElapsed = (System.nanoTime() - timeStarted)/1000000;
		if(timeElapsed % 1000 == 0){
			control.clearMessages();
			int elapsedMinutes = (int) Math.floor(timeElapsed/(1000*60)) ;
			int elapsedSeconds = (int) Math.round(timeElapsed/1000) % 60;
			String formatTimeElapsed = (elapsedMinutes == 0) ? elapsedSeconds + "s ": elapsedMinutes + "m " + elapsedSeconds + "s"; 
			control.println("Time Elapsed: " + formatTimeElapsed);
			control.println(minuteInfo);
			
			if(timeElapsed / 60000 > 0){
				DecimalFormat largeDecimal = new DecimalFormat("0.##E0");
				double mcsPerMinute = 60000* np.mcs / timeElapsed;
				double timeRemain = (totalMCS - np.mcs) / mcsPerMinute;
				int minutes = (int) Math.floor(timeRemain);
				int seconds = (int) Math.round((timeRemain - minutes) * 60); 
				minuteInfo = "MCS per minute: " + largeDecimal.format(mcsPerMinute)
						+"\nTime remaining: " + minutes + "m " + seconds + "s";
			}
		}

	}

	/**
	 * Resets the colloid polymer mixture model to its default state.
	 */
	public void reset() {
		enableStepsPerDisplay(true);
		control.setValue("N Polymers", 1);
		control.setValue("N Nano", 216);
		control.setValue("Polymer colloid ratio", 5);
		control.setValue("Nanoparticle Volume Fraction",.2);
		//control.setValue("Lattice constant", 1.378234);
		control.setValue("x", 0.01);
		control.setValue("y", 0.01);
		control.setValue("z", 0.01);
		//CHANGE// Calculate epsilon in terms of 
		control.setValue("Penetration free energy (epsilon):", 4.40724/Math.pow(np.q,1.29932)); // change to appropriate value for q and phi_n
		//control.setValue("Penetration free parameter:", 6.8426); // more accurate
		//control.setValue("Penetration free parameter:", 6.839); // less accurate
		//control.setValue("Penetration free parameter:", 3);		
		control.setValue("Tolerance", 0.1);
		control.setValue("Rotation tolerance", 0.1);
		control.setValue("Shape tolerance", 0.001);
		control.setValue("Exact overlap", true);
		control.setValue("Energy profile", false);
		//control.setValue("Energy profile", true);
		control.setValue("Initial configuration", "square");
		control.setValue("Trial moves per MCS", 1);
		control.setAdjustableValue("Visualization on", true);
		control.setValue("Snapshot interval", 1000);
		control.setValue("Number of datapoints", 10000);
		control.setValue("Number of runs", 5);
		control.setValue("Penetration energy", true);
		control.setValue("Write Mode", 4);
		control.setAdjustableValue("Save", false);
		i=0;
		initialize();
	}

	public void stop() {
		double averageIntersections = totalIntersections / np.mcs;
		if(control.getInt("N Polymers")> 0){
		double volSpheres = np.nP * (4./3.) * Math.PI * np.polymers[0].getrX() * np.polymers[0].getrY() * np.polymers[0].getrZ() * np.nN;
		double volFract = volSpheres / (np.Lx * np.Ly * np.Lz);
		control.println("Average no. of Intersections: " + averageIntersections);
		control.println("Expected average no. of intersections with Ep = 0: " + volFract);
		control.println("Root fails: " + Polymer.rootFailCount);
		control.println("Root fail (% of total possible comparions): " + (Polymer.rootFailCount / np.mcs)*np.nN*np.nP * 100);
		control.println("Root fail (per mcs): " + (Polymer.rootFailCount / np.mcs));
		}
/*
		double volPolymers = 0;
		for(Polymer poly : np.polymers){
			volPolymers += (4./3.) * Math.PI * poly.getrX() * poly.getrY() * poly.getrZ();
		}
		double volSpheres = volPolymers * np.nN;
*/
		
		// close streams
		if(control.getBoolean("Save")){
			if(writeMode != WriteModes.WRITE_NONE){
				control.setAdjustableValue("Save", false);
				if(writeMode == WriteModes.WRITE_RADIAL){ // Radial Distribution Function 
					dataFiles[0].record(rdf.distributionData());
					System.out.println(rdf.nrData());
					dataFiles[0].write();
				}
				
				if(writeMode == WriteModes.WRITE_ALL){ // Radial Distribution Function 
					dataFiles[3].record(rdf.distributionData());
					System.out.println(rdf.nrData());
					dataFiles[3].write();
				}
				
				for(DataFile df : dataFiles){
					df.record("#Number of data points: " + dataPoints);
					df.close();
				}
			}
		}
	}
	
	/**
	 * Starts the Java application.
	 * 
	 * @param args
	 *            command line parameters
	 */
	public static void main(String[] args) { // set up animation control
												// structure using this class
		@SuppressWarnings("unused")
		SimulationControl control = SimulationControl.createApp(new CPMApp());
	}
}
