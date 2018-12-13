package main.java;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author Elaina Cole
 */

public enum Datum {
	WGS84("WGS84",  6378137.0, 6356752.314, 0.00167922, 0.081819191, 6367449.146,
	        new double[]{0.000837732, 0.000000760853, 0.00000000119765, 
	        0.00000000000242917, 0.00000000000000571182, 0.0000000000000000148, 
	        0.0000000000000000000410769, 0.000000000000000000000119999, 
	        0.000000000000000000000000364726, 0.00000000000000000000000000115606},
	        
	        new double[]{0.000837732, 0.0000000590587, 0.000000000167348, 
	        0.00000000000021648, 0.000000000000000378793, 0.000000000000000000723677,
	        0.00000000000000000000149345, 0.00000000000000000000000325384,
	        0.00000000000000000000000000739137, 0.0000000000000000000000000000173816}
	        );
	
	private String datum;
    private final double equatorialRadius;
    private final double polarRadius;
    private final double meridianRadius;
    private final double flattening3D;
    private final double eccentricity;
    private final double[] alphaSeries;
    private final double[] betaSeries;
    
    private static final String[] datumNameArray = {Datum.WGS84.datum};
    
    
    public static final Set<String> DATUMS = new HashSet<>(Arrays.asList(datumNameArray));
       
    
    private Datum(String datum, double equatorialRadius, double polarRadius,
            double flattening3D, double eccentricity, double meridianRadius, 
            double[] alphaSeries, double[] betaSeries){
        
        this.datum = datum;
        this.equatorialRadius = equatorialRadius;
        this.polarRadius = polarRadius;
        this.meridianRadius = meridianRadius;
        this.flattening3D =  flattening3D;
        this.eccentricity = eccentricity;
        this.alphaSeries = alphaSeries;
        this.betaSeries = betaSeries;
        
    }
    
    /**
     * Get the datum.
     * @return datum
     */
    public String getDatum(){
        return datum;
    }
    
    /**
     * Get the datum's equatorial radius.
     * @return equatorialRadius
     */
    public double getEquatorialRadius(){
        return equatorialRadius;
    }
    
    /**
     * Get the datum's polar radius.
     * @return polarRadius
     */
    public double getPolarRadius(){
        return polarRadius;
    }
    
    /**
     * Get the datum's meridian radius.
     * @return meridianRadius
     */
    public double getMeridianRadius() {
        return meridianRadius;
    }
    
    /**
     * Get the datum's flattening 3D.
     * @return flattening3D
     */
    public double getFlattening3D() {
        return flattening3D;
    }
    
    /**
     * Get the datum's eccentricity.
     * @return eccentricity
     */
    public double getEccentricity() {
        return eccentricity;
    }
    
    /**
     * Get the datum's alpha series numbers calculated using a Krueger series.
     * @return alphaSeries
     */
    public double[] getAlphaSeries(){
        return alphaSeries;
    }
    
    /**
     * Get the datum's beta series numbers calculated using a Krueger series.
     * @return betaSeries
     */
    public double[] getBetaSeries() {
        return betaSeries;
    }
    
    /**
     * See if enumeration contains a specific datum's information or not.
     * @param datumToFind
     * @return boolean
     */
    public static boolean containsDatum(String datumToFind) {
        boolean datumFound = false;
        
        if (DATUMS.contains(datumToFind.toUpperCase()))
            datumFound = true;
        
        return datumFound;
    }
    
    
}