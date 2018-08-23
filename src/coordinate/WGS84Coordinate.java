package coordinate;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;


/**
 * A class which implements the GPSCoordinate interface.
 * Written according to the EPSG:4326 projected coordinate system.
 * Google maps uses the (EPSG:4326) WGS84 projection system.
 * 
 * References: 
 * 1). http://lyzidiamond.com/posts/4326-vs-3857
 * The World Geodetic System of 1984 is the geographic coordinate system (the three-dimensional one)
 * used by GPS to express locations on the earth. WGS84 is the defined coordinate system for GeoJSON,
 * as longitude and latitude in decimal degrees. For the most part, when you describe a lon/lat 
 * coordinate location, it’s based on the EPSG:4326 coordinate system.
 * There is no way to visualize the WGS84 coordinate system on a two-dimensional plane (map), 
 * so most software programs project these coordinates using an equirectangular projection (Plate-Carrée)
 * 
 * The projected Pseudo-Mercator coordinate system takes the WGS84 coordinate system and projects it onto a square. 
 * (This projection is also called Spherical Mercator or Web Mercator.) 
 * But not all of it – the bounds of Pseudo-Mercator are limited to approximately 85.06º North and South latitude. 
 * This projection was first introduced by Google and is used in almost 100% of web maps, but it’s a strange one:
 * the projection uses the WGS84 coordinate system, which uses the WGS84 ellipsoid, but projects the coordinates onto a sphere.
 * 
 * For the most part, web maps rely on data stored with WGS84 coordinates (in some programs this is called “unprojected” data)
 * and then visualize the data using Pseudo-Mercator. 
 * 
 * 
 * 2).http://earth-info.nga.mil/GandG/wgs84/web_mercator/index.html
 * 
 * 3). When generating points along a line of the grid, incrementing the difference between p1 and pn
 * is what occurs in pyproj: 
 * https://github.com/jswhit/pyproj/blob/c5a4df3d2f716a2960af30c0ce865ac52c2372cb/_proj.pyx
 * 
 * 11/06/2018
 * Cannot figure out how to convert from UTM to WGS84 back to UTM without loss, come back to it!
*/
public class WGS84Coordinate extends CoordinateBase implements CartesianCoordinate{
	
	//convert everything to BigDecimal to ensure that operations give 
	//accurate results
	BigDecimal lat;
	BigDecimal lng;
	BigDecimal alt;
	
	public static BigDecimal LOWER_LAT_BOUND = new BigDecimal(-90);
	public static BigDecimal UPPER_LAT_BOUND = new BigDecimal(90);
	public static BigDecimal LOWER_LNG_BOUND = new BigDecimal(-180);
	public static BigDecimal UPPER_LNG_BOUND = new BigDecimal(180);
	//m above sea level
	public static BigDecimal lowerAltBound = new BigDecimal(0);
	//Assuming that the RAV won't be flying above 20Km...
	public static BigDecimal upperAltBound = new BigDecimal(20000);
	
	public WGS84Coordinate() throws Exception {
		this(0, 0);
	}
		
	public WGS84Coordinate(double lat, double lng) throws Exception {
		//default altitude of 100
		this(lat,lng,100.0);
	}
	public WGS84Coordinate(BigDecimal lat, BigDecimal lng) throws Exception {
		this(lat,lng, new BigDecimal(100.0));
	}
	
	
	public WGS84Coordinate(Double lat, Double lng, Double alt) throws Exception {
			this(new BigDecimal(lat), new BigDecimal(lng), new BigDecimal(alt));
	}
	
	public WGS84Coordinate(BigDecimal lat, BigDecimal lng, BigDecimal alt) throws Exception {
		super(lat, lng, alt);
		setLat(lat);
		setLng(lng);
		setAlt(alt);
	}
	
	
	
	@Override
	public boolean equals(Object other) {
		if(other == null) return false;
		if(other == this) return true;	
		
		if(other.getClass() != getClass()) {
			return false;
		}
		else {
			WGS84Coordinate otherCoord = (WGS84Coordinate) other;
			
			//return true if the difference is less than ~ 1m
//			if((getAlt() != null && otherCoord.getAlt() == null) || getAlt() == null && otherCoord.getAlt() != null) {
//				return false;
//			}
//			if(getAlt() == null && otherCoord.getAlt() == null) {
//				return (otherCoord.getLat().subtract(this.getLat())).abs().compareTo(new BigDecimal(Math.pow(10, -8)))<0 && (otherCoord.getLng().subtract(this.getLng()).abs().compareTo(new BigDecimal(Math.pow(10, -8)))<0);
//			}
//			else{
//				return ((otherCoord.getAlt().subtract(this.getAlt())).abs().compareTo(new BigDecimal(Math.pow(10, -8)))<0 && (otherCoord.getLat().subtract(this.getLat()).abs().compareTo(new BigDecimal(Math.pow(10, -8)))<0) && (otherCoord.getLng().subtract(this.getLng()).compareTo(new BigDecimal(Math.pow(10, -8)))<0));
//			}
			if((getAlt() != null && otherCoord.getAlt() == null) || getAlt() == null && otherCoord.getAlt() != null) {
				return false;
			}
			if(getAlt() == null && otherCoord.getAlt() == null) {
				return (getLat().equals(otherCoord.getLat()) && 
						getLng().equals(otherCoord.getLng()));
			}
			else{
				return (getAlt().equals(otherCoord.getAlt()) && 
						getLat().equals(otherCoord.getLat()) && 
						getLng().equals(otherCoord.getLng()));
			}
		}
	}
	
	@Override
	/**
	 * Create a deep copy of GPSCoordinate
	 */
	public WGS84Coordinate clone() {
		try {
			return new WGS84Coordinate(lat, lng, alt);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			//do nothing; if this object has already been created safely there
			//should be no issues creating a copy
			return this;
		}
	}
	
	@Override
	public String toString() {
		if(alt != null) return "" + lat + ", " + lng + ", " + alt + "";
		else return "" + lat + ", " + lng + "";		
	}
	
	
	public double getLatMetresToOther(WGS84Coordinate otherCoord) throws Exception {
		return WGS84CoordinateUtils.getDistanceMetresLatToOther(this, otherCoord);
	}
	
	/**
	 * gets the distance in metres to another coordinate that has the same longitude as this
	 * @param otherCoord
	 * @return
	 * @throws Exception
	 */
	public double getLngMetresToOther(WGS84Coordinate otherCoord) throws Exception {
		return WGS84CoordinateUtils.getDistanceMetresLngToOther(this, otherCoord);
	}
	
	public Exception throwRangeException(String type, double lowerBound, double upperBound) {
		return new Exception(type + " must lie in the range (" + lowerBound + ", " + upperBound + ")");
	}
	
	public double getMetresToOther(WGS84Coordinate otherCoord) throws Exception {
		return WGS84CoordinateUtils.getDistanceMetresBetweenWGS84(this, otherCoord);
	}
		
	public static boolean verifyLat(BigDecimal lat2) {
		return(LOWER_LAT_BOUND.compareTo(lat2)<0 && lat2.compareTo(UPPER_LAT_BOUND)<0);
	}
	public static boolean verifyLat(Double lat2) {
		return(LOWER_LAT_BOUND.compareTo(new BigDecimal(lat2))<0 && new BigDecimal(lat2).compareTo(UPPER_LAT_BOUND)<0);
	}
	
	public static boolean verifyLng(BigDecimal lng) {
		return(LOWER_LNG_BOUND.compareTo(lng)<0 && lng.compareTo(UPPER_LNG_BOUND)<0);
	}
	
	public static boolean verifyLng(double lng) {
		return(LOWER_LNG_BOUND.compareTo(new BigDecimal(lng))<0 && new BigDecimal(lng).compareTo(UPPER_LNG_BOUND)<0);
	}
	
	public static boolean verifyAlt(BigDecimal alt2) {
		if(alt2 == null) {
			//null is valid for alt
			return true;
		}
		else {
			return(lowerAltBound.compareTo(alt2)<0 && alt2.compareTo(upperAltBound)<0);
		}
	}
	public static boolean verifyAlt(Double alt2) {
		if(alt2 == null) {
			//null is valid for alt
			return true;
		}
		else {
			return(lowerAltBound.compareTo(new BigDecimal(alt2))<0 && new BigDecimal(alt2).compareTo(upperAltBound)<0);
		}
	}
	
	public static boolean verifyAlt(Integer altitude) {
		//need == instead of .equals to avoid null pointer
		if(altitude == null) {
			return true;
		}
		else {
			return verifyAlt(new BigDecimal(altitude));
		}
	}
	
	public BigDecimal getLat() {
		return lat;
	}

	public void setLat(BigDecimal lat2) throws Exception {
		if(verifyLat(lat2)) {
			this.lat = lat2;
		}
		else {
			throw throwRangeException("Latitude", LOWER_LAT_BOUND.doubleValue(), UPPER_LAT_BOUND.doubleValue());
		}
	}

	public BigDecimal getLng() {
		return lng;
	}

	public void setLng(BigDecimal lng2) throws Exception {
		if(verifyLng(lng2)) {
			this.lng = lng2;
		}
		else {
			throw throwRangeException("Longitude", LOWER_LNG_BOUND.doubleValue(), UPPER_LNG_BOUND.doubleValue());
		}
	}

	public BigDecimal getAlt() {
		return alt;
	}

	public void setAlt(BigDecimal alt2) throws Exception {
		if(alt2 == null || verifyAlt(alt2)) {
			this.alt = alt2;
		}
		else {
			throw throwRangeException("Altitude", lowerAltBound.doubleValue(), upperAltBound.doubleValue());
		}
	}
	
	public void setAlt(Double alt2) throws Exception {
		if(alt2 == null || verifyAlt(alt2)) {
			this.alt = new BigDecimal(alt2);
		}
		else {
			throw throwRangeException("Altitude", lowerAltBound.doubleValue(), upperAltBound.doubleValue());
		}
	}
	
	public void setAlt(Integer newAlt) throws NumberFormatException, Exception {
		if(newAlt != null) {
			setAlt(new BigDecimal(newAlt));
		}
		else {
			BigDecimal nullBigDecimal = null;
			setAlt(nullBigDecimal);
		}
	}
	
	public WGS84Coordinate add(WGS84Coordinate otherCoord) throws Exception {
		WGS84Coordinate returnCoord;
//		if(getAlt()!=null && otherCoord.getAlt()!=null) {
//			returnCoord = new GPSCoordinate(new BigDecimal(0), new BigDecimal(0), new BigDecimal(0.0));
//			returnCoord.setAlt(getAlt().add(otherCoord.getAlt()));
//		}
//		else {
			returnCoord = new WGS84Coordinate(0, 0);
//		}
		returnCoord.setLat(getLat().add(otherCoord.getLat()));
		returnCoord.setLng(getLng().add(otherCoord.getLng()));
		return returnCoord;
	}
	
	public WGS84Coordinate multiply(BigDecimal number) throws Exception {
		if(number.equals(BigDecimal.ZERO)) {
			return new WGS84Coordinate(0,0);
		}
		WGS84Coordinate returnCoord = clone();
		returnCoord.setLat(getLat().multiply(number));
		returnCoord.setLng(getLng().multiply(number));
//		for(int counter = 1; counter < number; counter++) {
//			returnCoord = returnCoord.add(this);
//		}
		return returnCoord;
	}
	
	public WGS84Coordinate multiply(int number) throws Exception {
		if(number == 0) {
			return new WGS84Coordinate(0,0);
		}
		WGS84Coordinate returnCoord = clone();
		returnCoord.setLat(getLat().multiply(new BigDecimal(number)));
		returnCoord.setLng(getLng().multiply(new BigDecimal(number)));
//		for(int counter = 1; counter < number; counter++) {
//			returnCoord = returnCoord.add(this);
//		}
		return returnCoord;
	}
	
	public WGS84Coordinate divide(int divisor) throws Exception {
		if(divisor == 0) {
			return new WGS84Coordinate(0,0);
		}
		WGS84Coordinate returnCoord = clone();
		returnCoord.setLat(getLat().divide(new BigDecimal(divisor), MathContext.DECIMAL32));
		returnCoord.setLng(getLng().divide(new BigDecimal(divisor), MathContext.DECIMAL32));
		return returnCoord;
	}

	public WGS84Coordinate subtract(WGS84Coordinate otherCoord) throws Exception {
		WGS84Coordinate returnCoord;
//		if(getAlt()!=null && otherCoord.getAlt()!=null) {
//			returnCoord = new GPSCoordinate();
//			returnCoord.setAlt(getAlt().subtract(otherCoord.getAlt()));
//		}	
//		else {
			returnCoord = new WGS84Coordinate(0, 0);
//		}
		returnCoord.setLat(getLat().subtract(otherCoord.getLat()));
		returnCoord.setLng(getLng().subtract(otherCoord.getLng()));
			
		return returnCoord;
	}

	
	public static BigDecimal getLowerLatBound() {
		return LOWER_LAT_BOUND;
	}


	public static BigDecimal getUpperLatBound() {
		return UPPER_LAT_BOUND;
	}


	public static BigDecimal getLowerLngBound() {
		return LOWER_LNG_BOUND;
	}


	public static BigDecimal getUpperLngBound() {
		return UPPER_LNG_BOUND;
	}


	public static BigDecimal getLowerAltBound() {
		return lowerAltBound;
	}


	public static BigDecimal getUpperAltBound() {
		return upperAltBound;
	}

	@Override
	public BigDecimal getX() {
		// TODO Auto-generated method stub
		return getLng();
	}

	@Override
	public BigDecimal getY() {
		// TODO Auto-generated method stub
		return getLat();
	}

	@Override
	public BigDecimal getZ() {
		// TODO Auto-generated method stub
		return getAlt();
	}

	@Override
	public BigDecimal getMetresToOther(CartesianCoordinate other) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public CartesianCoordinate add(CartesianCoordinate other) {
		// TODO Auto-generated method stub
		return this.add(other);
	}

	@Override
	public CartesianCoordinate addX(BigDecimal x) throws Exception {
		// TODO Auto-generated method stub
		return new WGS84Coordinate(getLat(), getLng().add(x), getAlt());
	}
			
	@Override
	public CartesianCoordinate subtractY(BigDecimal y) throws Exception {
		// TODO Auto-generated method stub
		return new WGS84Coordinate(getLat().subtract(y), getLng(), getAlt());
	}
			
	@Override
	public CartesianCoordinate subtractX(BigDecimal x) throws Exception {
		// TODO Auto-generated method stub
		return new WGS84Coordinate(getLat(), getLng().subtract(x), getAlt());
	}
	
	@Override
	public CartesianCoordinate subtract(CartesianCoordinate other) {
		// TODO Auto-generated method stub
		return this.add(other);
	}

	@Override
	public CartesianCoordinate addY(BigDecimal y) throws Exception {
		// TODO Auto-generated method stub
		return new WGS84Coordinate(getLat().add(y), getLng(), getAlt());
	}
			
}
