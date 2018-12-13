package main.java;

//https://msi.nga.mil/msisitecontent/staticfiles/calculators/degree.html
//inspect degree.html
public class WGS84CoordinateUtils {
	
	//For a posssible future Haversine implementation
	protected static double m1 = 111132.92;
	protected static double m2 = -559.82;
	protected static double m3 = 1.175;
	protected static double m4 = -0.0023;
	
	protected static double p1 = 111412.84;
	protected static double p2 = -93.5;
	protected static double p3 = 0.118;
	
	//WGS defining parameters found here: https://en.wikipedia.org/wiki/Geodetic_datum#World_Geodetic_System_1984_(WGS_84)
	
	//semi-major axis
	protected static final double WGSa = 6378137.0;
	//semi-minor axis
	protected static final double WGSb = 6356752.314245;
	//reciprocal of flattening
	protected static final double WGSf = 1 / 298.257223563; 
	
	private final static double TwoPi = 2.0 * Math.PI;
	
	public static double degToRadian(double degrees) {
		return Math.toRadians(degrees);
	}
	public static double radianToDeg(double radians) {
		return Math.toDegrees(radians);
	}
	
	/**
	 * @see <a href="https://www.movable-type.co.uk/scripts/latlong-vincenty.html">vincenty reference formula</a>
	 * @param coord - The coordinate from which to extend a distance given in metres
	 * @param distance - The distance in metres to extend from coord
	 * @param initialBearing - The initial azimuth (In land navigation, azimuth is usually denoted alpha and defined as a 
	 * horizontal angle measured clockwise from a north base line or meridian). i.e. in this case azimuth is the angle 
	 * between the start point and the desired latitude
	 * @return A WGS84 coordinate which is distance m away from coord, in the specified bearing
	 * @throws Exception 
	 */
	private static WGS84Coordinate vincentyGeodesicDirect(WGS84Coordinate coord, double distance, double initialBearing) throws Exception {
	    double Phi1 = Math.toRadians(coord.getLat().doubleValue()), Lambda1 = Math.toRadians(coord.getLng().doubleValue());
	    double Alpha1 = Math.toRadians(initialBearing);
	    double s = distance;

		double a = WGSa;
		double b = WGSb;
		double f = WGSf;

		double sinAlpha1 = Math.sin(Alpha1);
		double cosAlpha1 = Math.cos(Alpha1);

		double tanU1 = (1-f) * Math.tan(Phi1), cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1)), sinU1 = tanU1 * cosU1;
		double Sigma1 = Math.atan2(tanU1, cosAlpha1);
		double sinAlpha = cosU1 * sinAlpha1;
		double cosSqAlpha = 1 - sinAlpha*sinAlpha;
		double uSq = cosSqAlpha * (a*a - b*b) / (b*b);
		double A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
		double B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

	    double  cos2SigmaM, sinSigma, cosSigma, DeltaSigma;

	    double Sigma = s / (b*A), SigmaPrime, iterations = 0;
	    do {
	        cos2SigmaM = Math.cos(2*Sigma1 + Sigma);
	        sinSigma = Math.sin(Sigma);
	        cosSigma = Math.cos(Sigma);
	        DeltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
	            B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));
	        SigmaPrime = Sigma;
	        Sigma = s / (b*A) + DeltaSigma;
	    } while (Math.abs(Sigma-SigmaPrime) > 1e-12 && ++iterations<100);
	    if (iterations >= 100) throw new Exception("Formula failed to converge"); // not possible!

	    double x = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
	    double Phi2 = Math.atan2(sinU1*cosSigma + cosU1*sinSigma*cosAlpha1, (1-f)*Math.sqrt(sinAlpha*sinAlpha + x*x));
	    double Lambda = Math.atan2(sinSigma*sinAlpha1, cosU1*cosSigma - sinU1*sinSigma*cosAlpha1);
	    double C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
	    double L = Lambda - (1-C) * f * sinAlpha *
	        (Sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
	    double Lambda2 = (Lambda1+L+3*Math.PI)%(2*Math.PI) - Math.PI;  // normalise to -180..+180

	    double Alpha2 = Math.atan2(sinAlpha, -x);
	    Alpha2 = (Alpha2 + 2*Math.PI) % (2*Math.PI); // normalise to 0..360

	    return new WGS84Coordinate(Math.toDegrees(Phi2), Math.toDegrees(Lambda2));
//	    return {
//	        point:        new LatLon(Phi2.toDegrees(), Lambda2.toDegrees(), this.datum),
//	        finalBearing: Alpha2.toDegrees(),
//	        iterations:   iterations,
//	    }
	}
	
	//private 
	private static double vincentyGeodesicInverse(WGS84Coordinate from, WGS84Coordinate to) throws Exception {
		
//		double a = WGSa;
//		double b = WGSb;
//		double f = WGSf;
		
		double a = 	WGSa;
		double b = WGSb;
		double f = WGSf;
		
		double phiOne = Math.toRadians(from.getLat().doubleValue());
		double phiTwo = Math.toRadians(to.getLat().doubleValue());

		double lambdaOne = Math.toRadians(from.getLng().doubleValue());
		double lambdaTwo = Math.toRadians(to.getLng().doubleValue());
		
		double L = lambdaTwo - lambdaOne;
		
		double tanU1 = (1-f) * Math.tan(phiOne);
		double cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1));
		double sinU1 = tanU1 * cosU1;
		
		double tanU2 = (1-f) * Math.tan(phiTwo);
		double cosU2 = 1 / Math.sqrt((1 + tanU2*tanU2));
		double sinU2 = tanU2 * cosU2;
		
		double sinLambda;
		double cosLambda;
		double sinSqSigma;
		double sinSigma=0;
		double cosSigma=0;
		double Sigma=0;
		double sinAlpha;
		double cosSqAlpha=0;
		double cos2SigmaM=0;
		double C;
		

		double Lambda = L, LambdaPrime, iterations = 0;
		//System.out.println("Lambda: " + L);
		boolean antimeridian = Math.abs(L) > Math.PI;
		//System.out.println("anitmeridian: " + antimeridian);
	    do {
	        sinLambda = Math.sin(Lambda);
	        cosLambda = Math.cos(Lambda);
	        sinSqSigma = (cosU2*sinLambda) * (cosU2*sinLambda) + (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda);
	        if (sinSqSigma == 0) break; // co-incident points
	        sinSigma = Math.sqrt(sinSqSigma);
	        cosSigma = sinU1*sinU2 + cosU1*cosU2*cosLambda;
	        Sigma = Math.atan2(sinSigma, cosSigma);
	        sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
	        cosSqAlpha = 1 - sinAlpha*sinAlpha;
	        cos2SigmaM = (cosSqAlpha != 0) ? (cosSigma - 2*sinU1*sinU2/cosSqAlpha) : 0; // equatorial line: cosSqAlpha=0 (§6)
	        C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
	        LambdaPrime = Lambda;
	        Lambda = L + (1-C) * f * sinAlpha * (Sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
	        double iterationCheck = antimeridian ? Math.abs(Lambda)-Math.PI : Math.abs(Lambda);
	        //System.out.println("iterationCheck: " + iterationCheck);
	        if (iterationCheck > Math.PI) throw new Exception("Lambda > pi after " + iterations + " iterations.");
	    } while (Math.abs(Lambda-LambdaPrime) > 1e-12 && ++iterations<1000);
	    
	    
	    if (iterations >= 1000) throw new Exception("Formula failed to converge");

		double uSq = cosSqAlpha * (a*a - b*b) / (b*b);
		double A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
		double B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
		double DeltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)-
		        B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+4*cos2SigmaM*cos2SigmaM)));

		double s = b*A*(Sigma-DeltaSigma);

//		double Alpha1 = Math.atan2(cosU2*sinLambda,  cosU1*sinU2-sinU1*cosU2*cosLambda);
//		double Alpha2 = Math.atan2(cosU1*sinLambda, -sinU1*cosU2+cosU1*sinU2*cosLambda);
//
//		Alpha1 = (Alpha1 + 2*Math.PI) % (2*Math.PI); // normalise to 0..360
//		    Alpha2 = (Alpha2 + 2*Math.PI) % (2*Math.PI); // normalise to 0..360
		return s;
	}	
	
	//private 
//	private static double vincentyGeodesic(GPSCoordinate start, GPSCoordinate end) {
//		//taken from http://www.gavaghan.org/blog/2007/11/16/java-gps-receivers-and-geocaching-vincentys-formula/
//		//returns the distance in metres between two GPS coordinates defined according to the WGS84 standard
//		// get constants
//		double a = 	6378137.0;
//		double b = 6356752.314245;
//		double f = 298.257223563;
//		// get parameters as radians
//		double phi1 = Math.toRadians(start.getLat().doubleValue());
//		double lambda1 = Math.toRadians(start.getLng().doubleValue());
//		double phi2 = Math.toRadians(end.getLat().doubleValue());
//		double lambda2 = Math.toRadians(end.getLng().doubleValue());
//		// calculations
//		double a2 = a * a;
//		double b2 = b * b;
//		double a2b2b2 = (a2 - b2) / b2;
//		double omega = lambda2 - lambda1;
//		double tanphi1 = Math.tan(phi1);
//		double tanU1 = (1.0 - f) * tanphi1;
//		double U1 = Math.atan(tanU1);
//		double sinU1 = Math.sin(U1);
//		double cosU1 = Math.cos(U1);
//		double tanphi2 = Math.tan(phi2);
//		double tanU2 = (1.0 - f) * tanphi2;
//		double U2 = Math.atan(tanU2);
//		double sinU2 = Math.sin(U2);
//		double cosU2 = Math.cos(U2);
//		double sinU1sinU2 = sinU1 * sinU2;
//		double cosU1sinU2 = cosU1 * sinU2;
//		double sinU1cosU2 = sinU1 * cosU2;
//		double cosU1cosU2 = cosU1 * cosU2;
//		// eq. 13
//		double lambda = omega;
//		// intermediates we’ll need to compute ‘s’
//		double A = 0.0;
//		double B = 0.0;
//		double sigma = 0.0;
//		double deltasigma = 0.0;
//		double lambda0;
//		boolean converged = false;
//		
//		for (int i = 0; i < 20; i++)
//		{
//			lambda0 = lambda;
//			double sinlambda = Math.sin(lambda);
//			double coslambda = Math.cos(lambda);
//			// eq. 14
//			double sin2sigma = (cosU2 * sinlambda * cosU2 * sinlambda)
//			+ (cosU1sinU2 - sinU1cosU2 * coslambda)
//			* (cosU1sinU2 - sinU1cosU2 * coslambda);
//			double sinsigma = Math.sqrt(sin2sigma);
//			// eq. 15
//			double cossigma = sinU1sinU2 + (cosU1cosU2 * coslambda);
//			// eq. 16
//			sigma = Math.atan2(sinsigma, cossigma);
//			// eq. 17 Careful! sin2sigma might be almost 0!
//			double sinalpha = (sin2sigma == 0) ? 0.0
//			: cosU1cosU2 * sinlambda / sinsigma;
//			double alpha = Math.asin(sinalpha);
//			double cosalpha = Math.cos(alpha);
//			double cos2alpha = cosalpha * cosalpha;
//			// eq. 18 Careful! cos2alpha might be almost 0!
//			double cos2sigmam = cos2alpha == 0.0 ? 0.0
//			: cossigma - 2 * sinU1sinU2 / cos2alpha;
//			double u2 = cos2alpha * a2b2b2;
//			double cos2sigmam2 = cos2sigmam * cos2sigmam;
//			// eq. 3
//			A = 1.0 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
//			// eq. 4
//			B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
//			// eq. 6
//			deltasigma = B * sinsigma * (cos2sigmam + B / 4
//			* (cossigma * (-1 + 2 * cos2sigmam2) - B / 6 * cos2sigmam
//			* (-3 + 4 * sin2sigma) * (-3 + 4 * cos2sigmam2)));
//			// eq. 10
//			double C = f / 16 * cos2alpha * (4 + f * (4 - 3 * cos2alpha));
//			// eq. 11 (modified)
//			lambda = omega + (1 - C) * f * sinalpha
//			* (sigma + C * sinsigma * (cos2sigmam + C * cossigma * (-1 + 2 *
//			cos2sigmam2)));
//			
//			// see how much improvement we got
//			double change = Math.abs((lambda - lambda0) / lambda);
//			if ((i > 1) && (change < 0.0000000000001)){
//				converged = true;
//				break;
//			}
//		}
//		// eq. 19
//		double s = b * A * (sigma - deltasigma);
//		double alpha1;
//		double alpha2;
//		
//		// didn’t converge? must be N/S
//		if (!converged){
//			if (phi1 > phi2){
//				alpha1 = 180.0;
//				alpha2 = 0.0;
//			}
//			else if (phi1 < phi2){
//				alpha1 = 0.0;
//				alpha2 = 180.0;
//			}
//			else{
//				alpha1 = Double.NaN;
//				alpha2 = Double.NaN;
//			}
//		}
//		// else, it converged, so do the math
//		else
//		{
//			double radians;
//			// eq. 20
//			radians = Math.atan2(cosU2 * Math.sin(lambda),
//			(cosU1sinU2 - sinU1cosU2 * Math.cos(lambda)));
//			if (radians < 0.0)
//				radians += TwoPi;
//			alpha1 = Math.toDegrees(radians);
//			// eq. 21
//			radians = Math.atan2(cosU1 * Math.sin(lambda),
//			(-sinU1cosU2 + cosU1sinU2 * Math.cos(lambda))) + Math.PI;
//			if (radians < 0.0)
//				radians += TwoPi;
//				alpha2 = Math.toDegrees(radians);
//		}
//		if (alpha1 >= 360.0)
//			alpha1 -= 360.0;
//		if (alpha2 >= 360.0) 
//			alpha2 -= 360.0;
//		return s;
//		
//	}
	/**
	 * 
	 * @param from
	 * @param to
	 * @return The distance in metres from WGS84 coordinate 'from' to WGS84 coordinate 'to' 
	 * @throws Exception
	 */
	
	public static double getDistanceMetresBetweenWGS84(WGS84Coordinate from, WGS84Coordinate to) throws Exception {
		//returns the length of the geodesic from coordinate from to coordinate to
		return vincentyGeodesicInverse(from, to);
	}
	
	/**
	 * 
	 * @param from
	 * @param to
	 * @return The latitudinal distance in metres from WGS84 coordinate 'from' to WGS84 coordinate 'to' 
	 * @throws Exception
	 */
	public static double getDistanceMetresLatToOther(WGS84Coordinate from, WGS84Coordinate to) throws Exception {
		//returns the length of the geodesic along a line of latitude of the projection of coordinate from to coordinate to (project to latitude space)
		return vincentyGeodesicInverse(from, new WGS84Coordinate(to.getLat(), from.getLng()));
	}
	
	/**
	 * 
	 * @param from
	 * @param to
	 * @return The longitudinal distance in metres from WGS84 coordinate 'from' to WGS84 coordinate 'to' 
	 * @throws Exception
	 */
	public static double getDistanceMetresLngToOther(WGS84Coordinate from, WGS84Coordinate to) throws Exception {
		//returns the length of the geodesic along a line of longitude of the projection of coordinate from to coordinate to (project to longitude space)
		return vincentyGeodesicInverse(from, new WGS84Coordinate(from.getLat(), to.getLng()));
	}
	
	/**
	 * @param from - The point from which to calculate the final GPS coordinate the given distance away
	 * @param distance - The distance at which the return point with be from the provided 'from' coordinate
	 * @param initialBearing - The latitude +- from the 'from' coordinate at which the destination coordiante will lie
	 * @return A GPSCoordinate that lies a distance 'distance' from the provided 'from' coordinate 'initialBearing' degrees lat above or below
	 * @throws Exception
	 */
	public static WGS84Coordinate getGPSCoordGivenDistanceInitialBearing(WGS84Coordinate from, double distance, double initialBearing) throws Exception {
		return vincentyGeodesicDirect(from, distance, initialBearing);
	}
	
	
	
	
	
	
//	//private
//	private static double getLenMetresOneDegreeLat(BigDecimal bigDecimal) {
//		bigDecimal = degToRadian(bigDecimal);
//		return (GPSCoordinateUtils.m1 + (GPSCoordinateUtils.m2 * Math.cos(2 * bigDecimal)) + (GPSCoordinateUtils.m3 * Math.cos(4 * bigDecimal)) +
//		(GPSCoordinateUtils.m4 * Math.cos(6.bigDecimal)));
//	}
	
//	//private
//	private static double getLenMetresOneDegreeLong(double lat) {
//		lat = degToRadian(lat);
//		return (GPSCoordinateUtils.p1 * Math.cos(lat)) + (GPSCoordinateUtils.p2 * Math.cos(3 * lat)) +
//		(GPSCoordinateUtils.p3 * Math.cos(5 * lat));
//	}
	
//	public static double convertLatDegreeDifferenceToMetres(double lat1, double lat2) {
//				// longitude calculation term 3
//		//given two degrees of latitude, calculates the distance in metres between them.
//		//assuming lat1 and lat2 are 'close'
//		return Math.abs(lat1 - lat2) * getLenMetresOneDegreeLat(lat1);
//	}
	
//	public static double convertLatDegreeDifferenceToMetresSigned(BigDecimal bigDecimal, BigDecimal bigDecimal2) {
//		// longitude calculation term 3
//	//given two degrees of latitude, calculates the distance in metres between them.
//	//assuming lat1 and lat2 are 'close'
//	return (bigDecimal.subtract(bigDecimal2).multiply(getLenMetresOneDegreeLat(bigDecimal)));
//	}
	
//	public static double convertLongDegreeDifferenceToMetres(double long1, double long2, double lat) {
//		// longitude calculation term 3
//		//given two degrees of latitude, calculates the distance in metres between them.
//		//assuming lat1 and lat2 are 'close'
//		return Math.abs(long1 - long2) * getLenMetresOneDegreeLong(lat);
//	}
	
//	public static double convertLongDegreeDifferenceToMetresSigned(double long1, double long2, double lat) {
//		// longitude calculation term 3
//		//given two degrees of latitude, calculates the distance in metres between them.
//		//assuming lat1 and lat2 are 'close'
//		return (long1 - long2) * getLenMetresOneDegreeLong(lat);
//	}
	
//	public static double convertMetresLatToDegrees(double metres, double lat) {
//		//given a distance in metres at a certain latitude, convert to a degree delta
//		return metres / getLenMetresOneDegreeLat(lat);
//		
//	}
	
//	public static double convertMetresLongToDegrees(double metres, double lat) {
//		return metres/getLenMetresOneDegreeLong(lat);
//	}

//	public static double getAcuteAngle(GPSCoordinate p1, GPSCoordinate p2, GPSCoordinate p3) throws Exception {
	//This doesn't make sense if using WGS84
//		//first check p1!=p2, p2!=p3, p1!=p3
//		
//		//Assuming all lats positive/negative and all longs positive/negative
//		//angle calculated is /_p1 p2 p3
//		//translate p2 back to 0
//		//rotate p3 to x-axis
//		//find angle between vector p1 and x-axis
//		//first get a translator from p2 back to 0
//		GPSCoordinateTranslator t = p2.getTranslatorFromThisTo(new GPSCoordinate(0,0));		
//		p1 = t.translate(p1);
//		p3 = t.translate(p3);
//		
//		//set p3 as the angle closer 0 deg from positive x - axis
//		if(p1.getAngleRelativeToOriginXAxis() < p3.getAngleRelativeToOriginXAxis()) {
//			GPSCoordinate temp = new GPSCoordinate(p1.getLat(), p1.getLng(), p1.getAlt());
//			p1 = new GPSCoordinate(p3.getLat(), p3.getLng(), p3.getAlt());
//					//p3.clone();
//			p3 = new GPSCoordinate(temp.getLat(), temp.getLng(), temp.getAlt());
//					//temp.clone();
//		}
////		System.out.println("p1 Angle relative to originx axis: " + p1.getAngleRelativeToOriginXAxis());
////		System.out.println("p3 Angle relative to originx axis: " + p3.getAngleRelativeToOriginXAxis());
//		double returnAngle = p1.getAngleRelativeToOriginXAxis() - p3.getAngleRelativeToOriginXAxis();
//		if(returnAngle > 180) return (360 - returnAngle);
//		else return returnAngle;
//	}
	
	//ArrayList<Integer> que = new ArrayList<Integer>(); 
	
//	private void breadthfirstTraverse(Integer node)  {    
//	}

}
