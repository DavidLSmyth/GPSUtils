package coordinate;

import java.math.BigDecimal;

/**
 * 
 * @author David Smyth
 * A base class for a given coordinate system. May need adjustment as I find out more about reference coordinate systems of earth
 */
public abstract class CoordinateBase {
	protected BigDecimal lat;
	protected BigDecimal lng; 
	protected BigDecimal alt;
	
	public CoordinateBase(BigDecimal lat, BigDecimal lng) throws Exception {
		this(lat, lng, null);
	}
	
	public CoordinateBase(BigDecimal lat, BigDecimal lng, BigDecimal alt) throws Exception {
		setLat(lat);
		setLng(lng);
		setAlt(alt);
	}
	
	public BigDecimal getLat() {
		return lat;
	}

	public void setLat(BigDecimal lat) throws Exception{
		this.lat = lat;
	}

	public BigDecimal getLng() {
		return lng;
	}

	public void setLng(BigDecimal lng) throws Exception{
		this.lng = lng;
	}

	public BigDecimal getAlt() {
		return alt;
	}

	public void setAlt(BigDecimal alt) throws Exception{
		this.alt = alt;
	}
	//force subclasses to implement this so that user knows which 
	//class they are working with
	public abstract String toString();
}
