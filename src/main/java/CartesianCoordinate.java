package main.java;

import java.math.BigDecimal;

/**
 * A specification of the familiar 3D Cartesian Coordinate system.
 * Assumes that unit length is 1m.
 * 
 * @author 13383861
 *
 */
public interface CartesianCoordinate {
	
	public BigDecimal getX();
	public BigDecimal getY();
	public BigDecimal getZ();
	public BigDecimal getMetresToOther(CartesianCoordinate other);
	
	//must be able to add and subtract coordinates to each other to form a grid
	public CartesianCoordinate add(CartesianCoordinate other);
	public CartesianCoordinate addX(BigDecimal x)  throws Exception;	
	public CartesianCoordinate addY(BigDecimal y)  throws Exception;	
		
	public CartesianCoordinate subtract(CartesianCoordinate other);
	public CartesianCoordinate subtractY(BigDecimal y)  throws Exception;
	public CartesianCoordinate subtractX(BigDecimal x)  throws Exception;
}
