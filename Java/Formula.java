package rad.sys;

public class Formula
{
	private static final double PI = Math.PI;
	
	private static final double A = 6378.137000;
	private static final double F = 1/298.257223563;
	private static final double E_ = F * (2 - F);
	
	private static final double T_0C = 15;
	
	private static final double T_RCPF = 0.0019812;
	private static final double T_TRC = -56.5;
	
	private static final double P_0H = 29.92126;
	
	private static final double RHO_0SF = 0.002376892;
	
	private static final double CS_0 = 661.4786;
	
	/**
	 * Converts given angular value into radian.
	 * @param degree Angular value in degree
	 * @return Radian value of angular value.
	 */
	public double angleToRadian(double degree)
	{
		return (PI / 180) * degree;
	}
	
	/**
	 * Converts given angular value into degree.
	 * @param radian Angular value in radian
	 * @return Degree value of angular value.
	 */
	public double angleToDegree(double radian)
	{
		return (180 / PI) * radian;
	}
	
	/**
	 * Converts any nautical distance into radian distance.
	 * @param nautical Nautical distance value
	 * @return Radian value of nautical distance.
	 */
	public double distanceToRadian(double nautical)
	{
		return (PI / (60 * 180)) * nautical;
	}
	
	/**
	 * Converts any radian distance into nautical distance.
	 * @param radian Radian distance value.
	 * @return Nautical value of radian distance.
	 */
	public double distanceToNautical(double radian)
	{
		return ((180 * 60) / PI) * radian;
	}
	
	/**
	 * Calculates great circle distance between two coordinate points.
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return Radian distance value between two points.
	 */
	public double twoPointGCD(double lat1, double lon1, double lat2, double lon2)
	{
		return Math.acos(Math.sin(lat1) * Math.sin(lat2) + Math.cos(lat1) * Math.cos(lat2) * Math.cos(lon1 - lon2));
	}
	
	/**
	 * Calculates great circle distance between two coordinate points. <br>
	 * This method is for shorter distance.
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return Radian distance value between two points.
	 */
	public double twoPointGCDShort(double lat1, double lon1, double lat2, double lon2)
	{
		return 2 * Math.asin(Math.sqrt((Math.sin((lat1-lat2)/2) * Math.sin((lat1-lat2)/2)) + Math.cos(lat1) * Math.cos(lat2) * (Math.sin((lon1 - lon2) / 2) * Math.sin((lon1 - lon2) / 2))));
	}
	
	/**
	 * Calculates true course between two points (Initial course bearing, not straight point-to-point course!)
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return Radian angle value of initial course bearing.
	 */
	public double courseBearing(double lat1, double lon1, double lat2, double lon2)
	{
		double d = twoPointGCD(lat1, lon1, lat2, lon2);
		
		if(Math.sin(lon2 - lon1) < 0)
			return Math.acos((Math.sin(lat2) - Math.sin(lat1) * Math.cos(d)) / (Math.sin(d) * Math.cos(lat1)));
		else
			return 2 * PI - Math.acos((Math.sin(lat2) - Math.sin(lat1) * Math.cos(d)) / (Math.sin(d) * Math.cos(lat1)));
	}
	
	/**
	 * Calculates true course between two points, without GCD value (Initial course bearing, not straight point-to-point course!)
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return Radian angle value of initial course bearing.
	 */
	public double courseBearingNoGCD(double lat1, double lon1, double lat2, double lon2)
	{
		return getMod(Math.atan2(Math.sin(lon1 - lon2) * Math.cos(lat2), Math.cos(lat1) * Math.sin(lat2) - Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon1 - lon2)), 2 * PI);
	}
	
	/**
	 * Calculates intermediate latitude position which lies on a great circle distance path.
	 * @param lon Longitude of any given points that lies on a great circle distance
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return Radian angle value of latitude point from given longitude.
	 */
	public double getLatPointOnGCD(double lon, double lat1, double lon1, double lat2, double lon2)
	{
		if(Math.sin(lon1 - lon2) == 0)
			return 0;
		
		return Math.atan((Math.sin(lat1) * Math.cos(lat2) * Math.sin(lon - lon2) - Math.sin(lat2) * Math.cos(lat1) * Math.sin(lon - lon1)) / (Math.cos(lat1) * Math.cos(lat2) * Math.sin(lon1 - lon2)));
	}
	
	/**
	 * Calculates polar coordinate of given points from given center.<br>
	 * Note that this bearing is from center, not from the given target. 
	 * @param lat Latitude of target point. Must be in radian.
	 * @param lon Longitude of target point. Must be in radian.
	 * @param lat1 Latitude of center point. Must be in radian.
	 * @param lon1 Longitude of center point. Must be in radian.
	 * @return An array of double with two elements, {bearing, distance}, both in radian value.
	 */
	public double[] getCoordinatePolar(double lat, double lon, double lat1, double lon1)
	{
		double[] coord = new double[2];
		
		coord[0] = courseBearing(lat1, lon1, lat, lon);
		coord[1] = twoPointGCD(lat1, lon1, lat, lon);
	    
	    return coord;
	}
	
	/**
	 * Calculates Lat/Lon point using given distance and course bearing, and center position.
	 * @param d Distance from center to target. Must be in radian value.
	 * @param tc Course bearing from center to target. Must be in radian value.
	 * @param lat1 Latitude of center point. Must be in radian.
	 * @param lon1 Longitude of center point. Must be in radian.
	 * @return An array of double with two elements, {Latitude, Longitude}, both in radian value.
	 */
	public double[] getCoordinate(double d, double tc, double lat1, double lon1)
	{
		double[] coord = new double[2];
		
		coord[0] = Math.asin(Math.sin(lat1) * Math.cos(d) + Math.cos(lat1) * Math.sin(d) * Math.cos(tc));
	    double dlon = Math.atan2(Math.sin(tc) * Math.sin(d) * Math.cos(lat1), Math.cos(d) - Math.sin(lat1) * Math.sin(coord[0]));
	    coord[1] = getMod(lon1 - dlon + PI, 2 * PI) - PI;
	    
	    return coord;
	}
	
	/**
	 * Calculates third coordinate measured from first and second point.<br>
	 * It uses given two initial course bearing to calculate a thrid point. If these two points does not intersects, returns nothing.
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @param crs13 Course bearing from first point to third point. Must be in radian.
	 * @param crs23 Course bearing from second point to third point. Must be in radian.
	 * @return An array of double with two elements, {Lat, Lon}, both in radian value. Null if there are no intersecting points.
	 */
	public double[] getIntersectRadial(double lat1, double lon1, double lat2, double lon2, double crs13, double crs23)
	{
		double[] coord = new double[2];
		
		double dst12 = twoPointGCD(lat1, lon1, lat2, lon2);
		double crs12 = courseBearing(lat1, lon1, lat2, lon2);
		double crs21 = courseBearing(lat2, lon2, lat1, lon1);
		
		double ang1 = getMod(crs13 - crs12 + PI, 2. * PI) - PI;
		double ang2 = getMod(crs21 - crs23 + PI, 2. * PI) - PI;
		
		if(Math.sin(ang1) == 0 && Math.sin(ang2) == 0)
			return null;
		else if(Math.sin(ang1) * Math.sin(ang2) < 0)
			return null;
		else
		{
			ang1 = Math.abs(ang1);
			ang2 = Math.abs(ang2);
			double ang3 = Math.acos(-Math.cos(ang1) * Math.cos(ang2) + Math.sin(ang1) * Math.sin(ang2) * Math.cos(dst12));
			double dst13 = Math.atan2(Math.sin(dst12) * Math.sin(ang1) * Math.sin(ang2), Math.cos(ang2) + Math.cos(ang1) * Math.cos(ang3));
			coord = getCoordinate(dst13, ang3, lat1, lon1);
			
			return coord;
		}
	}
	
	/**
	 * Calculates two passing longitude coordinate with given parallel latitude coordinate which passes through great circle distance path.
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @param lat3 Latitude of thrist point which crosses GCD and is parallel. Must be in radian.
	 * @return An array of double with two elements, {lon3_0, lon3_1}, both in radian. Returns null if not crossing GCD.
	 */
	public double[] getXParallel(double lat1, double lon1, double lat2, double lon2, double lat3)
	{
		double[] loncoord = new double[2];
		
		double lon12 = lon1 - lon2;
		double A = Math.sin(lat1) * Math.cos(lat2) * Math.cos(lat3) * Math.sin(lon12);
		double B = Math.sin(lat1) * Math.cos(lat2) * Math.cos(lat3) * Math.cos(lon12) - Math.cos(lat1) * Math.sin(lat2) * Math.cos(lat3);
		double C = Math.cos(lat1) * Math.cos(lat2) * Math.sin(lat3) * Math.sin(lon12);
		
		double lon = Math.atan2(B, A);
		
		if(Math.abs(C) > Math.sqrt((A * A) + (B * B)))
			return null;
		else
		{
			double dlon = Math.acos(C / Math.sqrt((A * A) + (B * B)));
			loncoord[0] = getMod(lon1 + dlon + lon + PI, 2 * PI) - PI;
			loncoord[1] = getMod(lon1 - dlon + lon + PI, 2 * PI) - PI;
			
			return loncoord;
		}
	}
	
	/**
	 * Calculates intermediate point that lies on a GCD, with given fractional value.
	 * @param f Fractional value between 0 to 1. 0 is close to the first point, 1 is close to the second point.
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return An array of double with two elements, {Lat, Lon}, both in radian value.
	 */
	public double[] getIntermedPointOnGCD(double f, double lat1, double lon1, double lat2, double lon2)
	{
		double[] coord = new double[2];
		
		double d = twoPointGCD(lat1, lon1, lat2, lon2);
		
		double A = Math.sin((1 - f) * d) / Math.sin(d);
		double B = Math.sin(f * d) / Math.sin(d);
		
		double x = A * Math.cos(lat1) * Math.cos(lon1) + B * Math.cos(lat2) * Math.cos(lon2);
		double y = A * Math.cos(lat1) * Math.sin(lon1) + B * Math.cos(lat2) * Math.sin(lon2);
		double z = A * Math.sin(lat1) + B * Math.sin(lat2);
		
		coord[0] = Math.atan2(z, Math.sqrt((x * x) + (y * y)));
		coord[1] = Math.atan2(y, x);
		
		return coord;
	}
	
	/**
	 * Calculates cross-track error distance.
	 * @param xd Distance between starting point to arrived position. In Radian.
	 * @param tx Course bearing between starting point to arrived position. In Radian.
	 * @param tc Initial course bearing between starting point to actual destination position. In Radian. 
	 * @return Radian value of cross-track error distance.
	 */
	public double getXTD(double xd, double tx, double tc)
	{
		return Math.asin(Math.sin(xd) * Math.sin(tx - tc));
	}
	
	/**
	 * Calculates along-track distance, the distance between starting point to final point to the point right-angle to actual arrival point
	 * @param xd Distance between starting point to arrived potision. In Radian.
	 * @param xtd Cross-track error distance between start point to arrived position. In radian.
	 * @return Radian value of along-track distance.
	 */
	public double getATD(double xd, double xtd)
	{
		return Math.acos(Math.cos(xd) / Math.cos(xtd));
	}
	
	/**
	 * Calculates along-track distance, the distance between starting point to final point to the point right-angle to actual arrival point
	 * @param xd Distance between starting point to arrived potision. In Radian.
	 * @param xtd Cross-track error distance between start point to arrived position. In radian.
	 * @return Radian value of along-track distance.
	 */
	public double getATDShort(double xd, double xtd)
	{
		return Math.asin(Math.sqrt((Math.sin(xd) * Math.sin(xd)) - (Math.sin(xtd) * Math.sin(xtd)) ) / Math.cos(xtd));
	}
	
	/**
	 * Calculates a coordinate with given initial point, target point, and an offset points.<br>
	 * This uses give distance from the offset point and calculates any possible point that lies between great circle distance path and the offset point. 
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @param lat3 Latitude of third point. Must be in radian. (Offset Point)
	 * @param lon3 Longitude of third point. Must be in radian. (Offset Point)
	 * @param d Distance from third point and great circle distance path.
	 * @param dir -1 or 1.
	 * @return An array of double with two elements, {Lat, Lon}, both in radian.
	 */
	public double[] getCoordOffset(double lat1, double lon1, double lat2, double lon2, double lat3, double lon3, double d, double dir)
	{
		double crsA = courseBearing(lat1, lon1, lat3, lon3) - courseBearing(lat1, lon1, lat2, lon2);
		double dst13 = twoPointGCD(lat1, lon1, lat3, lon3);
		double XTD = Math.pow((Math.cos(dst13) * Math.cos(dst13)) + (Math.sin(dst13) * Math.sin(dst13)) * (Math.cos(crsA) * Math.cos(crsA)), 1/2);
		double ATD = Math.atan2(Math.sin(dst13) * Math.cos(crsA), Math.cos(dst13));
		
		if((Math.cos(d) * Math.cos(d)) > (XTD * XTD))
			return null;
		else
		{
			double dst = (ATD * dir) + Math.acos(Math.cos(d) / XTD);
			double tc = courseBearing(lat1, lon1, lat3, lon3);
			double[] coord = getCoordinate(dst, tc, lat1, lon1);
			
			return coord;
		}
	}
	
	/**
	 * Calculates polar coordinate between given position.
	 * @param lat1 Latitude of first point. Must be in radian.
	 * @param lon1 Longitude of first point. Must be in radian.
	 * @param lat2 Latitude of second point. Must be in radian.
	 * @param lon2 Longitude of second point. Must be in radian.
	 * @return An array of double with two elements, {Bearing, Distance}, both in radian.
	 */
	public double[] getApproxFlatPolarCoord(double lat1, double lon1, double lat2, double lon2)
	{
		double[] coord = new double[2];
		
		double R1 = A * (1 - E_) / Math.pow(1 - E_ * (Math.sin(lat1) * Math.sin(lat1)), 3/2);
		double R2 = A / Math.sqrt(1 - E_ * (Math.sin(lat1) * Math.sin(lat1)));
		
		double dstN = R1 * (lat2 - lat1);
		double dstE = R2 * Math.cos(lat1) * (lon2 - lon1);
		
		coord[0] = getMod(Math.atan2(dstE, dstN), 2 * PI);
		coord[1] = Math.sqrt((dstN * dstN) + (dstE * dstE));
		
		return coord;
	}
	
	/**
	 * Calculates wind heading and speed with given course bearing, heading, true air speed and ground speed.
	 * @param crs Course bearing. Must be in radian.
	 * @param hdg Heading. Must be in radian.
	 * @param tas True air speed. In knots.
	 * @param gs Ground speed. In knots.
	 * @return An array of double with two elements, {Bearing, Speed}, in radian and knots 
	 */
	public double[] getWindData(double crs, double hdg, double tas, double gs)
	{
		double[] data = new double[2];
		
		data[0] = crs + Math.atan2(tas * Math.sin(hdg - crs), tas * Math.cos(hdg - crs) - gs);
		data[0] = getMod(data[0], 2 * PI);
		data[1] = Math.sqrt(Math.pow(tas - gs, 2) + 4 * tas * gs * Math.pow(Math.sin((hdg - crs) / 2), 2));
		
		return data;
	}
	
	/**
	 * Calculates heading and ground speed with given course bearing, true air speed, wind speed and wind direction.
	 * @param crs Course bearing. Must be in radian.
	 * @param tas True air spee. In knots.
	 * @param ws Wind speed, in knots.
	 * @param wd Wind direction. Must be in radian. 
	 * @return An array of double with two elements, {Heading, Speed}, in radian and knots.
	 */
	public double[] getObjectHDGGSWind(double crs, double tas, double ws, double wd)
	{
		double swc = (ws / tas) * Math.sin(wd - crs);
		
		if(Math.abs(swc) > 1)
			return null;
		else
		{
			double[] data = new double[2];
			data[0] = crs + Math.asin(swc);
			data[0] = getMod(data[0], 2 * PI);
			
			data[1] = tas * Math.sqrt(1 - (swc * swc)) - ws * Math.cos(wd - crs);
			if(data[1] < 0)
				return null;
			else return data;
		}
	}
	
	/**
	 * Calculates ground speed and course bearing with given heading, true air speed, wind speed and wind direction.
	 * @param hdg Heading. Must be in radian.
	 * @param tas True air speed. In knots.
	 * @param ws Wind speed. In knots.
	 * @param wd Wind direction. Must be in radian.
	 * @return An array of double with two elements, {Ground speed, Course Bearing}, in knots and radian.
	 */
	public double[] getObjectGSCRSWind(double hdg, double tas, double ws, double wd)
	{
		double[] data = new double[2];
		
		data[0] = Math.sqrt((ws * ws) + (tas * tas) - 2 * ws * tas * Math.cos(hdg - wd));
		double wca = Math.atan2(ws * Math.sin(hdg - wd), tas - wd * Math.cos(hdg - wd));
		if(ws < tas)
			wca = Math.asin((ws / data[0]) * Math.sin(hdg - wd));
		data[1] = getMod(hdg + wca, 2 * PI);
		
		return data;
	}
	
	/**
	 * Calculates head-wind and cross-wind.
	 * @param ws Wind speed. In knots.
	 * @param wd Wind direction. Must be in radian.
	 * @param hdg Heading. Must be in radian.
	 * @return An array of double with two elements, {Head-wind, Cross-wind}, both in knots.
	 */
	public double[] getHWXW(double ws, double wd, double hdg)
	{
		double[] data = new double[2];
		
		data[0] = ws * Math.cos(wd - hdg);
		data[1] = ws * Math.sin(wd - hdg);
		
		return data;
	}
	
	/**
	 * Calculates true air speed and wind speed with given three positions.
	 * @param v1 Ground speed of first object. In knots.
	 * @param v2 Ground speed of second object. In knots.
	 * @param v3 Ground speed of third object. In knots.
	 * @return An array of double with two elements, {True air speed, Wind speed}, both in knots
	 */
	public double[] getGPSTASWSData(double v1, double v2, double v3)
	{
		double vms = (Math.pow(v1, 2) + Math.pow(v2, 2) + Math.pow(v3, 2)) / 3;
		double a1 = (v1 * v1)/ (vms - 1);
		double a2 = (v2 * v2)/ (vms - 1);
		double a3 = (v3 * v3)/ (vms - 1);
		double mu = (Math.pow(a1, 2) + Math.pow(a2, 2) + Math.pow(a3, 2)) / 6;
		
		double bp = (1 / 2) + Math.sqrt(1 / (4 - mu));
		double bm = mu / bp;
		
		double[] data = new double[2];
		data[0] = Math.sqrt(vms * bp);
		data[1] = Math.sqrt(vms * bm);
		
		if(data[0] < data[1])
		{
			data[0] = Math.sqrt(vms * bm);
			data[1] = Math.sqrt(vms * bp);
		}
		
		return data;
	}
	
	/**
	 * Calculates temperature with given altitude.
	 * @param h Altitude in feet.
	 * @return Temperature in celsius.
	 */
	public double getTemp(double h)
	{
		if(h < 36089.24)
			return T_0C - (T_RCPF * h);
		else
			return T_TRC;
	}
	
	/**
	 * Calculates pressure at given altitude.
	 * @param h Altitude in feet.
	 * @return Pressure in Hg.
	 */
	public double getPressure(double h)
	{
		if(h < 36089.24)
			return P_0H * Math.pow(1 - 6.8755856E-6 * h, 5.2558797);
		else
		{
			double p_tr = 0.2233609 * P_0H;
			return p_tr * Math.exp(-4.806346E-5 * (h - 36089.24));
		}
	}
	
	/**
	 * Calculates air density at given altitude.
	 * @param h Altitude in feet.
	 * @return Air density in slugs/feet^3
	 */
	public double getRho(double h)
	{
		if(h < 36089.24)
			return RHO_0SF * Math.pow((1. - 6.8755856E-6 * h), 4.2558797);
		else
		{
			double rho_Tr = 0.2970756 * RHO_0SF;
			return rho_Tr * Math.exp(-4.806346E-5 * (h - 36089.24));
		}
	}
	
	/**
	 * Calculates corrected altitude with given pressure and altitude.
	 * @param p Pressure in Hg.
	 * @param h Altitude in feet.
	 * @return Corrected altitude in feet.
	 */
	public double getCorrAlt(double p, double h)
	{
		double corr = 145442.2 * (1 - Math.pow((p / 29.92126),0.190261));
		return h + corr;
	}
	
	/**
	 * Calculates altitude based on air density.
	 * @param ph Altitude corrected by pressure. In feet.
	 * @param th Temperature at given altitude. In celsius.
	 * @param t Actual temperature at given altitude. In celcius.
	 * @return Corrected altitude in feet.
	 */
	public double getRhoAlt(double ph, double th, double t)
	{
		th += 273.15;
		t += 273.15;
		
		return ph + (th / T_RCPF) * (1. - Math.pow((th / t), 0.2349690));
	}
	
	/**
	 * Calculates mach number with given true air speed and speed of sound.
	 * @param tas True air speed, in knots.
	 * @param cs Speed of sound, in knots.
	 * @return Mach number.
	 */
	public double getMachNumber(double tas, double cs)
	{
		return tas / cs;
	}
	
	/**
	 * Calculates speed of sound at given temperature
	 * @param t Temperature in celsius.
	 * @return Speed of sound in knots.
	 */
	public double getCS(double t)
	{
		return 38.967854 * Math.sqrt(t + 273.15);
	}
	
	/**
	 * Calculates indicated air temperature with given temperature and true air speed.
	 * @param t Temperature in celsius.
	 * @param tas True air speed in knots.
	 * @return Indicated air temperature in celsius.
	 */
	public double getIAT(double t, double tas)
	{
		return t + (0.95 * (tas * tas)) / 7597;
	}
	
	/**
	 * Calculates ouside actual temperature with given indicated air temperature and mach number.
	 * @param iat Indicated air temperature in celsius.
	 * @param m Mach number.
	 * @return outside actual temperature in celsius.
	 */
	public double getOAT(double iat, double m)
	{
		return (iat + 273.15) / (1 + 0.2 * 0.95 * (m * m)) - 273.15;
	}
	
	/**
	 * Calculates true air speed at low mach number (M < 0.3)
	 * @param cas Calibrated air speed in knots.
	 * @param rho Air density in slug/ft^3
	 * @param dh Altitude based on given air density.
	 * @return True air speed in knots.
	 */
	public double getTASLowMach(double cas, double rho, double dh)
	{
		if(dh < 36089.24)
			return cas / Math.pow((1 - 6.8755856E-6 * dh), 2.127940);
		else return cas * Math.pow((RHO_0SF / rho), 0.5);
	}
	
	/**
	 * Calculates true air speed.
	 * @param ias Indicated air speed in knots.
	 * @param p Pressure at given altitude in Hg.
	 * @param cs Speed of sound at given altitude.
	 * @return
	 */
	public double getTASComps(double ias, double p, double cs)
	{
		double dp = P_0H * (Math.pow((1 + 0.2 * Math.pow((ias / CS_0),2)), 3.5) - 1);
		double m = Math.pow((5 * (Math.pow((dp / p + 1), (2 / 7)) - 1)), 0.5);
		
		return m * cs;
	}
	
	/**
	 * Calculates indicated air speed with given mach numbr and altitude.
	 * @param m Mach number.
	 * @param ph Altitude with given pressure.
	 * @return Inidicated air speed in knots.
	 */
	public double getIAS(double m, double ph)
	{
		double x = Math.pow((1 - 6.8755856E-6 * ph), 5.2558797);
		return CS_0 * Math.pow(5 * Math.pow(1 + x * Math.pow(1 + Math.pow(m, 2.5), 3.5), 2 / 7.) - 1, 0.5);
	}
	
	/**
	 * Calculates drift correction with given pressure, latitude, true air speed and distance.
	 * @param p1 Pressure of first position in Hg.
	 * @param p2 Pressure of second position in Hg.
	 * @param lat Average latitude between two points in radian.
	 * @param tas True air speed in knots.
	 * @param d Distance between first and second points.
	 * @return An array of double with two elements, {Drift distance, Drift angle}, both in radians.
	 */
	public double[] getDriftCorrection(double p1, double p2, double lat, double tas, double d)
	{
		double[] data = new double[2];
		
		data[0] = 21500 * (p2 - p1) / (Math.sin(lat) * tas);
		data[1] = 1230000 * (p2 - p1) / (Math.sin(lat) * tas * d);
		
		return data;
	}
	
	/**
	 * Calculates bank radius with given velocity and bank angle.
	 * @param v Velocity in knots.
	 * @param b Bank angle in radian.
	 * @return Bank radius in feet.
	 */
	public double getBankRadius(double v, double b)
	{
		b = angleToDegree(b);
		return (v * v) / (11.23 * Math.tan(0.01745 * b));
	}
	
	/**
	 * Calculates turn rate with given bank radius and velocity.
	 * @param R Bank radius in feet.
	 * @param v Velocity in knots
	 * @return Turn rate in degrees/sec.
	 */
	public double getTurnRate(double R, double v)
	{
		return (96.7 * v) / R;
	}
	
	/**
	 * Calculates standard rate turn bank angle with give velocity.
	 * @param v Velocity in knots
	 * @return Standard rate turn bank angle in degrees.
	 */
	public double getStandardBankAngle(double v)
	{
		return 57.3 * Math.atan(v / 362.1);
	}
	
	/**
	 * Calculates pivot altitude with given velocity.
	 * @param v Velocity in knots.
	 * @return Pivotal altitude in feet.
	 */
	public double getPivotAlt(double v)
	{
		return (v * v) / 11.23;
	}
	
	/**
	 * Calculates distance to horizon with given altitude, above ground.
	 * @param h Altitude in feet.
	 * @return Horizontal distance in nautical.
	 */
	public double getDistanceToHorizon(double h)
	{
		return 1.17 * Math.sqrt(h);
	}
	
	/**
	 * Converts polar coordinate into rectangular coordinate.
	 * @param polarCoord Polar coordinate, {Degree, Distance}.
	 * @return Rectangular form, {X, Y}.
	 */
	public double[] getRectangularCoord(double[] polarCoord)
	{
		double rad = angleToRadian(polarCoord[0]);
		
		double[] coord = new double[2];
		
		coord[0] = Math.cos(rad) * polarCoord[1];
		coord[1] = Math.sin(rad) * polarCoord[1];
		
		return coord;
	}
	
	/**
	 * Converts rectangular coordinate into polar coordinate.
	 * @param rectCoord Rectangular coordinate, {X, Y}.
	 * @return Polar form, {Degree, Distance}.
	 */
	public double[] getPolarCoord(double[] rectCoord)
	{
		double[] coord = new double[2];
		
		coord[0] = angleToDegree(Math.atan(rectCoord[0] / rectCoord[1]));
		coord[1] = Math.sqrt((rectCoord[0] * rectCoord[0]) + (rectCoord[1] * rectCoord[1]));
		
		return coord;
	}
	
	/**
	 * Converts Numerical degree coordinates into string form.
	 * @param lat Latitude of coordinate.
	 * @param lon Longitude of coordinate.
	 * @return String representation of {Lat, Lon}.
	 */
	public String[] getCoordToDMS(double lat, double lon)
	{
		double modLat = lat % 1.0;
		double modLon = lon % 1.0;
		
		String degLat = String.valueOf((int)lat);
		String degLon = String.valueOf((int)lon);
		
		lat = modLat * 60.;
		modLat = lat % 1.;
		lon = modLon * 60.;
		modLon = lon % 1.;
		
		String minLat = String.valueOf((int)lat);
		String minLon = String.valueOf((int)lon);
		
		lat = modLat * 60.;
		lon = modLon * 60.;
		
		String secLat = String.valueOf((int)lat);
		String secLon = String.valueOf((int)lon);
		
		String[] coord = new String[2];
		
		coord[0] = degLat + "º" + minLat + "'" + secLat + "\"";
		coord[1] = degLon + "º" + minLon + "'" + secLon + "\"";
		
		return coord;
	}
	
	/**
	 * Modified modulus function.
	 * @param y
	 * @param x
	 * @return
	 */
	public double getMod(double y, double x)
	{
		return y - x * Math.floor(y / x);
	}
}
