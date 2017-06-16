/* Aviation Formulae Arduino Edition
	Initial Verion: 2017 June 16
	Last Edited:	2017 June 16*/

#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

typedef struct Constants
{
	double A = 6378.137000;
	double F = 1 / 298.257223563;
	double E_ = F * (2 - F);

	double T_0C = 15;

	double T_RCPF = 0.0019812;
	double T_TRC = -56.5;

	double P_0H = 29.92126;

	double RHO_0SF = 0.002376892;

	double CS_0 = 661.4786;
}Constant;

class AVF_Ard
{
public:
	Constant earth;
	/**
	* Modified modulus function.
	* @param y
	* @param x
	* @return
	*/
	double GetModulus(double x, double y);

	/**
	* Converts given angular value into radian.
	* @param degree Angular value in degree
	* @return Radian value of angular value.
	*/
	double AngleToRadian(double degree);

	/**
	* Converts given angular value into degree.
	* @param radian Angular value in radian
	* @return Degree value of angular value.
	*/
	double AngleToDegree(double radian);

	/**
	* Converts any nautical distance into radian distance.
	* @param nautical Nautical distance value
	* @return Radian value of nautical distance.
	*/
	double DistanceToRadian(double nautical);

	/**
	* Converts any radian distance into nautical distance.
	* @param radian Radian distance value.
	* @return Nautical value of radian distance.
	*/
	double DistanceToNautical(double radian);
	
	/**
	* Calculates great circle distance between two coordinate points.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @return Radian distance value between two points.
	*/
	double TwoPointGCD(double lat1, double lon1, double lat2, double lon2);

	/**
	* Calculates great circle distance between two coordinate points. <br>
	* This method is for shorter distance.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @return Radian distance value between two points.
	*/
	double TwoPointGDCShort(double lat1, double lon1, double lat2, double lon2);

	/**
	* Calculates true course between two points (Initial course bearing, not straight point-to-point course!)
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @return Radian angle value of initial course bearing.
	*/
	double CourseBearing(double lat1, double lon1, double lat2, double lon2);

	/**
	* Calculates true course between two points, without GCD value (Initial course bearing, not straight point-to-point course!)
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @return Radian angle value of initial course bearing.
	*/
	double CourseBearingNoGCD(double lat1, double lon1, double lat2, double lon2);

	/**
	* Calculates intermediate latitude position which lies on a great circle distance path.
	* @param lon Longitude of any given points that lies on a great circle distance
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @return Radian angle value of latitude point from given longitude.
	*/
	double GetLatPointOnGCD(double lon, double lat1, double lon1, double lat2, double lon2);

	/**
	* Calculates polar coordinate of given points from given center.<br>
	* Note that this bearing is from center, not from the given target.
	* @param lat Latitude of target point. Must be in radian.
	* @param lon Longitude of target point. Must be in radian.
	* @param lat1 Latitude of center point. Must be in radian.
	* @param lon1 Longitude of center point. Must be in radian.
	* @param bearing double pointer for bearing. This is the return value.
	* @param distance double pointer for distance. This is the return value.
	*/
	void GetCoordinatePolar(double lat, double lon, double lat1, double lon1, double* bearing, double* distance);

	/**
	* Calculates Lat/Lon point using given distance and course bearing, and center position.
	* @param d Distance from center to target. Must be in radian value.
	* @param tc Course bearing from center to target. Must be in radian value.
	* @param lat1 Latitude of center point. Must be in radian.
	* @param lon1 Longitude of center point. Must be in radian.
	* @param lat Double pointer for latitude. This is the return value.
	* @param lon Double pointer for longitutde. This is the return value.
	* @return An array of double with two elements, {Latitude, Longitude}, both in radian value.
	*/
	void GetCoordinate(double d, double tc, double lat1, double lon1, double* lat, double* lon);

	/**
	* Calculates third coordinate measured from first and second point.<br>
	* It uses given two initial course bearing to calculate a thrid point. If these two points does not intersects, returns nothing.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @param crs13 Course bearing from first point to third point. Must be in radian.
	* @param crs23 Course bearing from second point to third point. Must be in radian.
	* @param lat Double pointer for latitude. This is the return value.
	* @param lon Double pointer for longitutde. This is the return value.
	* 
	* Note: If two points does not intersects, lat and lon will return -1
	*/
	void GetIntersectRadial(double lat1, double lon1, double lat2, double lon2, double crs13, double crs23, double* lat, double* lon);

	/**
	* Calculates two passing longitude coordinate with given parallel latitude coordinate which passes through great circle distance path.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @param lat3 Latitude of thrist point which crosses GCD and is parallel. Must be in radian.
	* @param lon31 Double pointer for first longitutde. This is the return value.
	* @param lon32 Double pointer for second longitutde. This is the return value.
	*
	* Note: If given latitude coordinate does not pass through GCD path, two points return as -1.
	*/
	void GetXParallel(double lat1, double lon1, double lat2, double lon2, double lat3, double* lon31, double* lon32);

	/**
	* Calculates intermediate point that lies on a GCD, with given fractional value.
	* @param f Fractional value between 0 to 1. 0 is close to the first point, 1 is close to the second point.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @param lat Double pointer for latitude. This is the return value.
	* @param lon Double pointer for longitutde. This is the return value.
	*/
	void GetIntermedPointOnGCD(double f, double lat1, double lon1, double lat2, double lon2, double* lat, double* lon);

	/**
	* Calculates cross-track error distance.
	* @param xd Distance between starting point to arrived position. In Radian.
	* @param tx Course bearing between starting point to arrived position. In Radian.
	* @param tc Initial course bearing between starting point to actual destination position. In Radian.
	* @return Radian value of cross-track error distance.
	*/
	double GetXTD(double xd, double tx, double tc);

	/**
	* Calculates along-track distance, the distance between starting point to final point to the point right-angle to actual arrival point
	* @param xd Distance between starting point to arrived potision. In Radian.
	* @param xtd Cross-track error distance between start point to arrived position. In radian.
	* @return Radian value of along-track distance.
	*/
	double GetATD(double xd, double xtd);

	/**
	* Calculates along-track distance, the distance between starting point to final point to the point right-angle to actual arrival point
	* @param xd Distance between starting point to arrived potision. In Radian.
	* @param xtd Cross-track error distance between start point to arrived position. In radian.
	* @return Radian value of along-track distance.
	*/
	double GetATDShort(double xd, double xtd);

	/**
	* Calculates a coordinate with given initial point, target point, and an offset points.<br>
	* This uses given distance from the offset point and calculates any possible point that lies between great circle distance path and the offset point.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @param lat3 Latitude of third point. Must be in radian. (Offset Point)
	* @param lon3 Longitude of third point. Must be in radian. (Offset Point)
	* @param d Distance from third point and great circle distance path.
	* @param dir -1 or 1. (1 is forward, -1 is backward)
	* @param lat Double pointer for latitude. This is the return value.
	* @param lon Double pointer for longitutde. This is the return value.
	*
	* Note: If point does not lie within GCD, lat and lon returns -1.
	*/
	void GetCoordOffset(double lat1, double lon1, double lat2, double lon2, double lat3, double lon3, double d, double dir, double* lat, double* lon);

	/**
	* Calculates polar coordinate between given position.
	* @param lat1 Latitude of first point. Must be in radian.
	* @param lon1 Longitude of first point. Must be in radian.
	* @param lat2 Latitude of second point. Must be in radian.
	* @param lon2 Longitude of second point. Must be in radian.
	* @param bearing Double pointer for bearing. This is the return value.
	* @param distance Double pointer for bearing. This is the return value.
	*
	* Radian value is returned
	*/
	void GetApproxFlatPolarCoord(double lat1, double lon1, double lat2, double lon2, double* bearing, double* distance);

	/**
	* Calculates wind heading and speed with given course bearing, heading, true air speed and ground speed.
	* @param crs Course bearing. Must be in radian.
	* @param hdg Heading. Must be in radian.
	* @param tas True air speed. In knots.
	* @param gs Ground speed. In knots.
	* @param lat Double pointer for latitude. This is the return value. (Radian)
	* @param lon Double pointer for longitutde. This is the return value. (Radian)
	*/
	void GetWindData(double crs, double hdg, double tas, double gs, double* lat, double* lon);

	/**
	* Calculates heading and ground speed with given course bearing, true air speed, wind speed and wind direction.
	* @param crs Course bearing. Must be in radian.
	* @param tas True air spee. In knots.
	* @param ws Wind speed, in knots.
	* @param wd Wind direction. Must be in radian.
	* @param heading Double pointer for latitude. This is the return value. (Radian)
	* @param speed Double pointer for longitutde. This is the return value. (Knots)
	*/
	void GetObjectHDGGSWind(double crs, double tas, double ws, double wd, double* heading, double* speed);

	/**
	* Calculates ground speed and course bearing with given heading, true air speed, wind speed and wind direction.
	* @param hdg Heading. Must be in radian.
	* @param tas True air speed. In knots.
	* @param ws Wind speed. In knots.
	* @param wd Wind direction. Must be in radian.
	* @param gs Double pointer for ground speed in knots (Return value)
	* @param cb Double pointer for course bearing in radian (Return value)
	*/
	void GetObjectGSCRSWind(double hdg, double tas, double ws, double wd, double* gs, double* cb);

	/**
	* Calculates head-wind and cross-wind.
	* @param ws Wind speed. In knots.
	* @param wd Wind direction. Must be in radian.
	* @param hdg Heading. Must be in radian.
	* @param hw Head wind, in knots (return value)
	* @param xw Cross wind, in knots (return value)
	*/
	void GetHWXW(double ws, double wd, double hdg, double* hw, double* xw);

	/**
	* Calculates true air speed and wind speed with given three positions.
	* @param v1 Ground speed of first object. In knots.
	* @param v2 Ground speed of second object. In knots.
	* @param v3 Ground speed of third object. In knots.
	* @param tas True air speed in knots (Return value)
	* @param ws Wind speed in knots (Return value)
	*/
	void GetGPSTASWSData(double v1, double v2, double v3, double* tas, double* ws);

	/**
	* Calculates temperature with given altitude.
	* @param h Altitude in feet.
	* @return Temperature in celsius.
	*/
	double GetTemp(double h);

	/**
	* Calculates pressure at given altitude.
	* @param h Altitude in feet.
	* @return Pressure in Hg.
	*/
	double GetPressure(double h);

	/**
	* Calculates air density at given altitude.
	* @param h Altitude in feet.
	* @return Air density in slugs/feet^3
	*/
	double GetRho(double h);

	/**
	* Calculates corrected altitude with given pressure and altitude.
	* @param p Pressure in Hg.
	* @param h Altitude in feet.
	* @return Corrected altitude in feet.
	*/
	double GetCorrAlt(double p, double h);

	/**
	* Calculates altitude based on air density.
	* @param ph Altitude corrected by pressure. In feet.
	* @param th Temperature at given altitude. In celsius.
	* @param t Actual temperature at given altitude. In celcius.
	* @return Corrected altitude in feet.
	*/
	double GetRhoAlt(double ph, double th, double t);

	/**
	* Calculates mach number with given true air speed and speed of sound.
	* @param tas True air speed, in knots.
	* @param cs Speed of sound, in knots.
	* @return Mach number.
	*/
	double GetMachNumber(double tas, double cs);

	/**
	* Calculates speed of sound at given temperature
	* @param t Temperature in celsius.
	* @return Speed of sound in knots.
	*/
	double GetCS(double t);

	/**
	* calculates indicated air temperature with given temperature and true air speed.
	* @param t temperature in celsius.
	* @param tas true air speed in knots.
	* @return indicated air temperature in celsius.
	*/
	double Getiat(double t, double tas);

	/**
	* Calculates ouside actual temperature with given indicated air temperature and mach number.
	* @param iat Indicated air temperature in celsius.
	* @param m Mach number.
	* @return outside actual temperature in celsius.
	*/
	double GetOAT(double iat, double m);

	/**
	* Calculates true air speed at low mach number (M < 0.3)
	* @param cas Calibrated air speed in knots.
	* @param rho Air density in slug/ft^3
	* @param dh Altitude based on given air density.
	* @return True air speed in knots.
	*/
	double GetTASLowMach(double cas, double rho, double dh);

	/**
	* Calculates true air speed.
	* @param ias Indicated air speed in knots.
	* @param p Pressure at given altitude in Hg.
	* @param cs Speed of sound at given altitude.
	* @return True air speed
	*/
	double GetTASComps(double ias, double p, double cs);

	/**
	* Calculates indicated air speed with given mach numbr and altitude.
	* @param m Mach number.
	* @param ph Altitude with given pressure.
	* @return Inidicated air speed in knots.
	*/
	double GetIAS(double m, double ph);

	/**
	* Calculates drift correction with given pressure, latitude, true air speed and distance.
	* @param p1 Pressure of first position in Hg.
	* @param p2 Pressure of second position in Hg.
	* @param lat Average latitude between two points in radian.
	* @param tas True air speed in knots.
	* @param d Distance between first and second points.
	* @param dd Double pointer for drift distance in radians (Return value)
	* @param da Double pointer for drift angle in radians (Return value)
	*/
	void GetDriftCorrection(double p1, double p2, double lat, double tas, double d, double* dd, double* da);

	/**
	* Calculates bank radius with given velocity and bank angle.
	* @param v Velocity in knots.
	* @param b Bank angle in radian.
	* @return Bank radius in feet.
	*/
	double GetBankRadius(double v, double b);

	/**
	* Calculates turn rate with given bank radius and velocity.
	* @param R Bank radius in feet.
	* @param v Velocity in knots
	* @return Turn rate in degrees/sec.
	*/
	double GetTurnRate(double R, double v);

	/**
	* Calculates standard rate turn bank angle with give velocity.
	* @param v Velocity in knots
	* @return Standard rate turn bank angle in degrees.
	*/
	double GetStandardBankAngle(double v);

	/**
	* Calculates pivot altitude with given velocity.
	* Pivotal altitude is the altitude with fixed bank angle, and pivot point
	* is stationary.
	* Note that this formula works with speed less than 250 knots
	* @param v Velocity in knots.
	* @return Pivotal altitude in feet.
	*/
	double GetPivotAlt(double v);

	/**
	* Calculates distance to horizon with given altitude, above ground.
	* @param h Altitude in feet.
	* @return Horizontal distance in nautical.
	*/
	double GetDistanceToHorizon(double h);

	/**
	* Converts polar coordinate into rectangular coordinate.
	* @param polarCoord Polar coordinate, {Degree, Distance}.
	* @param degree Double pointer for degree
	* @param dist Double pointer for distance
	* @param x Double poitner for X coord (Return Value)
	* @param y Double pointer for Y coord (Return Value)
	*/
	void GetRectangularCoord(double* degree, double* dist, double* x, double* y);

	/**
	* Converts rectangular coordinate into polar coordinate.
	* @param x Double pointer for X coord
	* @param y Double pointer for Y coord
	* @param degree Double pointer for degree (Return Value)
	* @param dist Double pointer for distance (Return Value)
	*/
	void GetPolarCoord(double* x, double* y, double* degree, double* dist);

	/**
	* Converts Numerical degree coordinates into string form.
	* @param lat Latitude of coordinate.
	* @param lon Longitude of coordinate.
	* @return String representation of {Lat, Lon}.
	*/
	void GetCoordToDMS(double lat, double lon, char* slat, char* slon);
};