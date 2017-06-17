/* Aviation Formulae Arduino Edition
Initial Verion: 2017 June 16
Last Edited:	2017 June 16*/

#if ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

#include <AVF_Ard.h>
#include <Math.h>
#include <string.h>

double AVF_Ard::GetModulus(double x, double y)
{
	return y - x * floor(y / x);
}

double AVF_Ard::AngleToRadian(double degree)
{
	return (M_PI / 180.0) * degree;
}

double AVF_Ard::AngleToDegree(double radian)
{
	return (180.0 / M_PI) * radian;
}

double AVF_Ard::DistanceToRadian(double nautical)
{
	return (M_PI / (60.0 * 180.0)) * nautical;
}

double AVF_Ard::DistanceToNautical(double radian)
{
	return ((180.0 * 60.0) / M_PI) * radian;
}

double AVF_Ard::TwoPointGCD(double lat1, double lon1, double lat2, double lon2)
{
	return acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2));
}

double AVF_Ard::TwoPointGDCShort(double lat1, double lon1, double lat2, double lon2)
{
	return 2 * asin(sqrt((sin((lat1 - lat2) / 2) * sin((lat1 - lat2) / 2)) + cos(lat1) * cos(lat2) * (sin((lon1 - lon2) / 2) * sin((lon1 - lon2) / 2))));
}

double AVF_Ard::CourseBearing(double lat1, double lon1, double lat2, double lon2)
{
	double d = TwoPointGCD(lat1, lon1, lat2, lon2);

	if (sin(lon2 - lon1) < 0)
		return acos((sin(lat2) - sin(lat1) * cos(d)) / (sin(d) * cos(lat1)));
	else
		return 2 * M_PI - acos((sin(lat2) - sin(lat1) * cos(d)) / (sin(d) * cos(lat1)));
}

double AVF_Ard::CourseBearingNoGCD(double lat1, double lon1, double lat2, double lon2)
{
	return GetModulus(atan2(sin(lon1 - lon2) * cos(lat2), cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(lon1 - lon2)), 2 * M_PI);
}

double AVF_Ard::GetLatPointOnGCD(double lon, double lat1, double lon1, double lat2, double lon2)
{
	if (sin(lon1 - lon2) == 0)
		return 0;

	return atan((sin(lat1) * cos(lat2) * sin(lon - lon2) - sin(lat2) * cos(lat1) * sin(lon - lon1)) / (cos(lat1) * cos(lat2) * sin(lon1 - lon2)));
}

void AVF_Ard::GetCoordinatePolar(double lat, double lon, double lat1, double lon1, double* bearing, double* distance)
{
	*bearing = CourseBearing(lat1, lon1, lat, lon);
	*distance = TwoPointGCD(lat1, lon1, lat, lon);
}

void AVF_Ard::GetCoordinate(double d, double tc, double lat1, double lon1, double* lat, double* lon)
{
	*lat = asin(sin(lat1) * cos(d) + cos(lat1) * sin(d) * cos(tc));
	double dlon = atan2(sin(tc) * sin(d) * cos(lat1), cos(d) - sin(lat1) * sin(*lat));
	*lon = GetModulus(lon1 - dlon + M_PI, 2 * M_PI) - M_PI;
}

void AVF_Ard::GetIntersectRadial(double lat1, double lon1, double lat2, double lon2, double crs13, double crs23, double* lat, double* lon)
{
	double dst12 = TwoPointGCD(lat1, lon1, lat2, lon2);
	double crs12 = CourseBearing(lat1, lon1, lat2, lon2);
	double crs21 = CourseBearing(lat2, lon2, lat1, lon1);

	double ang1 = GetModulus(crs13 - crs12 + M_PI, 2. * M_PI) - M_PI;
	double ang2 = GetModulus(crs21 - crs23 + M_PI, 2. * M_PI) - M_PI;

	if (sin(ang1) == 0 && sin(ang2) == 0)
	{
		*lat = -1; *lon = -1;
	}
	else if (sin(ang1) * sin(ang2) < 0)
	{
		*lat = -1; *lon = -1;
	}
	else
	{
		ang1 = abs(ang1);
		ang2 = abs(ang2);
		double ang3 = acos(-cos(ang1) * cos(ang2) + sin(ang1) * sin(ang2) * cos(dst12));
		double dst13 = atan2(sin(dst12) * sin(ang1) * sin(ang2), cos(ang2) + cos(ang1) * cos(ang3));
		GetCoordinate(dst13, ang3, lat1, lon1, lat, lon);
	}
}

void AVF_Ard::GetXParallel(double lat1, double lon1, double lat2, double lon2, double lat3, double* lon31, double* lon32)
{
	double lon12 = lon1 - lon2;
	double A = sin(lat1) * cos(lat2) * cos(lat3) * sin(lon12);
	double B = sin(lat1) * cos(lat2) * cos(lat3) * cos(lon12) - cos(lat1) * sin(lat2) * cos(lat3);
	double C = cos(lat1) * cos(lat2) * sin(lat3) * sin(lon12);

	double lon = atan2(B, A);

	if (abs(C) > sqrt((A * A) + (B * B)))
	{
		*lon31 = -1; *lon32 = -1;
	}
	else
	{
		double dlon = acos(C / sqrt((A * A) + (B * B)));
		*lon31 = GetModulus(lon1 + dlon + lon + M_PI, 2 * M_PI) - M_PI;
		*lon32 = GetModulus(lon1 - dlon + lon + M_PI, 2 * M_PI) - M_PI;
	}
}

void AVF_Ard::GetIntermedPointOnGCD(double f, double lat1, double lon1, double lat2, double lon2, double* lat, double* lon)
{
	double d = TwoPointGCD(lat1, lon1, lat2, lon2);

	double A = sin((1 - f) * d) / sin(d);
	double B = sin(f * d) / sin(d);

	double x = A * cos(lat1) * cos(lon1) + B * cos(lat2) * cos(lon2);
	double y = A * cos(lat1) * sin(lon1) + B * cos(lat2) * sin(lon2);
	double z = A * sin(lat1) + B * sin(lat2);

	*lat = atan2(z, sqrt((x * x) + (y * y)));
	*lon = atan2(y, x);
}

double AVF_Ard::GetXTD(double xd, double tx, double tc)
{
	return asin(sin(xd) * sin(tx - tc));
}

double AVF_Ard::GetATD(double xd, double xtd)
{
	return acos(cos(xd) / cos(xtd));
}

double AVF_Ard::GetATDShort(double xd, double xtd)
{
	return asin(sqrt((sin(xd) * sin(xd)) - (sin(xtd) * sin(xtd))) / cos(xtd));
}

void AVF_Ard::GetCoordOffset(double lat1, double lon1, double lat2, double lon2, double lat3, double lon3, double d, double dir, double* lat, double* lon)
{
	double crsA = CourseBearing(lat1, lon1, lat3, lon3) - CourseBearing(lat1, lon1, lat2, lon2);
	double dst13 = TwoPointGCD(lat1, lon1, lat3, lon3);
	double XTD = pow((cos(dst13) * cos(dst13)) + (sin(dst13) * sin(dst13)) * (cos(crsA) * cos(crsA)), 1.0 / 2.0);
	double ATD = atan2(sin(dst13) * cos(crsA), cos(dst13));

	if ((cos(d) * cos(d)) > (XTD * XTD))
	{
		*lat = -1; *lon = -1;
	}
	else
	{
		double dst = (ATD * dir) + acos(cos(d) / XTD);
		double tc = CourseBearing(lat1, lon1, lat3, lon3);
		GetCoordinate(dst, tc, lat1, lon1, lat, lon);
	}
}

void AVF_Ard::GetApproxFlatPolarCoord(double lat1, double lon1, double lat2, double lon2, double* bearing, double* distance)
{
	double R1 = earth.A * (1 - earth.E_) / pow(1 - earth.E_ * (sin(lat1) * sin(lat1)), 3.0 / 2.0);
	double R2 = earth.A / sqrt(1 - earth.E_ * (sin(lat1) * sin(lat1)));

	double dstN = R1 * (lat2 - lat1);
	double dstE = R2 * cos(lat1) * (lon2 - lon1);

	*bearing = GetModulus(atan2(dstE, dstN), 2 * M_PI);
	*distance = sqrt((dstN * dstN) + (dstE * dstE));
}

void AVF_Ard::GetWindData(double crs, double hdg, double tas, double gs, double* lat, double* lon)
{
	*lat = crs + atan2(tas * sin(hdg - crs), tas * cos(hdg - crs) - gs);
	*lat = GetModulus(*lat, 2 * M_PI);
	*lon = sqrt(pow(tas - gs, 2) + 4 * tas * gs * pow(sin((hdg - crs) / 2), 2));
}

void AVF_Ard::GetObjectHDGGSWind(double crs, double tas, double ws, double wd, double* heading, double* speed)
{
	double swc = (ws / tas) * sin(wd - crs);

	if (abs(swc) > 1)
	{
		*heading = -1; *speed = -1;
	}
	else
	{
		*heading = crs + asin(swc);
		*heading = GetModulus(*heading, 2 * M_PI);

		*speed = tas * sqrt(1 - (swc * swc)) - ws * cos(wd - crs);
		if (*speed < 0)
		{
			*heading = -1; *speed = -1;
		}
	}
}

void AVF_Ard::GetObjectGSCRSWind(double hdg, double tas, double ws, double wd, double* gs, double* cb)
{
	*gs = sqrt((ws * ws) + (tas * tas) - 2.0 * ws * tas * cos(hdg - wd));
	double wca = atan2(ws * sin(hdg - wd), tas - wd * cos(hdg - wd));
	if (ws < tas)
		wca = asin((ws / *gs) * sin(hdg - wd));
	*cb = GetModulus(hdg + wca, 2.0 * M_PI);
}

void AVF_Ard::GetHWXW(double ws, double wd, double hdg, double* hw, double* xw)
{
	*hw = ws * cos(wd - hdg);
	*xw = ws * sin(wd - hdg);
}

void AVF_Ard::GetGPSTASWSData(double v1, double v2, double v3, double* tas, double* ws)
{
	double vms = (pow(v1, 2) + pow(v2, 2.0) + pow(v3, 2.0)) / 3.0;
	double a1 = (v1 * v1) / (vms - 1);
	double a2 = (v2 * v2) / (vms - 1);
	double a3 = (v3 * v3) / (vms - 1);
	double mu = (pow(a1, 2.0) + pow(a2, 2.0) + pow(a3, 2.0)) / 6.0;

	double bp = (1.0 / 2.0) + sqrt(1.0 / (4.0 - mu));
	double bm = mu / bp;

	*tas = sqrt(vms * bp);
	*ws = sqrt(vms * bm);

	if (*tas < *ws)
	{
		*tas = sqrt(vms * bm);
		*ws = sqrt(vms * bp);
	}
}

double AVF_Ard::GetTemp(double h)
{
	if (h < 36089.24)
		return earth.T_0C - (earth.T_RCPF * h);
	else
		return earth.T_TRC;
}

double AVF_Ard::GetPressure(double h)
{
	if (h < 36089.24)
		return earth.P_0H * pow(1 - 6.8755856E-6 * h, 5.2558797);
	else
	{
		double p_tr = 0.2233609 * earth.P_0H;
		return p_tr * exp(-4.806346E-5 * (h - 36089.24));
	}
}

double AVF_Ard::GetRho(double h)
{
	if (h < 36089.24)
		return earth.RHO_0SF * pow((1. - 6.8755856E-6 * h), 4.2558797);
	else
	{
		double rho_Tr = 0.2970756 * earth.RHO_0SF;
		return rho_Tr * exp(-4.806346E-5 * (h - 36089.24));
	}
}

double AVF_Ard::GetCorrAlt(double p, double h)
{
	double corr = 145442.2 * (1.0 - pow((p / 29.92126), 0.190261));
	return h + corr;
}

double AVF_Ard::GetRhoAlt(double ph, double th, double t)
{
	th += 273.15;
	t += 273.15;

	return ph + (th / earth.T_RCPF) * (1.0 - pow((th / t), 0.2349690));
}

double AVF_Ard::GetMachNumber(double tas, double cs)
{
	return tas / cs;
}

double AVF_Ard::GetCS(double t)
{
	return 38.967854 * sqrt(t + 273.15);
}

double AVF_Ard::Getiat(double t, double tas)
{
	return t + (0.95 * (tas * tas)) / 7597.0;
}

double AVF_Ard::GetOAT(double iat, double m)
{
	return (iat + 273.15) / (1.0 + 0.2 * 0.95 * (m * m)) - 273.15;
}

double AVF_Ard::GetTASLowMach(double cas, double rho, double dh)
{
	if (dh < 36089.24)
		return cas / pow((1 - 6.8755856E-6 * dh), 2.127940);
	else return cas * pow((earth.RHO_0SF / rho), 0.5);
}

double AVF_Ard::GetTASComps(double ias, double p, double cs)
{
	double dp = earth.P_0H * (pow((1.0 + 0.2 * pow((ias / earth.CS_0), 2.0)), 3.5) - 1.0);
	double m = pow((5.0 * (pow((dp / p + 1.0), (2.0 / 7.0)) - 1.0)), 0.5);

	return m * cs;
}

double AVF_Ard::GetIAS(double m, double ph)
{
	double x = pow((1 - 6.8755856E-6 * ph), 5.2558797);
	return earth.CS_0 * pow(5.0 * pow(1.0 + x * pow(1.0 + pow(m, 2.5), 3.5), 2.0 / 7.0) - 1.0, 0.5);
}

void AVF_Ard::GetDriftCorrection(double p1, double p2, double lat, double tas, double d, double* dd, double* da)
{
	*dd = 21500.0 * (p2 - p1) / (sin(lat) * tas);
	*da = 1230000.0* (p2 - p1) / (sin(lat) * tas * d);
}

double AVF_Ard::GetBankRadius(double v, double b)
{
	b = AngleToDegree(b);
	return (v * v) / (11.23 * tan(0.01745 * b));
}

double AVF_Ard::GetTurnRate(double R, double v)
{
	return (96.7 * v) / R;
}

double AVF_Ard::GetStandardBankAngle(double v)
{
	return 57.3 * atan(v / 362.1);
}

double AVF_Ard::GetM_PIvotAlt(double v)
{
	return (v * v) / 11.23;
}

double AVF_Ard::GetDistanceToHorizon(double h)
{
	return 1.17 * sqrt(h);
}

void AVF_Ard::GetRectangularCoord(double* degree, double* dist, double* x, double* y)
{
	double rad = AngleToRadian(*degree);

	*x = cos(rad) * (*dist);
	*y = sin(rad) * (*dist);
}

void AVF_Ard::GetPolarCoord(double* x, double* y, double* degree, double* dist)
{
	*degree = AngleToDegree(atan((*x) / (*y)));
	*dist = sqrt(((*x) * (*x)) + ((*y) * (*y)));
}

void AVF_Ard::GetCoordToDMS(double lat, double lon, char* slat, char* slon)
{
	double modLat = GetModulus(lat, 1.0);
	double modLon = GetModulus(lon, 1.0);

	const char degLat[3];
	const char degLon[3];
	itoa((int)lat, degLat, 10);
	itoa((int)lon, degLon, 10);

	lat = modLat * 60.0;
	modLat = GetModulus(lat, 1.0);
	lon = modLon * 60.0;
	modLon = GetModulus(lon, 1.0);

	const char minLat[3];
	const char minLon[3];
	itoa((int)lat, minLat, 10);
	itoa((int)lon, minLon, 10);

	lat = modLat * 60.0;
	lon = modLon * 60.0;

	const char secLat[3];
	const char secLon[3];
	itoa((int)lat, secLat, 10);
	itoa((int)lon, secLon, 10);

	//String[] coord = new String[2];

	//slat = degLat + "�" + minLat + "'" + secLat + "\"";
	strcat(slat, degLat); strcat(slat, "�");
	strcat(slat, minLat); strcat(slat, "\'");
	strcat(slat, secLat); strcat(slat, "\"");
	//slon = degLon + "�" + minLon + "'" + secLon + "\"";
	strcat(slon, degLat); strcat(slon, "�");
	strcat(slon, minLat); strcat(slon, "\'");
	strcat(slon, secLat); strcat(slon, "\"");
}
