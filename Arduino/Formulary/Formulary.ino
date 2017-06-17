#include <AVF_Ard.h>

/* Consider this as sample sketch */
static AVF_Ard avf;

void setup()
{
  Serial.begin(9600);
  Serial.println("Arduino Aviation Formulary");
  
  /* AngleToRadian returns degree angle to radian
      Here, we will convert 60 to degree
      Returns 1.05*/
  Serial.print("AngleToRadian(60.0): ");
  Serial.println(avf.AngleToRadian(60.0));

  /* AngleToDegree returns radian angle to degree
      Here, we will convert 1.05 to degree
      Returns 60.16*/
  Serial.print("AngleToDegree(1.05): ");
  Serial.println(avf.AngleToDegree(1.05));

  /* DistanceToRadian returns radian distance from
      nautical mile. Here, we will convert 2144 nmi
      to radians returns 0.62 radians*/
  Serial.print("DistanceToRadian(2144): ");
  Serial.println(avf.DistanceToRadian(2144));

  /* DistanceToNautical returns nautical mile distance
      from radian distance. Here, we will convert 0.62
      radians to nautical miles returns 2131.40*/
  Serial.print("DistanceToNautical(0.62): ");
  Serial.println(avf.DistanceToNautical(0.62));
}

void loop()
{
}

