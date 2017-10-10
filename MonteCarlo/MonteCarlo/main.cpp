#define _USE_MATH_DEFINES

#include <matplotlibcpp.h>
#include <glm\glm.hpp>
#include <glm\gtc\random.hpp>

#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>

namespace plt = matplotlibcpp;
using namespace glm;

const dvec3 surfaceNormal = { 0, 0, 1 };

//double f(double thetai, double phii, double thetar, double phir, double n)
//{
//	return sin(thetai) * sin(phii);
//}

// Blinn Phong.
// Expects normal vectors as arguments.
double blinnPhong(dvec3 lightDirIn, dvec3 lightDirOut, double specPow)
{
	dvec3 halfVector = normalize(lightDirIn + lightDirOut);
	double nDotH = dot(surfaceNormal, halfVector);
	return pow(nDotH, specPow);
}

double fActual(double n) 
{
	return 8 * M_PI * (n + pow(2, -n / 2)) / ((n + 4) * (n + 2));
}

int main() 
{
	srand(time(0));

	double a = 0;
	double b = M_PI / 2;
	double c = 0;
	double d = 2 * M_PI;
	double specPow = 10;

	size_t inSamples = 3000;
	size_t reflSamples = 10000;

	std::vector<double> incomingAnglesTheta(inSamples), 
		incomingAnglesPhi(inSamples), 
		irradianceRatio(inSamples);
	
	double expectedValue = 1;
	double integralEst = 1;
	double maxSoFar = 0;
	double maxFromThetai = 0;
	double maxFromPhii = 0;
	for (size_t i = 1; i <= inSamples; ++i) {

		dvec3 lightDirIn = glm::sphericalRand(1.0); // Magnitude 1 vector on sphere
		lightDirIn.z = std::abs(lightDirIn.z); // Reflect to get only top hemisphere
		lightDirIn = normalize(lightDirIn);
		lightDirIn = { 0, 0, 1 };

		double sum = 0;
		for (size_t j = 1; j <= reflSamples; ++j) {
			dvec3 lightDirOut = glm::sphericalRand(1.0); // Magnitude 1 vector on sphere
			lightDirOut.z = std::abs(lightDirOut.z); // Reflect to get only top hemisphere
			lightDirOut = normalize(lightDirOut);

			float cosThetaR = dot(surfaceNormal, lightDirOut);
			float sinThetaR = length(cross(surfaceNormal, lightDirOut));
			sum += blinnPhong(lightDirIn, lightDirOut, specPow) * cosThetaR * sinThetaR;

			if (j == reflSamples) {
				expectedValue = sum / j;
				integralEst = expectedValue * (b - a) * (d - c);

				double thetaI = std::acos(lightDirIn.z);
				double phiI = std::atan2(lightDirIn.y, lightDirIn.x);
				incomingAnglesTheta.at(i - 1) = thetaI;
				incomingAnglesPhi.at(i - 1) = phiI;
				irradianceRatio.at(i - 1) = integralEst;

				if (integralEst > maxSoFar) {
					maxSoFar = integralEst;
					maxFromThetai = thetaI;
					maxFromPhii = phiI;
				}

				printf("\nMonte Carlo E[f(x)] estimate after %zu samples is : %f\n", j, expectedValue);
				printf("Theta incoming is : %f\n", thetaI);
				printf("phi incoming is : %f\n", phiI);
				printf("Monte Carlo integral(f(x)dx) estimate after %zu samples is : %f\n", j, integralEst);
				printf("Pressumed Max : %f\n", fActual(specPow));
				printf("Max so far : %f, from Zenith angle : %f, Azimuth angle : %f\n", maxSoFar, maxFromThetai, maxFromPhii);
			}
		}
	}

	plt::plot(incomingAnglesTheta, irradianceRatio, "ro");
	plt::ylabel("Ratio of outgoing to incoming irradiance");
	plt::xlabel("Incoming light angle (Zenith)");
	plt::show();

	plt::plot(incomingAnglesPhi, irradianceRatio, "ro");
	plt::ylabel("Ratio of outgoing to incoming irradiance");
	plt::xlabel("Incoming light angle (Azimuth)");
	plt::show();
}