#define _USE_MATH_DEFINES

#include <matplotlibcpp.h>
#include <glm\glm.hpp>

#include <vector>
#include <cmath>
#include <random>
#include <iostream>

namespace plt = matplotlibcpp;
using namespace glm;

const dvec3 surfaceNormal = { 0.0, 0.0, 1.0 };
const double thetaMin = 0.0;
const double thetaMax = M_PI / 2;
const double phiMin = 0.0;
const double phiMax = 2.0 * M_PI;
const double specPow = 10.0;
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 randomEngine(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<double> THETA(thetaMin, thetaMax); // Random variable Theta
std::uniform_real_distribution<double> PHI(phiMin, phiMax);   // Random variable Phi

//double f(double thetai, double phii, double thetar, double phir, double n)
//{
//	return sin(thetai) * sin(phii);
//}

// Blinn Phong.
// Expects normal vectors as arguments.
double blinnPhong(const dvec3& lightDirIn, const dvec3& lightDirOut)
{
	dvec3 halfVector = normalize(lightDirIn + lightDirOut);
	double nDotH = dot(surfaceNormal, halfVector);
	return std::pow(nDotH, specPow);
}

double lambert(const dvec3& lightDirIn, const dvec3& lightDirOut)
{
	return 1;
}

double fActual(double n) 
{
	return 8.0 * M_PI * (n + std::pow(2.0, -n / 2.0)) / ((n + 4.0) * (n + 2.0));
}

int main() 
{
	size_t inSamples = 5000;
	size_t reflSamples = 25000;

	std::vector<double> incomingAnglesTheta(inSamples), 
		incomingAnglesPhi(inSamples), 
		irradianceRatio(inSamples);
	
	double integralEst = 1;
	double maxSoFar = 0;
	double maxFromThetai = 0;
	double maxFromPhii = 0;
	double integralEstSum = 0;
	for (size_t i = 1; i <= inSamples; ++i) {
		double thetaI = THETA(randomEngine);
		double phiI = PHI(randomEngine);
		dvec3 lightDirIn;
		lightDirIn.x = sin(thetaI) * cos(phiI);
		lightDirIn.y = sin(thetaI) * sin(phiI);
		lightDirIn.z = cos(thetaI);

		double EVSum = 0;
		size_t j = 0;
		for (; j < reflSamples; ++j) {
			double thetaR = THETA(randomEngine);
			double phiR = PHI(randomEngine);
			dvec3 lightDirOut;
			lightDirOut.x = sin(thetaR) * cos(phiR);
			lightDirOut.y = sin(thetaR) * sin(phiR);
			lightDirOut.z = cos(thetaR);

			EVSum += blinnPhong(lightDirIn, lightDirOut) * cos(thetaR) * sin(thetaR);
		}

		double expectedValue = EVSum / j;
		integralEst = expectedValue * (thetaMax - thetaMin) * (phiMax - phiMin);
		integralEstSum += integralEst;
		double expectedIntegralValue = integralEstSum / i;

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
		printf("Expected integral value : %f\n", expectedIntegralValue);
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