#define _USE_MATH_DEFINES

#include <matplotlibcpp.h>
#include <glm\glm.hpp>

#include <vector>
#include <cmath>
#include <random>
#include <iostream>

namespace plt = matplotlibcpp;
using namespace glm;

//double f(double thetai, double phii, double thetar, double phir, double n)
//{
//	return sin(thetai) * sin(phii);
//}

double blinnPhongModified(double thetai, double phii, double thetar, double phir, double n) 
{
	return pow((cos(thetai) + cos(thetar))
	          / sqrt(2 * (1 + sin(thetai) * cos(phii) * sin(thetar) * cos(phir)
	                        + cos(thetai) * cos(thetar)
	                        + sin(thetai) * sin(phii) * sin(thetar) * sin(phir))),
	       n);
}

double fActual(double n) 
{
	return 8 * M_PI * (n + pow(2, -n / 2)) / ((n + 4) * (n + 2));
}

int main() 
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 randomEngine(rd()); //Standard mersenne_twister_engine seeded with rd()

	double a = 0;
	double b = M_PI / 2;
	double c = 0;
	double d = 2 * M_PI;
	double n = 10;
	std::uniform_real_distribution<double> THETA(a, b);
	std::uniform_real_distribution<double> PHI(c, d);

	std::cout << "\nIntegrating f(x) from 0 to pi...\n\n";

	size_t inSamples = 1000;
	size_t reflSamples = 50000;

	std::vector<double> incomingAnglesTheta(inSamples), incomingAnglesPhi(inSamples), irradianceRatio(inSamples);
	
	double expectedValue = 1;
	double integralEst = 1;
	double maxSoFar = 0;
	double maxFromThetai = 0;
	double maxFromPhii = 0;
	for (size_t i = 1; i <= inSamples; ++i) {

		double thetai = THETA(randomEngine);
		double phii = PHI(randomEngine);
		
		double sum = 0;
		for (size_t j = 1; j <= reflSamples; ++j) {
			double thetar = THETA(randomEngine);
			double phir = PHI(randomEngine);
			sum += blinnPhongModified(thetai, phii, thetar, phir, n) * cos(thetar) * sin(thetar);

			if (j == reflSamples) {
				expectedValue = sum / j;
				integralEst = expectedValue * (b - a) * (d - c);

				incomingAnglesTheta.at(i - 1) = thetai;
				incomingAnglesPhi.at(i - 1) = phii;
				irradianceRatio.at(i - 1) = integralEst;

				if (integralEst > maxSoFar) {
					maxSoFar = integralEst;
					maxFromThetai = thetai;
					maxFromPhii = phii;
				}

				printf("\nMonte Carlo E[f(x)] estimate after %zu samples is : %f\n", j, expectedValue);
				printf("Theta incoming is : %f\n", thetai);
				printf("phi incoming is : %f\n", phii);
				printf("Monte Carlo integral(f(x)dx) estimate after %zu samples is : %f\n", j, integralEst);
				printf("Pressumed Max : %f\n", fActual(n));
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