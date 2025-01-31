#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

// Define constants
const double d = 0.50; // Delta
const double w = 1.00; // Omega
const double g = 0.25; // Gamma
const double F= 0.26;

double forced(double t, double x, double v, double F) {
	double result;
	
	result = - g * v + d * x * (1 - 4 * x * x) + F * cos(w * t);
	return result;
}

void rk4_forced(double& t, double& y, double& z, double h, double F) {
	double k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z;
	
	k1y = h * z;
	k1z = h * forced(t, y, z, F);
	k2y = h * (z + k1z / 2);
	k2z = h * forced(t + h / 2, y + k1y / 2, z + k1z / 2, F);
	k3y = h * (z + k2z / 2);
	k3z = h * forced(t + h / 2, y + k2y / 2, z + k2z / 2, F);
	k4y = h * (z + k3z);
	k4z = h * forced(t + h, y + k3y, z + k3z, F);
	y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
	z += (k1z + 2 * k2z + 2 * k3z + k4z) / 6;
	t += h;
}

int main() {
	std::ofstream vXx_file("vXx.txt");
	std::cout << "Starting main function.\n" << std::endl;
	double t = 0;
	double x = -0.5, v = 0.5; // Initial conditions
	double h = 0.01 * (2 * M_PI / w); // Initial h
	for (int i = 0; i < 200000; ++i) {
		rk4_forced(t, x, v, h, F);
	}
	
	h = 0.001 * (2 * M_PI / w); // Update h
	for (int i = 0; i < 20000; ++i) {
          	// Advance one period
            	for (int j = 0; j < 1000; ++j) {
                	rk4_forced(t, x, v, h, F);
            	}
            	// Save the point on the Poincaré section
            	vXx_file << v << " " << x << "\n";
        }
	std::cout << "\nFinished!" << std::endl;
    	vXx_file.close();
    	std::cout << "Bifurcation data saved to vXx.txt\n";
    	return 0;
}
