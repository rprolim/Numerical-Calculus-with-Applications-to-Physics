#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>

// Define constants
const double h = 0.01;

double g(double t, double y, double z) {
	double g_func;
	
	g_func = z + y - (t * t * t) - 3 * (t * t) + 7 * t + 1;
	return g_func;
}

double y_analytical(double t) {
	double y_func;
	
	y_func = (t * t * t) - t;
	return y_func;
}

double z_analytical(double t) {
	double z_func;
	
	z_func = 3 * (t * t) - 1;
	return z_func;
}

void euler(double t, double& y, double& z, double h) {
	double y_prime, z_prime;
	
	y_prime = y;
	z_prime = z;
	y = y_prime + h * z_prime;
	z = z_prime + h * g(t, y_prime, z_prime);
}

void rk4(double t, double& y, double& z, double h) {
	double k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z;
	
	k1y = h * z;
	k1z = h * g(t, y, z);
	k2y = h * (z + k1z / 2);
	k2z = h * g(t + h / 2, y + k1y / 2, z + k1z / 2);
	k3y = h * (z + k2z / 2);
	k3z = h * g(t + h / 2, y + k2y / 2, z + k2z / 2);
	k4y = h * (z + k3z);
	k4z = h * g(t + h, y + k3y, z + k3z);
	y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
	z += (k1z + 2 * k2z + 2 * k3z + k4z) / 6;
}

int main() {
	double y = 0, z = -1; // Initial conditions
	std::cout << "#------ Euler Method ------#" << std::endl;
	for (double t = 0; t <= 5; t += h) {
		euler(t, y, z, h);
	}
	// Print in double precision (default)
	std::cout << "#------ Double Precision ------#" << std::endl;
        std::cout << " y(5): " << std::setprecision(15) << std::fixed << "Euler Method: " << y << " | " << "Analytical: " << y_analytical(5) << "\n"
                  << " dy/dt(5): " << std::setprecision(15) << std::fixed << "Euler Method: " << z << " | " << "Analytical: " << z_analytical(5) << std::endl;

        // Print in single precision (float)
        std::cout << "#------ Single Precision ------#" << std::endl;
        std::cout << " y(5): " << std::setprecision(7) << std::fixed << "Euler Method: " << static_cast<float>(y) << " | " << "Analytical: " << y_analytical(5) << "\n"
                  << " dy/dt(5): " << std::setprecision(7) << std::fixed << "Euler Method: " << static_cast<float>(z) << " | " << "Analytical: " << z_analytical(5) << std::endl;
        
        y = 0, z = -1;          
        std::cout << "#------ RK4 Method ------#" << std::endl;
	for (double t = 0; t <= 5; t += h) {
		rk4(t, y, z, h);
	}
	// Print in double precision (default)
	std::cout << "#------ Double Precision ------#" << std::endl;
        std::cout << " y(5): " << std::setprecision(15) << std::fixed << "Euler Method: " << y << " | " << "Analytical: " << y_analytical(5) << "\n"
                  << " dy/dt(5): " << std::setprecision(15) << std::fixed << "Euler Method: " << z << " | " << "Analytical: " << z_analytical(5) << std::endl;

        // Print in single precision (float)
        std::cout << "#------ Single Precision ------#" << std::endl;
        std::cout << " y(5): " << std::setprecision(7) << std::fixed << "Euler Method: " << static_cast<float>(y) << " | " << "Analytical: " << y_analytical(5) << "\n"
                  << " dy/dt(5): " << std::setprecision(7) << std::fixed << "Euler Method: " << static_cast<float>(z) << " | " << "Analytical: " << z_analytical(5) << std::endl;
	return 0;
}
