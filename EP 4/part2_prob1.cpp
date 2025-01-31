#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

// Define constants
const double w = 1.00; // Omega

double func(double t, double x, double v, double F, double g, int item) {
	double result;
	if (item == 1) {
		result = 0.5 * x * (1 - 4 * x * x);
	}
	if (item == 2) {
		result = 0.5 * x * (1 - 4 * x * x) - g * v;
	}
	if (item == 3) {
		result = F * cos(w * t) + 0.5 * x * (1 - 4 * x * x) - 0.25 * v;
	}
	return result;
}

void rk4(double& t, double& y, double& z, double h, double F, double g, int item) {
	double k1y, k1z, k2y, k2z, k3y, k3z, k4y, k4z;
	
	k1y = h * z;
	k1z = h * func(t, y, z, F, g, item);
	k2y = h * (z + k1z / 2);
	k2z = h * func(t + h / 2, y + k1y / 2, z + k1z / 2, F, g, item);
	k3y = h * (z + k2z / 2);
	k3z = h * func(t + h / 2, y + k2y / 2, z + k2z / 2, F, g, item);
	k4y = h * (z + k3z);
	k4z = h * func(t + h, y + k3y, z + k3z, F, g, item);
	y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
	z += (k1z + 2 * k2z + 2 * k3z + k4z) / 6;
	t += h;
}

int main() {
	// Item a)
	std::vector<std::ofstream> files_a(3);
	double v_vec[3] = {0.1, 0.25, 0.5};
	for (int i = 0; i < 3; ++i) {
		std::ostringstream filename;
        	filename << "item_a-" << i + 1 << ".txt";
        	files_a[i].open(filename.str());
		
		double t = 0.0, x = -0.5; // Initial conditions
		double v = v_vec[i];
		double h = 0.001;
		double F = 0, g = 0;
		for (int j = 0; j < 1000000; ++j) {
			rk4(t, x, v, h, F, g, 1);
			files_a[i] << x << " " << v << "\n";
		}	
		files_a[i].close();
	}
    	
    	// Item b)
    	std::vector<std::ofstream> files_b(2);
	double g_vec[3] = {0.25, 0.8};
	for (int i = 0; i < 2; ++i) {
		std::ostringstream filename;
        	filename << "item_b-" << i + 1 << ".txt";
        	files_b[i].open(filename.str());
        	
		double t = 0.0, x = -0.5, v = 0.5; // Initial conditions
		double g = g_vec[i];
		double h = 0.001;
		double F = 0;
		for (int j = 0; j < 100000; ++j) {
			rk4(t, x, v, h, F, g, 2);
			files_b[i] << x << " " << v << "\n";
		}
		files_b[i].close();	
	}
	
	// Item c)
	std::vector<std::ofstream> files_c(4);
	double F_vec[5] = {0.11, 0.115, 0.14, 0.35};
	for (int i = 0; i < 4; ++i) {
		std::ostringstream filename;
        	filename << "item_c-" << i + 1 << ".txt";
        	files_c[i].open(filename.str());
        	
		double t = 0.0, x = -0.5, v = 0.5; // Initial conditions
		double F = F_vec[i];
		double h = 0.01; // Initial h
		double g = 0.0;
		for (int j = 0; j < 200000; ++j) {
			rk4(t, x, v, h, F, g, 3);
		}
		h = 0.001; // Update h
		for (int j = 0; j < 200000; ++j) {
			rk4(t, x, v, h, F, g, 3);
			files_c[i] << x << " " << v << "\n";
		}
		files_c[i].close();	
	}
    	return 0;
}
