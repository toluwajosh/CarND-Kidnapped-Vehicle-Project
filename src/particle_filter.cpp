/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

//using namespace std;

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	
	// Number of particles
	num_particles = 1000;

	particles.resize(num_particles);

	// Initialize weights to 1
	weights = std::vector<double>(num_particles);
	
	std::fill(weights.begin(), weights.end(), 1);

	// Standard deviations for x, y, and psi
	double std_x, std_y, std_theta; 
	std_x = 2;
	std_y = 2;
	std_theta = 0.05;
	//
	std[0] = std_x;
	std[1] = std_y;
	std[2] = std_theta;

	// Create a normal (Gaussian) distribution for x, y and theta.
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	// Initialize all particles
	for (int i = 0; i < num_particles; i++)
	{
		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine gen;

	

	for ( int i = 0; i < num_particles; i++)
	{
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;

		std::normal_distribution<double> dist_x(x, std_pos[0]);
		std::normal_distribution<double> dist_y(y, std_pos[1]);
		std::normal_distribution<double> dist_theta(theta, std_pos[2]);

		particles[i].x = dist_x(gen) + (velocity / yaw_rate)*(sin(theta + yaw_rate * delta_t) - sin(theta));
		particles[i].y = dist_y(gen) + (velocity / yaw_rate)*(cos(theta) - cos(theta + yaw_rate + delta_t));
		particles[i].theta = dist_theta(gen) + yaw_rate * delta_t;
	}


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	for (int i = 0; i < predicted.size(); i++)
	{
		for (int j = 0; j < observations.size(); j++)
		{
			double distance = dist(predicted[i].x, predicted[i].y, observations[i].x, observations[i].y);
			int choice_id = j;
		}
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
