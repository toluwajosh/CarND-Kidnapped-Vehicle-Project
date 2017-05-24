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

#include "particle_filter.h"

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Number of particles
	num_particles = 100;
	particles.resize(num_particles);

	// Create a normal (Gaussian) distribution for x, y and theta.
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(0, std[0]);
	std::normal_distribution<double> dist_y(0, std[1]);
	std::normal_distribution<double> dist_theta(0, std[2]);

	// Initialize all particles
	for (int i = 0; i < num_particles; i++)
	{
		particles[i].id = i;
		particles[i].x = x + dist_x(gen);
		particles[i].y = y + dist_y(gen);
		particles[i].theta = theta + dist_theta(gen);
		// set weights to 1
		particles[i].weight = 1.0;
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// Added noise using std::normal_distribution and std::default_random_engine.
	// Ref:
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/

	// declare random engine
	std::default_random_engine gen;

	for (int i = 0; i < num_particles; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		
		double v_yaw_rate = velocity / yaw_rate;
		double yaw_delta_t = yaw_rate * delta_t;
		double v_dt = velocity * delta_t;

		// prediction
		if (fabs(yaw_rate) > 1e-5){
			x += v_yaw_rate*(sin(theta + yaw_delta_t) - sin(theta));
			y += v_yaw_rate*(cos(theta) - cos(theta + yaw_delta_t));
			theta += yaw_delta_t;
		}
		else{
			x += v_dt*sin(theta);
			y += v_dt*cos(theta);
			// no change in theta since car is moving straight
		}

		// Add noise
		std::normal_distribution<double> dist_x(0, std_pos[0]);
		std::normal_distribution<double> dist_y(0, std_pos[1]);
		std::normal_distribution<double> dist_theta(0, std_pos[2]);

		particles[i].id = i;
		particles[i].x = x + dist_x(gen);
		particles[i].y = y + dist_y(gen);
		particles[i].theta = theta + dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and,
	// assign the observed measurement to this particular landmark.

	// Note:
	// predicted(here) is observation in filter
	// observation(here) is landmarks within range in filter
	// we want to return IDs of the landmark that corresponds to each observation.

	for (int j = 0; j < observations.size(); j++) {
		int corresponding_id = -1; // initialize corresponding ID
		double min_distance = std::numeric_limits<double>::max(); // very high value for initial minimum distance
		for (unsigned int i = 0; i < predicted.size(); i++) {
			double distance = dist(observations[j].x, observations[j].y, predicted[i].x, predicted[i].y);
			if (distance < min_distance){
				min_distance = distance;
				corresponding_id = i;
			}
		}
		observations[j].id = corresponding_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
	std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. 
	// Ref: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. 
	// Particles are located according to the MAP'S coordinate system. 
	// We need to transform between the two systems.
	// Ref:
	// https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	// http://planning.cs.uiuc.edu/node99.html

	// standard deviations for observations
	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	// for each particle;
	for (unsigned int p = 0; p < num_particles; p++){
		// get particles properties
		double px = particles[p].x;
		double py = particles[p].y;
		double ptheta = particles[p].theta;

		// chose only landmarks within sensor range
		std::vector<LandmarkObs> landmarks_in_range;
		int landmark_size = map_landmarks.landmark_list.size();
		for (unsigned int mk = 0; mk < landmark_size; mk++){
			float lx = map_landmarks.landmark_list[mk].x_f;
			float ly = map_landmarks.landmark_list[mk].y_f;
			int l_id = map_landmarks.landmark_list[mk].id_i;

			// find landmark distance from particle
			double lndmrk_dist = dist(px, py, lx, ly);

			if (lndmrk_dist < sensor_range){
				landmarks_in_range.push_back({ l_id, lx, ly });
			}
		}

		// transform vehicle coordinate observations into map's coordinate
		std::vector<LandmarkObs> observations_m;
		for (unsigned int ob = 0; ob < observations.size(); ob++){
			double ox = observations[ob].x;
			double oy = observations[ob].y;

			double x = ox * cos(ptheta) - oy * sin(ptheta) + px;
			double y = ox * sin(ptheta) + oy * cos(ptheta) + py;
			observations_m.push_back({ observations[ob].id, x, y });
		}

		// perform data association
		dataAssociation(landmarks_in_range, observations_m);

		// re-initialize particle weight
		particles[p].weight = 1.0; 

		for (unsigned int i = 0; i < observations_m.size(); i++){
			int obs_index = observations_m[i].id;

			double o_x, o_y, obs_x, obs_y;
			o_x = observations_m[i].x;
			o_y = observations_m[i].y;

			double var_x = (landmarks_in_range[obs_index].x - o_x);
			double var_y = (landmarks_in_range[obs_index].y - o_y);
			double exponent = pow(var_x, 2) / (2.0 * pow(std_x, 2)) + pow(var_y, 2) / (2.0 * pow(std_y, 2));

			particles[p].weight *= (1.0 / (2.0*M_PI*std_x*std_y)) * exp(-exponent);
		}

	}
}

void ParticleFilter::resample() {
	// We resample particles with replacement with probability proportional to their weight. 
	// using std::discrete_distribution according to;
	// http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Produces random integers on the interval [0, n], 
	// where the probability of each individual integer i is defined as
	// w_(i) / S, 
	// that is the weight of the ith integer divided by the sum of all n weights.
	
	// integer random number generator
	std::random_device rd; 
	std::mt19937 gen(rd());

	// get all of the current weights
	std::vector<double> weights;
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}

	// create distribution
	std::discrete_distribution<> d(weights.begin(), weights.end());

	// placeholder for resampled particles
	std::vector<Particle> resampled_particles;
	resampled_particles.resize(num_particles);

	// sample new particles
	for (int n = 0; n < num_particles; n++) {
		int new_index = d(gen);
		resampled_particles[n] = particles[new_index];
	}
	particles = resampled_particles;
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
