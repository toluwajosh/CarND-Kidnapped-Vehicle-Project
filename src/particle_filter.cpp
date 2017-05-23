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
	num_particles = 50;
	particles.resize(num_particles);

	// Standard deviations for x, y, and theta
	double std_x, std_y, std_theta;
	std_x = std[0];
	std_y = std[1];
	std_theta = std[2];
	//

	// Create a normal (Gaussian) distribution for x, y and theta.
	std::default_random_engine gen;
	std::normal_distribution<double> dist_x(0, std_x);
	std::normal_distribution<double> dist_y(0, std_y);
	std::normal_distribution<double> dist_theta(0, std_theta);

	// Initialize all particles
	for (int i = 0; i < num_particles; i++)
	{
		particles[i].id = i;
		particles[i].x = x + dist_x(gen);
		particles[i].y = y + dist_y(gen);
		particles[i].theta = theta + dist_theta(gen);
		
		// Initialize weights to 1
		particles[i].weight = 1;

		//std::cout << "x: " << particles[i].x << " y: " << particles[i].y << " theta: " << particles[i].theta << std::endl;
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std::default_random_engine gen;

	for (int i = 0; i < num_particles; i++)
	{
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta * std::_Pi / 180.0;

		const double v_yaw_rate = velocity / yaw_rate;
		const double yaw_delta_t = yaw_rate * delta_t;
		const double v_dt = velocity * delta_t;

		//
		//std::cout << "Before: x: " << x << " y: " << y << " theta: " << theta << std::endl;

		// Noise
		std::normal_distribution<double> dist_x(0, std_pos[0]);
		std::normal_distribution<double> dist_y(0, std_pos[1]);
		std::normal_distribution<double> dist_theta(0, std_pos[2]);
		if (yaw_rate>1e-3){
			x += dist_x(gen) + v_yaw_rate*(sin(theta + yaw_delta_t) - sin(theta));
			y += dist_y(gen) + v_yaw_rate*(cos(theta) - cos(theta + yaw_rate * delta_t));
			theta += dist_theta(gen) + yaw_delta_t;
		}
		else
		{
			x += dist_x(gen) + v_dt*sin(theta);
			y += dist_y(gen) + v_dt*cos(theta);
			theta += dist_theta(gen);
		}
		particles[i].x = x;
		particles[i].y = y;
		particles[i].theta = theta;
		//std::cout << "After: x: " << particles[i].x << " y: " << particles[i].y << " theta: " << particles[i].theta << std::endl;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	//for (int i = 0; i < predicted.size(); i++) {
	//	std::cout << "observations: " << predicted[i].x << std::endl;
	//}

	//for (int i = 0; i < observations.size(); i++) {
	//	std::cout << "landmarks: " << observations[i].x << std::endl;
	//}

	for (int j = 0; j < observations.size(); j++) {
		int choice_id = 0;
		double old_distance = 1000; // a high value for initial distance
		for (int i = 0; i < predicted.size(); i++) {
			double distance = dist(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
			if (distance < old_distance){
				old_distance = distance;
				choice_id = i;
				//std::cout << "distance: " << distance << std::endl;
			}
		}
		observations[j].id = choice_id;
	}
	
	//// for debug
	//for (int i = 0; i < observations.size(); i++) {
	//	std::cout << observations[i].x << " corresponds to: " << predicted[observations[i].id].x << std::endl;
	//	double x1, y1, x2, y2;
	//	x1 = observations[i].x;
	//	y1 = observations[i].y;
	//	x2 = predicted[observations[i].id].x;
	//	y2 = predicted[observations[i].id].y;
	//	double distance = dist(x1, y1, x2, y2);
	//	std::cout << "Distance between: " << distance << std::endl;
	//}
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

	double std_x = std_landmark[0];
	double std_y = std_landmark[1];

	for (int p = 0; p < num_particles; p++) {
		
		// Transform observations to map cordinate
		for (int i = 0; i < observations.size(); i++) {
			double x = observations[i].x;
			double y = observations[i].y;
			double xt = particles[p].x;
			double yt = particles[p].y;
			double theta = particles[p].theta * std::_Pi / 180.0;
			//double theta = particles[p].theta;
			observations[i].x = x * cos(theta) - y * sin(theta) + xt;
			observations[i].y = x * sin(theta) + y * cos(theta) + yt;
		}

		// use landmarks within sensor range
		std::vector<LandmarkObs> lndmks_in_range;
		LandmarkObs obc;
		//std::cout << "number of landmarks: " << map_landmarks.landmark_list.size() << std::endl;
		int no_of_lndmrks = map_landmarks.landmark_list.size();
		for (int l = 0; l < no_of_lndmrks; l++){
			double lndmrk_dist = dist(particles[p].x, particles[p].y, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].x_f);
			if (lndmrk_dist < sensor_range) {
				obc.x = map_landmarks.landmark_list[l].x_f;
				obc.y = map_landmarks.landmark_list[l].y_f;
				lndmks_in_range.push_back(obc);
			}
		}
		//std::cout << "number of landmarks: " << lndmks_in_range.size() << std::endl;

		// Make data association
		dataAssociation(observations, lndmks_in_range);

		//for (int i = 0; i < observations.size(); i++) {
		//	std::cout << "observations: " << observations[i].x << std::endl;
		//}

		//for (int i = 0; i < lndmks_in_range.size(); i++) {
		//	std::cout << "landmarks: " << lndmks_in_range[i].x << std::endl;
		//}

		// Update weight
		double weight = particles[p].weight;
		for (int i = 0; i < lndmks_in_range.size(); i++) {
			//std::cout << "weight before: " << weight << std::endl;
			//std::cout << lndmks_in_range[i].x <<" corresponds to: " << observations[lndmks_in_range[i].id].x << std::endl;
			double var_x = (lndmks_in_range[i].x - observations[lndmks_in_range[i].id].x)*(lndmks_in_range[i].x - observations[lndmks_in_range[i].id].x);
			double var_y = (lndmks_in_range[i].y - observations[lndmks_in_range[i].id].y)*(lndmks_in_range[i].y - observations[lndmks_in_range[i].id].y);
			double distdist = var_x + var_y;

			//std::cout << "var_x: " << var_x << std::endl;
			//std::cout << "var_y: " << var_y << std::endl;
			//std::cout << "distdist: " << distdist << std::endl;
			if (distdist<50)
			{
				double inter_weight = (1.0 / (2.0 * std::_Pi*std_x*std_y))*exp(-(var_x / (2.0 * std_x*std_x)) - (var_y / (2.0 * std_y*std_y)));

				//std::cout << "inter weight: " << inter_weight << "\n" << std::endl;
				weight *= inter_weight;
			}
			//std::cout << "weight after: " << weight <<"\n"<< std::endl;
		}
		particles[p].weight = weight;
		//std::cout<<"landmarks observations: "<< observations.size() << std::endl;
	}

	// Normalize weights
	double weight_sum = 0;
	for (int i = 0; i < particles.size(); i++){
		weight_sum += particles[i].weight;
	}
	for (int i = 0; i < particles.size(); i++){
		particles[i].weight /= weight_sum;
		std::cout << "weight: " << particles[i].weight << std::endl;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::random_device rd;
	std::mt19937 gen(rd());

	int index = std::rand() % particles.size();
	double beta = 0.0;
	double max_weight = 0.0;
	for (int i = 0; i < particles.size(); i++){
		if (particles[i].weight>max_weight){
			max_weight = particles[i].weight;
		}
	}

	for (int i = 0; i < particles.size(); i++) {
		beta += (std::rand() / 32767.0) *2.0 * max_weight;
		while (beta>particles[index].weight){
			beta -= beta > particles[index].weight;
			index = (index + 1) % particles.size();
		}
		particles[i].x = particles[index].x;
		particles[i].y = particles[index].y;
		particles[i].weight = particles[index].weight;
	}
	//std::cout << "random number: " << std::rand() << std::endl;

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
