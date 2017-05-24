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

//using namespace std;

// declare random engine:
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Number of particles
	num_particles = 100;
	particles.resize(num_particles);

	// Create a normal (Gaussian) distribution for x, y and theta.
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
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

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
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// predicted(here) is observation
	// observation(here) is landmarks in range
	// we want to return ids of the landmark that corresponds to the observation

	for (int j = 0; j < observations.size(); j++) {
		int choice_id = -1;
		double min_distance = std::numeric_limits<double>::max(); // very a high value for initial minimum distance
		for (unsigned int i = 0; i < predicted.size(); i++) {
			double distance = dist(observations[j].x, observations[j].y, predicted[i].x, predicted[i].y);
			if (distance < min_distance){
				min_distance = distance;
				choice_id = predicted[i].id;
			}
		}
		observations[j].id = choice_id;
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

	// for each particle
	for (int p = 0; p < num_particles; p++){
		// get particles properties
		double px = particles[p].x;
		double py = particles[p].y;
		double ptheta = particles[p].theta;

		// chose landmarks within sensor range
		std::vector<LandmarkObs> landmarks_in_range;
		int landmark_size = map_landmarks.landmark_list.size();
		for (unsigned int mk = 0; mk < landmark_size; mk++){
			float lx = map_landmarks.landmark_list[mk].x_f;
			float ly = map_landmarks.landmark_list[mk].y_f;
			int l_id = map_landmarks.landmark_list[mk].id_i;

			//std::cout << "lx :" << lx << " ly: " << ly << " l_id: " << l_id << std::endl;

			double lndmrk_dist = dist(px, py, lx, ly);
			if (lndmrk_dist < sensor_range){
				landmarks_in_range.push_back({ l_id, lx, ly });
			}

		}

		// transform vehicle coordinate observations into map coordinates
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

		double std_x = std_landmark[0];
		double std_y = std_landmark[1];

		particles[p].weight = 1.0;
		for (unsigned int i = 0; i < observations_m.size(); i++){
			int obs_index = observations_m[i].id;

			double o_x, o_y, obs_x, obs_y;
			o_x = observations_m[i].x;
			o_y = observations_m[i].y;

			// get the x,y coordinates of the prediction associated with the current observation
			for (unsigned int k = 0; k < landmarks_in_range.size(); k++) {
				if (landmarks_in_range[k].id == obs_index) {
					obs_x = landmarks_in_range[k].x;
					obs_y = landmarks_in_range[k].y;
				}
			}

			double var_x = (obs_x - o_x);
			double var_y = (obs_y - o_y);
			double exponent = pow(var_x, 2) / (2.0 * pow(std_x, 2)) + pow(var_y, 2) / (2.0 * pow(std_y, 2));

			particles[p].weight *= (1.0 / (2.0*M_PI*std_x*std_y)) * exp(-exponent);
		}

	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	std::random_device rd;
	std::mt19937 gen(rd());
	std::vector<Particle> new_particles;

	// get all of the current weights
	std::vector<double> weights;
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}

	// generate random starting index for resampling wheel
	std::uniform_int_distribution<int> uniintdist(0, num_particles - 1);
	auto index = uniintdist(gen);

	// get max weight
	double max_weight = *max_element(weights.begin(), weights.end());

	// uniform random distribution [0.0, max_weight)
	std::uniform_real_distribution<double> unirealdist(0.0, max_weight);

	double beta = 0.0;

	// spin the resample wheel!
	for (int i = 0; i < num_particles; i++) {
		beta += unirealdist(gen) * 2.0;
		while (beta > weights[index]) {
			beta -= weights[index];
			index = (index + 1) % num_particles;
			
		}
		new_particles.push_back(particles[index]);
		//std::cout << "chosen index: "<< index << std::endl;
	}

	//std::cout << "number of particles new: "<<particles.size() << std::endl;
	particles = new_particles;

	//std::discrete_distribution<> d(weights.begin(), weights.end());

	//std::cout << "Number of particles \t" << num_particles << std::endl;
	//
	//
	//std::cout << "Distribution \t" << d << std::endl;

	//std::vector<Particle> resampled_particles;
	//resampled_particles.resize(num_particles);

	//for (int n = 0; n < num_particles; n++) {

	//	// cout << "Resample step: \t" << n << endl;
	//	int new_index = d(gen);
	//	 std::cout << new_index << std::endl ;
	//	resampled_particles[n] = particles[new_index];
	//}
	//// cout << "Resample complete" << endl;
	//particles = resampled_particles;
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
