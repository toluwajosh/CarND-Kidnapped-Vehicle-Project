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

// random engine:
std::default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//std::cout << "original:: x: " << x << " y: " << y << " theta: " << theta << std::endl;

	// Number of particles
	num_particles = 50;
	particles.resize(num_particles);

	// Create a normal (Gaussian) distribution for x, y and theta.
	
	std::normal_distribution<double> dist_x(0, std[0]);
	std::normal_distribution<double> dist_y(0, std[1]);
	std::normal_distribution<double> dist_theta(0, std[2]);

	// Initialize all particles
	//std::cout << "\nSamples:\n" << std::endl;
	for (int i = 0; i < num_particles; i++)
	{
		particles[i].id = i;
		particles[i].x = x + dist_x(gen);
		particles[i].y = y + dist_y(gen);
		particles[i].theta = theta + dist_theta(gen);

		// set weights to 1
		particles[i].weight = 1.0;

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

	for (int i = 0; i < num_particles; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		
		//std::cout << "Before: x: " << x << " y: " << y << " theta: " << theta << std::endl;
		
		const double v_yaw_rate = velocity / yaw_rate;
		const double yaw_delta_t = yaw_rate * delta_t;
		const double v_dt = velocity * delta_t;

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

		// Noise
		std::normal_distribution<double> dist_x(0, std_pos[0]);
		std::normal_distribution<double> dist_y(0, std_pos[1]);
		std::normal_distribution<double> dist_theta(0, std_pos[2]);

		particles[i].id = i;
		particles[i].x = x + dist_x(gen);
		particles[i].y = y + dist_y(gen);
		particles[i].theta = theta + dist_theta(gen);
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

	// predicted is observation
	// observation is landmarks in range
	// we want to return ids of the landmark that corresponds to the observation
	std::vector<LandmarkObs> new_observation;
	const double threshold = 10;

	for (int j = 0; j < observations.size(); j++) {
		int choice_id = 0;
		double min_distance = std::numeric_limits<double>::max(); // very a high value for initial minimum distance
		for (int i = 0; i < predicted.size(); i++) {
			double distance = dist(observations[j].x, observations[j].y, predicted[i].x, predicted[i].y);
			if (distance < min_distance){
				min_distance = distance;
				choice_id = predicted[i].id;
				//std::cout << "distance: " << distance << std::endl;
			}
		}
		observations[j].id = choice_id;
		//new_observation.push_back(observations[j]);
	}

	//observations = new_observation;

	//// for debug
	//for (int i = 0; i < observations.size(); i++) {
	//	
	//	double x1, y1, x2, y2;
	//	x1 = observations[i].x;
	//	y1 = observations[i].y;
	//	x2 = predicted[observations[i].id].x;
	//	y2 = predicted[observations[i].id].y;
	//	double distance = dist(x1, y1, x2, y2);
	//	std::cout <<"["<< x1<<",\t"<<y1<<"]" << " \tcorresponds to: " << "[" << x2 << ",\t" << y2 << "]" << " \twith distance: \t" << distance << std::endl;

	//}
	//std::cout << "\n " << std::endl;
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

	const double std_x = std_landmark[0];
	const double std_y = std_landmark[1];

	long double weight_sum = 0;
	for (int p = 0; p < num_particles; p++){

		const double px = particles[p].x;
		const double py = particles[p].y;
		const double theta = particles[p].theta;// *M_PI / 180.0;

		// Transform observations to map cordinate
		for (int i = 0; i < observations.size(); i++) {
			const double ox = observations[i].x;
			const double oy = observations[i].y;
			
			//double theta = particles[p].theta;
			observations[i].x = ox * cos(theta) - oy * sin(theta) + px;
			observations[i].y = ox * sin(theta) + oy * cos(theta) + py;
			observations[i].id = i;
			//std::cout << observations[i].id << std::endl;
		}

		// use landmarks within sensor range
		std::vector<LandmarkObs> lndmks_in_range;
		LandmarkObs l_obs;

		int no_of_lndmrks = map_landmarks.landmark_list.size();
		for (int l = 0; l < no_of_lndmrks; l++) {

			const double lndmrk_x = map_landmarks.landmark_list[l].x_f;
			const double lndmrk_y = map_landmarks.landmark_list[l].y_f;
			const double lndmrk_id = map_landmarks.landmark_list[l].id_i;

			// distance of landmark from the particle
			const double lndmrk_dist = dist(px, py, lndmrk_x, lndmrk_y);

			if (lndmrk_dist < sensor_range){
				l_obs.id = lndmrk_id;
				l_obs.x = lndmrk_x;
				l_obs.y = lndmrk_y;

				lndmks_in_range.push_back(l_obs);
			}
		}

		//std::cout << "number of landmarks: " << lndmks_in_range.size() << std::endl;
		//std::cout << "number of observations: " << observations.size() << std::endl;

		// Make data association
		dataAssociation(observations, lndmks_in_range);

		//// for debug: check for correspondence after data association
		//for (int i = 0; i < observations.size(); i++) {
		//	
		//	double x1, y1, x2, y2;
		//	x1 = lndmks_in_range[i].x;
		//	y1 = lndmks_in_range[i].y;
		//	x2 = observations[lndmks_in_range[i].id].x;
		//	y2 = observations[lndmks_in_range[i].id].y;
		//	double distance = dist(x1, y1, x2, y2);
		//	std::cout << observations[i].id << std::endl;
		//	std::cout <<"["<< x1<<",\t"<<y1<<"]" << " \tcorresponds to: " << "[" << x2 << ",\t" << y2 << "]" << " \twith distance: \t" << distance << std::endl;

		//}
		//std::cout << "\n " << std::endl;

		// Update weight
		//long double weight = 1.0; // reinitialize weight
		particles[p].weight = 1.0;
		for (int i = 0; i < lndmks_in_range.size(); i++) {
			//std::cout << "weight before: " << weight << std::endl;
			//std::cout << lndmks_in_range[i].x <<" corresponds to: " << observations[lndmks_in_range[i].id].x << std::endl;
			int obs_index = lndmks_in_range[i].id;

			double obs_x;
			double obs_y;


			// get the x,y coordinates of the prediction associated with the current observation
			for (unsigned int k = 0; k < observations.size(); k++) {
				if (observations[k].id == obs_index) {
					obs_x = observations[k].x;
					obs_y = observations[k].y;
				}
			}

			double var_x = (obs_x - lndmks_in_range[i].x);
			double var_y = (obs_y - lndmks_in_range[i].y);
			double exponent = pow(var_x,2) / (2.0 * pow(std_x,2)) + pow(var_y,2) / (2.0 * pow(std_y,2));

			//std::cout << "var_x: " << var_x << std::endl;
			//std::cout << "var_y: " << var_y << std::endl;
			//std::cout << "exponent: " << exponent << std::endl;

			///////////////
			//if (exponent < 350)
			//{
			//	double inter_weight = (1.0 / (2.0*M_PI*std_x*std_y)) * exp(-exponent);
			//	//std::cout << "inter weight: " << inter_weight << "\n" << std::endl;
			//	weight *= inter_weight;
			//}

			////////////
			//const long double inter_weight = (1.0 / (2.0*M_PI*std_x*std_y)) * exp(-exponent);
			//std::cout << "inter weight: " << inter_weight << "\n" << std::endl;
			particles[p].weight *= (1.0 / (2.0*M_PI*std_x*std_y)) * exp(-exponent);

		}
		//particles[p].weight = weight;
		//weight_sum += particles[p].weight;
		////std::cout<<"landmarks observations: "<< observations.size() << std::endl;
		//weights.push_back(particles[p].weight);
	}

	//// normalize weights
	//for (int i = 0; i < particles.size(); i++) {
	//	//std::cout << "weight: " << particles[i].weight << std::endl;
	//	weights[i] /= weight_sum;
	//}

	//for (int i = 0; i < num_particles; i++) {
	//std::cout << "Weights \t" << particles[i].weight << std::endl;
	//}
	//std::cout << "\n " << std::endl;
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
		std::cout << index << std::endl;
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


	std::cout << "\n " << std::endl;
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
