/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles
  std::random_device mch;
  std::default_random_engine gen(mch());

  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  this->particles.clear();

  std::cout << "Standard Deviations for GPS: {" << std[0] << ", " << std[1] << ", " << std[2] << "}" << std::endl;
  std::cout << "Given GPS Position (x, y, theta) : (" << x << ", " << y  << ", " << theta << ")" << std::endl;

  float mean_x = 0;
  float mean_y = 0;

  for (int i = 0; i < num_particles; i++) {
     Particle p = {};
     p.id = i;
     p.x = dist_x(gen);
     p.y = dist_y(gen);
     p.theta = dist_theta(gen);
     p.weight = 1;
     this->particles.push_back(p);
     mean_x += p.x;
     mean_y += p.y;
  }
  this->is_initialized = true;
  std::cout << "Initialised particles mean (x, y) : (" << mean_x / num_particles << ", " << mean_y / num_particles  << ")" << std::endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::random_device mch;
  std::default_random_engine gen(mch());
  std::normal_distribution<float> dist_x(0.0, std_pos[0]);
  std::normal_distribution<float> dist_y(0.0, std_pos[1]);
  std::normal_distribution<float> dist_theta(0.0, std_pos[2]);
  float x, y, theta; // current particle position;
  float xp, yp, theta_p, dx, err_x; // predicted particle position;
  for (int i = 0; i < num_particles; i++) {
     x = this->particles[i].x;
     theta = this->particles[i].theta;
     if (yaw_rate == 0.0) {
        yaw_rate = 0.0001;
     }
     dx = velocity / yaw_rate * (sin(theta + yaw_rate * delta_t)  - sin(theta));
     xp = x + dx;
     err_x = dist_x(gen);
     xp += err_x;

     y = this->particles[i].y;
     yp = y + velocity / yaw_rate * (cos(theta)  - cos(theta  + yaw_rate * delta_t));
     yp += dist_y(gen);

     theta_p = theta + yaw_rate * delta_t;
     theta_p = theta_p + dist_theta(gen);

     particles[i].x = xp;
     particles[i].y = yp;
     particles[i].theta = theta_p;

   }

  std::cout << "Speed & Yaw, & time: (" << velocity << ", " << yaw_rate << ", " << delta_t << ")" << std::endl;
  std::cout << "Moved (x, dx, err_x, y, theta) -> (p1) is (" << x << ", " << dx << ", " << err_x << ", " << y << ", " << theta << ") -> (" << xp << ", " << yp << ", " << theta_p << ")" << std::endl;

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

   float distance = 0.0;
   float min_distance;
   for (int i = 0; i < observations.size(); i++) {
      min_distance = -1;
      for (int j = 0; j < predicted.size(); j++) {

         distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
         if (distance < min_distance || min_distance < 0) {
            min_distance = distance;
            // link observation to the nearest landmark
            observations[i].id = j;
         }
      }
   }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

   float x_sigma = std_landmark[0];
   float y_sigma = std_landmark[1];

   float prob; // probability of the observation-landmark match
   float total_weight = 0;

   float px; // particle x
   float py; // particle y
   float ptheta; // particle theta

   float x, y; // observation coordinates
   float mx, my; // position of the landmark in global coordinates
   float p; // probability of the landmark / observation combination


   vector<LandmarkObs> landmarks_predicted;
   vector<LandmarkObs> observs;



   for (int i = 0; i < particles.size(); i ++) {
     px = particles[i].x;
     py = particles[i].y;
     ptheta = particles[i].theta;


     for (int j = 0; j < observations.size(); j ++) {
       LandmarkObs observation = {}; // observation struct in global coordinates
       observation.x = px + observations[j].x * cos(ptheta) - observations[j].y * sin(ptheta);
       observation.y = py + observations[j].x * sin(ptheta) + observations[j].y * cos(ptheta);
       observation.id = observations[j].id;
       observs.push_back(observation);
     }

     // Vector of map landmarks limited to sensor_range
     for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
        /* Going through all landmarks and transforming them to car coordinates */
        LandmarkObs landmark = {}; // individual landmark with position in local coordinates;
        landmark.id = map_landmarks.landmark_list[j].id_i;
        landmark.x = map_landmarks.landmark_list[j].x_f;
        landmark.y = map_landmarks.landmark_list[j].y_f;
        if (dist(px, py, landmark.x, landmark.y) <= sensor_range * 1.1) {
           landmarks_predicted.push_back(landmark);
        }
     }

     // Assign landmarks to the observations
     this->dataAssociation(landmarks_predicted, observs);

     // For each particle calculate probability of the observations
     prob = 1;
     for (int k = 0; k < observs.size(); k++) {
        // for each observation associated with landmark calculate probability
        // of that measurement
        x = observs[k].x;
        y = observs[k].y;

        mx = landmarks_predicted[observs[k].id].x;
        my = landmarks_predicted[observs[k].id].y;

        p = 1 / (2 * x_sigma * y_sigma * M_PI) * exp(-((pow(x-mx, 2)) / 2 / pow(x_sigma, 2) + (pow(y-my, 2)) / 2 / pow(y_sigma, 2)));
        prob *= p;
     }

     // Assign probability of all observations
     particles[i].weight = prob;

     // Calculate sum of weights of particles for normalisation
     total_weight += prob;
     landmarks_predicted.clear();
     observs.clear();
   }

   for (int i = 0; i < particles.size(); i++) {
      particles[i].weight /= total_weight;
   }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
   std::random_device mch;
   std::default_random_engine gen(mch());
   std::vector<float> weights;
   std::vector<Particle> new_particles;

   // initialise beta as the largest weight of all particles
   for (int j = 0; j < particles.size(); j++) {
      weights.push_back(particles[j].weight);
   }

   std::discrete_distribution<int> d(weights.begin(), weights.end());

   int i = 0;
   int j = 0;

   while (i < num_particles) {
      j = d(gen);
      new_particles.push_back(particles[j]);
      i += 1;
   }

   this->particles = new_particles; // update filter with the new set of particles
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
