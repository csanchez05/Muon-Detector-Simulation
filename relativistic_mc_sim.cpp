#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>

struct Particle {
    double x, y, z;
    double vx, vy, vz;
    double m, q;
};

int main() {
    // 1. SETUP RANDOMNESS
    std::random_device rd;
    std::mt19937 gen(rd());
    
    // Position Jitter (5cm standard deviation)
    std::normal_distribution<double> pos_jitter(0.0, 0.05);

    // Speed Variation (0.75c to 0.999c)
    double c = 2.9979e8; 
    std::uniform_real_distribution<double> speed_dist(0.75 * c, 0.999 * c);

    // Coin Flip (50% chance True, 50% chance False)
    std::bernoulli_distribution coin_flip(0.5);

    // 2. CONSTANTS
    double m_rest = 200*9.109e-31;
    double q_mag = 1.602e-19; // Magnitude of charge
    
    // Tiny time step for high accuracy
    double dt = 1e-13; 
    double screen_dist = 1; 

    // 3. INPUTS
    double Ey, Bz;
    std::cout << "--- Dual-Charge Relativistic Simulation ---" << std::endl;
    std::cout << "Enter Electric Field (Ey) (V/m) [Try 1000000]: ";
    std::cin >> Ey;
    std::cout << "Enter Magnetic Field (Bz) (T) [Try 0.01]: ";
    std::cin >> Bz;

    // 4. FILE SETUP
    std::ofstream file;
    file.open("thomson_mc_rel.csv");
    // Header now includes 'q' so we can color-code by charge
    file << "id,x,y,z,gamma,q" << std::endl; 

    int num_particles = 1000;

    std::cout << "Simulating " << num_particles << " particles..." << std::endl;

    for (int p_id = 0; p_id < num_particles; p_id++) {
            
        Particle p;
        p.m = m_rest;

        // --- THE 50/50 CHARGE FLIP ---
        // If coin_flip is true, use -q (Electron). Else use +q (Positron).
        if (coin_flip(gen)) {
            p.q = -q_mag; 
        } else {
            p.q = q_mag;
        }
        
        // --- RANDOMIZE VELOCITY ---
        p.vx = speed_dist(gen); 
        p.vy = 0.0; 
        p.vz = 0.0; 

        // --- RANDOMIZE POSITION ---
        p.x = 0;
        p.y = pos_jitter(gen);
        p.z = pos_jitter(gen);

        while (p.x < screen_dist) {

            // --- A. RELATIVISTIC FACTORS ---
            double v_sq = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
            
            // Safety Clamp
            if (v_sq >= c*c) v_sq = c*c * 0.999999;

            double gamma = 1.0 / std::sqrt(1.0 - v_sq / (c*c));
            double m_rel = p.m * gamma;

            // --- B. PHYSICS ENGINE ---
            double Fe_y = p.q * Ey;
            double Fm_x = p.q * (p.vy * Bz);
            double Fm_y = -p.q * (p.vx * Bz);

            // Acceleration uses Relativistic Mass
            double ax = Fm_x / m_rel;
            double ay = (Fe_y + Fm_y) / m_rel;

            p.vx += ax * dt;
            p.vy += ay * dt;    

            p.x += p.vx * dt;
            p.y += p.vy * dt;

            // Wall check
            if (std::abs(p.y) > 0.5) break; 
            
            // --- C. FAST OUTPUT ---
            // Using "\n" instead of std::endl for speed
            file << p_id << "," << p.x << "," << p.y << "," << p.z << "," << gamma << "," << p.q << "\n";
        }
        
        if (p_id % 10 == 0) {
            std::cout << "Particle " << p_id << " done." << std::endl;
        }
    }
    file.close();
    std::cout << "Done. Saved to thomson_mc_rel.csv" << std::endl;
    return 0;
}