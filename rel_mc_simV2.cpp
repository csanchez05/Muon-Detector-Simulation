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
    std::normal_distribution<double> pos_jitter(0.0, 0.02);

    // Speed Variation (0.75c to 0.999c)
    double c = 2.9979e8; 
    std::uniform_real_distribution<double> speed_dist(0.994 * c, 0.998 * c);

    // Coin Flip (50% chance True, 50% chance False)
    std::bernoulli_distribution coin_flip(0.5);

    // 2. CONSTANTS
    double m_rest = 200 * 9.109e-31; // Muon Mass
    double q_mag = 1.602e-19; 
    
    double dt = 1e-13; 
    
    // GEOMETRY SETUP
    double screen_dist = 1.2;   // The screen is far away (1 meter)
    double magnet_width = 0.20; // The magnet is only 20cm long

    // 3. INPUTS
    double Ey_input = 0.0;
    double Bz_input;
    std::cout << "--- Finite Magnet Simulation ---" << std::endl;
    // std::cout << "Enter Electric Field (Ey) (V/m): ";
    // std::cin >> Ey_input;
    std::cout << "Enter Magnetic Field (Bz) (T): ";
    std::cin >> Bz_input;

    // 4. FILE SETUP
    std::ofstream file;
    file.open("thomson_mc_relV2.csv");
    file << "id,x,y,z,gamma,q,vx,vy,vz" << std::endl;

    int num_particles = 1250;

    std::cout << "Simulating " << num_particles << " particles..." << std::endl;

    for (int p_id = 0; p_id < num_particles; p_id++) {
            
        Particle p;
        p.m = m_rest;

        // Charge Flip
        if (coin_flip(gen)) p.q = -q_mag; 
        else p.q = q_mag;
        
        // Random Velocity & Position
        p.vx = speed_dist(gen); 
        p.vy = 0.0; 
        p.vz = 0.0; 

        p.x = 0;
        p.y = pos_jitter(gen);
        p.z = pos_jitter(gen);

        while (p.x < screen_dist) {

            // --- A. RELATIVISTIC FACTORS ---
            double v_sq = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
            if (v_sq >= c*c) v_sq = c*c * 0.999999;
            double gamma = 1.0 / std::sqrt(1.0 - v_sq / (c*c));
            double m_rel = p.m * gamma;

            // --- B. FIELD REGION CHECK (The "Wall") ---
            double current_Ey = 0.0;
            double current_Bz = 0.0;

            // Only apply fields if we are inside the magnet (x < 0.20)
            if (p.x <= magnet_width) {
                current_Ey = Ey_input;
                current_Bz = Bz_input;
            } 
            // If p.x > 0.20, fields remain 0.0 (Drift Region)
            
            // --- C. PHYSICS ENGINE ---
            // Use 'current_' variables, not global inputs
            double Fe_y = p.q * current_Ey;
            double Fm_x = p.q * (p.vy * current_Bz);
            double Fm_y = -p.q * (p.vx * current_Bz);

            double ax = Fm_x / m_rel;
            double ay = (Fe_y + Fm_y) / m_rel;

            p.vx += ax * dt;
            p.vy += ay * dt;    

            p.x += p.vx * dt;
            p.y += p.vy * dt;

            // Wall check (Increased to 1.0m since screen is further back)
            if (std::abs(p.y) > 1.2) break; 
            
            // Fast Output
            file << p_id << "," << p.x << "," << p.y << "," << p.z << "," << gamma << "," << p.q << "," << p.vx << "," << p.vy << "," << p.vz << "\n";
        }
        
        if (p_id % 50 == 0) {
            std::cout << "Particle " << p_id << " done." << std::endl;
        }
    }
    file.close();
    std::cout << "Done. Saved to thomson_mc_relV2.csv" << std::endl;
    return 0;
}