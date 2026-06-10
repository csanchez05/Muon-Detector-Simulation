#include "DipoleField.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

DipoleField::DipoleField(G4double peakB, G4double magnetZLength) {
    fPeakB = peakB;
    fHalfLength = magnetZLength / 2.0;
    
    // Equivalent magnetic moment calculation to match fPeakB at the center
    // Assuming effective radius of a 2x2 inch block
    fEffectiveRadius = 25.4 * mm; 
    fDipoleMoment = fPeakB * std::pow(fEffectiveRadius, 3);
}

void DipoleField::GetFieldValue(const G4double Point[4], G4double* Bfield) const {
    G4double x = Point[0];
    G4double y = Point[1];
    G4double z = Point[2];
    
    // Distance from the center of the magnet
    G4double r2 = x*x + y*y + z*z;
    G4double r = std::sqrt(r2);

    // Inside the magnet gap: Uniform field in Y-direction
    if (std::abs(z) <= fHalfLength && std::abs(x) <= fEffectiveRadius && std::abs(y) <= fEffectiveRadius) {
        Bfield[0] = 0.0;
        Bfield[1] = fPeakB;
        Bfield[2] = 0.0;
    } 
    // Outside the magnet: 3D Vector Dipole Formulation
    else {
        // Avoid division by zero exactly at origin if logic fails
        if (r < 1e-6) r = 1e-6; 

        G4double r5 = r2 * r2 * r;
        
        // Dipole moment is aligned along the Y-axis: m = (0, m_y, 0)
        // B = (3*r*(m dot r) / r^5) - (m / r^3)
        // Here, (m dot r) = m_y * y
        G4double mDotR = fDipoleMoment * y;
        
        Bfield[0] = (3.0 * x * mDotR) / r5;
        Bfield[1] = ((3.0 * y * mDotR) / r5) - (fDipoleMoment / std::pow(r, 3));
        Bfield[2] = (3.0 * z * mDotR) / r5;
    }
}