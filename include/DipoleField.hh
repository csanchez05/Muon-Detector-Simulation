#ifndef DipoleField_h
#define DipoleField_h 1

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

class DipoleField : public G4MagneticField {
public:
    // Constructor takes the peak field strength and the physical Z-length of the magnet
    DipoleField(G4double peakB, G4double magnetZLength);
    ~DipoleField() override = default;

    // This is the core engine: Geant4 calls this at every single tracking step
    void GetFieldValue(const G4double Point[4], G4double* Bfield) const override;

private:
    G4double fPeakB;
    G4double fHalfLength;
    G4double fEffectiveRadius;
    G4double fDipoleMoment;
};

#endif // DipoleField_h