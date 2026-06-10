#ifndef COMSOL_MAGNETIC_FIELD_HH
#define COMSOL_MAGNETIC_FIELD_HH

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <array>
#include <string>
#include <vector>

class COMSOLMagneticField : public G4MagneticField {
public:
    explicit COMSOLMagneticField(const std::string& filename);
    ~COMSOLMagneticField() override = default;

    void GetFieldValue(const G4double point[4], G4double* Bfield) const override;

private:
    std::size_t Index(std::size_t ix, std::size_t iy, std::size_t iz) const;

    std::size_t NearestIndex(const std::vector<G4double>& grid,
                             G4double value,
                             const char* axisName) const;

    bool FindCell(const std::vector<G4double>& grid,
                  G4double value,
                  std::size_t& i0,
                  G4double& t) const;

    std::array<G4double, 3> FieldAt(std::size_t ix,
                                    std::size_t iy,
                                    std::size_t iz) const;

private:
    std::vector<G4double> fX;
    std::vector<G4double> fY;
    std::vector<G4double> fZ;

    // flattened array: ((ix * ny + iy) * nz + iz)
    std::vector<std::array<G4double, 3>> fB;

    G4double fXMin;
    G4double fXMax;
    G4double fYMin;
    G4double fYMax;
    G4double fZMin;
    G4double fZMax;
};

#endif