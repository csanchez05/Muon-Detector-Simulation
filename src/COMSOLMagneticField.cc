#include "COMSOLMagneticField.hh"

#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

constexpr G4double kCoordinateTolerance = 1.0e-9 * mm;

bool IsCommentOrEmpty(const std::string& line) {
    if (line.empty()) return true;

    for (char c : line) {
        if (std::isspace(static_cast<unsigned char>(c))) continue;
        return c == '%';
    }

    return true;
}

}

COMSOLMagneticField::COMSOLMagneticField(const std::string& filename) {
    std::ifstream input(filename);

    if (!input) {
        throw std::runtime_error("COMSOLMagneticField: could not open field file: " + filename);
    }

    struct Row {
        G4double x;
        G4double y;
        G4double z;
        G4double bx;
        G4double by;
        G4double bz;
    };

    std::vector<Row> rows;
    std::set<G4double> xSet;
    std::set<G4double> ySet;
    std::set<G4double> zSet;

    std::string line;

    while (std::getline(input, line)) {
        if (IsCommentOrEmpty(line)) continue;

        std::istringstream iss(line);

        G4double x_mm = 0.0;
        G4double y_mm = 0.0;
        G4double z_mm = 0.0;
        G4double bx_T = 0.0;
        G4double by_T = 0.0;
        G4double bz_T = 0.0;
        G4double normB_T = 0.0;

        if (!(iss >> x_mm >> y_mm >> z_mm >> bx_T >> by_T >> bz_T >> normB_T)) {
            throw std::runtime_error("COMSOLMagneticField: failed to parse line:\n" + line);
        }

        Row row;
        row.x = x_mm * mm;
        row.y = y_mm * mm;
        row.z = z_mm * mm;

        row.bx = bx_T * tesla;
        row.by = by_T * tesla;
        row.bz = bz_T * tesla;

        rows.push_back(row);

        xSet.insert(row.x);
        ySet.insert(row.y);
        zSet.insert(row.z);
    }

    if (rows.empty()) {
        throw std::runtime_error("COMSOLMagneticField: field file has no numeric rows: " + filename);
    }

    fX.assign(xSet.begin(), xSet.end());
    fY.assign(ySet.begin(), ySet.end());
    fZ.assign(zSet.begin(), zSet.end());

    if (fX.size() < 2 || fY.size() < 2 || fZ.size() < 2) {
        throw std::runtime_error("COMSOLMagneticField: grid must have at least 2 points along each axis.");
    }

    fXMin = fX.front();
    fXMax = fX.back();
    fYMin = fY.front();
    fYMax = fY.back();
    fZMin = fZ.front();
    fZMax = fZ.back();

    const std::size_t nTotal = fX.size() * fY.size() * fZ.size();
    fB.assign(nTotal, {0.0, 0.0, 0.0});

    std::vector<bool> filled(nTotal, false);

    for (const auto& row : rows) {
        const std::size_t ix = NearestIndex(fX, row.x, "x");
        const std::size_t iy = NearestIndex(fY, row.y, "y");
        const std::size_t iz = NearestIndex(fZ, row.z, "z");

        const std::size_t idx = Index(ix, iy, iz);

        fB[idx] = {row.bx, row.by, row.bz};
        filled[idx] = true;
    }

    std::size_t missing = 0;
    for (bool value : filled) {
        if (!value) ++missing;
    }

    if (missing != 0) {
        throw std::runtime_error(
            "COMSOLMagneticField: regular grid is incomplete. Missing points: " +
            std::to_string(missing)
        );
    }

    std::cout << "\nCOMSOLMagneticField loaded successfully.\n"
              << "  file: " << filename << "\n"
              << "  nx = " << fX.size()
              << ", ny = " << fY.size()
              << ", nz = " << fZ.size() << "\n"
              << "  x range [mm]: " << fXMin / mm << " to " << fXMax / mm << "\n"
              << "  y range [mm]: " << fYMin / mm << " to " << fYMax / mm << "\n"
              << "  z range [mm]: " << fZMin / mm << " to " << fZMax / mm << "\n"
              << std::endl;
}

std::size_t COMSOLMagneticField::Index(std::size_t ix,
                                       std::size_t iy,
                                       std::size_t iz) const {
    const std::size_t ny = fY.size();
    const std::size_t nz = fZ.size();

    return (ix * ny + iy) * nz + iz;
}

std::size_t COMSOLMagneticField::NearestIndex(const std::vector<G4double>& grid,
                                              G4double value,
                                              const char* axisName) const {
    auto it = std::lower_bound(grid.begin(), grid.end(), value);

    if (it == grid.end()) {
        const std::size_t idx = grid.size() - 1;
        if (std::abs(grid[idx] - value) < kCoordinateTolerance) return idx;
    }

    if (it != grid.end() && std::abs(*it - value) < kCoordinateTolerance) {
        return static_cast<std::size_t>(std::distance(grid.begin(), it));
    }

    if (it != grid.begin()) {
        auto prev = it - 1;
        if (std::abs(*prev - value) < kCoordinateTolerance) {
            return static_cast<std::size_t>(std::distance(grid.begin(), prev));
        }
    }

    std::ostringstream msg;
    msg << "COMSOLMagneticField: coordinate value not found on "
        << axisName << " grid: " << value / mm << " mm";

    throw std::runtime_error(msg.str());
}

bool COMSOLMagneticField::FindCell(const std::vector<G4double>& grid,
                                   G4double value,
                                   std::size_t& i0,
                                   G4double& t) const {
    if (value < grid.front() || value > grid.back()) {
        return false;
    }

    if (value == grid.back()) {
        i0 = grid.size() - 2;
        t = 1.0;
        return true;
    }

    auto upper = std::upper_bound(grid.begin(), grid.end(), value);

    if (upper == grid.begin() || upper == grid.end()) {
        return false;
    }

    i0 = static_cast<std::size_t>(std::distance(grid.begin(), upper) - 1);

    const G4double x0 = grid[i0];
    const G4double x1 = grid[i0 + 1];

    t = (value - x0) / (x1 - x0);

    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    return true;
}

std::array<G4double, 3> COMSOLMagneticField::FieldAt(std::size_t ix,
                                                     std::size_t iy,
                                                     std::size_t iz) const {
    return fB[Index(ix, iy, iz)];
}

void COMSOLMagneticField::GetFieldValue(const G4double point[4],
                                        G4double* Bfield) const {
    const G4double x = point[0];
    const G4double y = point[1];
    const G4double z = point[2];

    std::size_t ix = 0;
    std::size_t iy = 0;
    std::size_t iz = 0;

    G4double tx = 0.0;
    G4double ty = 0.0;
    G4double tz = 0.0;

    const bool inside =
        FindCell(fX, x, ix, tx) &&
        FindCell(fY, y, iy, ty) &&
        FindCell(fZ, z, iz, tz);

    if (!inside) {
        Bfield[0] = 0.0;
        Bfield[1] = 0.0;
        Bfield[2] = 0.0;
        return;
    }

    const auto c000 = FieldAt(ix,     iy,     iz);
    const auto c100 = FieldAt(ix + 1, iy,     iz);
    const auto c010 = FieldAt(ix,     iy + 1, iz);
    const auto c110 = FieldAt(ix + 1, iy + 1, iz);

    const auto c001 = FieldAt(ix,     iy,     iz + 1);
    const auto c101 = FieldAt(ix + 1, iy,     iz + 1);
    const auto c011 = FieldAt(ix,     iy + 1, iz + 1);
    const auto c111 = FieldAt(ix + 1, iy + 1, iz + 1);

    for (int component = 0; component < 3; ++component) {
        const G4double c00 = c000[component] * (1.0 - tx) + c100[component] * tx;
        const G4double c10 = c010[component] * (1.0 - tx) + c110[component] * tx;
        const G4double c01 = c001[component] * (1.0 - tx) + c101[component] * tx;
        const G4double c11 = c011[component] * (1.0 - tx) + c111[component] * tx;

        const G4double c0 = c00 * (1.0 - ty) + c10 * ty;
        const G4double c1 = c01 * (1.0 - ty) + c11 * ty;

        Bfield[component] = c0 * (1.0 - tz) + c1 * tz;
    }
}