#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SimConfig.hh"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cmath>

namespace {

G4double Uniform(G4double halfWidth) {
    return (2.0 * G4UniformRand() - 1.0) * halfWidth;
}

G4double SamplePowerLawEnergy(G4double emin, G4double emax, G4double gamma) {
    // Samples f(E) ~ E^(-gamma), between emin and emax.
    // This is a simple design-stage energy sampler, not the final Gaisser model.
    const G4double u = G4UniformRand();
    const G4double a = 1.0 - gamma;
    return std::pow(std::pow(emin, a) + u * (std::pow(emax, a) - std::pow(emin, a)), 1.0 / a);
}

G4String ChooseMuonName() {
    if (SimConfig::kForceChargeMode > 0) return "mu+";
    if (SimConfig::kForceChargeMode < 0) return "mu-";
    return (G4UniformRand() < SimConfig::kMuPlusFraction) ? "mu+" : "mu-";
}

}

PrimaryGeneratorAction::PrimaryGeneratorAction(EventAction* eventAction)
    : G4VUserPrimaryGeneratorAction(),
      fParticleGun(new G4ParticleGun(1)),
      fEventAction(eventAction) {}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    auto particleTable = G4ParticleTable::GetParticleTable();
    auto particle = particleTable->FindParticle(ChooseMuonName());

    fParticleGun->SetParticleDefinition(particle);

    const G4double kineticEnergy = SimConfig::kUseFixedEnergy
        ? SimConfig::kFixedKineticEnergy
        : SamplePowerLawEnergy(SimConfig::kEnergyMin,
                               SimConfig::kEnergyMax,
                               SimConfig::kEnergyPowerLawGamma);
    fParticleGun->SetParticleEnergy(kineticEnergy);

    // Biased slit/gap generator:
    // Pick a point in the virtual slit/top-detector window, then pick a point in the desired
    // off-center magnetic-gap window. The ray is forced to pass through both.
    const G4ThreeVector topPoint(
        SimConfig::kSlitCenterX + Uniform(SimConfig::kSlitHalfX),
        SimConfig::kSlitCenterY + Uniform(SimConfig::kSlitHalfY),
        SimConfig::kTopZ
    );

    const G4ThreeVector gapPoint(
        SimConfig::kGapTargetX + Uniform(SimConfig::kGapHalfX),
        SimConfig::kGapTargetY + Uniform(SimConfig::kGapHalfY),
        SimConfig::kGapZ
    );

    const G4ThreeVector direction = (gapPoint - topPoint).unit();

    // Back-project from the top point to the source plane.
    // direction.z() is negative, sourceZ > topZ, so t is negative.
    const G4double t = (SimConfig::kSourceZ - SimConfig::kTopZ) / direction.z();
    const G4ThreeVector sourcePosition = topPoint + t * direction;

    fParticleGun->SetParticleMomentumDirection(direction);
    fParticleGun->SetParticlePosition(sourcePosition);

    if (fEventAction) {
        fEventAction->SetPrimary(particle->GetPDGEncoding(),
                                 particle->GetPDGCharge(),
                                 kineticEnergy,
                                 sourcePosition,
                                 direction);
    }

    // Exactly one primary vertex per event. Do not add a second GeneratePrimaryVertex call.
    fParticleGun->GeneratePrimaryVertex(event);
}
