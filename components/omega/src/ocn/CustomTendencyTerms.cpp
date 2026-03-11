//===-- ocn/CustomTendencyTerms.cpp - Custom tendency terms -----*- C++ -*-===//
//
// The customized tendency terms can be added to the tendency terms based
// based on an option 'UseCustomTendency' in Tendencies Config group.
// This file contains functions for initializing customized tendency terms.
//
//===----------------------------------------------------------------------===//

#include "CustomTendencyTerms.h"
#include "Config.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "TimeStepper.h"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace OMEGA {

//===-----------------------------------------------------------------------===/
// Initialize the manufactured solution tendency terms.
//===-----------------------------------------------------------------------===/
void ManufacturedSolution::init() {
   Error Err; // error code

   // Get ManufacturedSolConfig group
   Config *OmegaConfig = Config::getOmegaConfig();
   Config ManufacturedSolConfig("ManufacturedSolution");
   Err += OmegaConfig->get(ManufacturedSolConfig);
   CHECK_ERROR_ABORT(
       Err,
       "ManufacturedSolution: ManufacturedSolution group not found in Config");

   // Get TendConfig group
   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   CHECK_ERROR_ABORT(
       Err, "ManufacturedSolution: Tendencies group not found in Config");

   // Get manufactured solution parameters from Config
   R8 WavelengthX;
   R8 WavelengthY;
   R8 Amplitude;

   Err += ManufacturedSolConfig.get("WavelengthX", WavelengthX);
   CHECK_ERROR_ABORT(
       Err,
       "ManufacturedSolution: WavelengthX not found in ManufacturedSolConfig");

   Err += ManufacturedSolConfig.get("WavelengthY", WavelengthY);
   CHECK_ERROR_ABORT(
       Err,
       "ManufacturedSolution: WavelengthY not found in ManufacturedSolConfig");

   Err += ManufacturedSolConfig.get("Amplitude", Amplitude);
   CHECK_ERROR_ABORT(
       Err,
       "ManufacturedSolution: Amplitude not found in ManufacturedSolConfig");

   // Get Tendendices parameters for del2 and del4 source terms
   Err += TendConfig.get("VelDiffTendencyEnable",
                         ManufacturedVelTend.VelDiffTendencyEnable);
   Err += TendConfig.get("VelHyperDiffTendencyEnable",
                         ManufacturedVelTend.VelHyperDiffTendencyEnable);
   Err += TendConfig.get("ViscDel2", ManufacturedVelTend.ViscDel2);
   Err += TendConfig.get("ViscDel4", ManufacturedVelTend.ViscDel4);

   CHECK_ERROR_ABORT(
       Err,
       "ManufacturedSolution: Could not find del2, del4 parameters in Config");

   // Get the reference time to compute the model elapsed time
   /// Get model clock from time stepper
   TimeStepper *DefStepper             = TimeStepper::getDefault();
   Clock *ModelClock                   = DefStepper->getClock();
   ManufacturedThickTend.ReferenceTime = ModelClock->getCurrentTime();
   ManufacturedVelTend.ReferenceTime   = ManufacturedThickTend.ReferenceTime;

   // Get BottomDepth for the resting thickness
   /// This test case assumes that the restingThickness is horizontally uniform
   /// and that only one vertical level is used so only one set of indices is
   /// used here.
   HorzMesh *DefHorzMesh = HorzMesh::getDefault();
   R8 H0                 = DefHorzMesh->BottomDepthH(0);

   // Define and compute common constants
   R8 Kx      = TwoPi / WavelengthX;                      // Wave in X-dir
   R8 Ky      = TwoPi / WavelengthY;                      // Wave in Y-dir
   R8 AngFreq = sqrt(H0 * Gravity * (Kx * Kx + Ky * Ky)); // Angular frequency

   // Assign constants for thickness tendency function
   ManufacturedThickTend.H0      = H0;
   ManufacturedThickTend.Eta0    = Amplitude;
   ManufacturedThickTend.Kx      = Kx;
   ManufacturedThickTend.Ky      = Ky;
   ManufacturedThickTend.AngFreq = AngFreq;

   // Assign constants for velocity tendency function
   ManufacturedVelTend.Grav    = Gravity;
   ManufacturedVelTend.Eta0    = Amplitude;
   ManufacturedVelTend.Kx      = Kx;
   ManufacturedVelTend.Ky      = Ky;
   ManufacturedVelTend.AngFreq = AngFreq;

} // end ManufacturedSolution init

//===--------------------------------------------------------------------===/
// Manufactured tendency term for the thickness equation
//===--------------------------------------------------------------------===/
void ManufacturedSolution::ManufacturedThicknessTendency::operator()(
    Array2DReal ThicknessTend, const OceanState *State,
    const AuxiliaryState *AuxState, int ThickTimeLevel, int VelTimeLevel,
    TimeInstant Time) const {

   // Get elapsed time since reference time
   R8 ElapsedTimeSec;
   TimeInterval ElapsedTimeInterval = Time - ReferenceTime;
   ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

   auto *Mesh       = HorzMesh::getDefault();
   auto NVertLayers = ThicknessTend.extent_int(1);

   Array1DReal XCell = Mesh->XCell;
   Array1DReal YCell = Mesh->YCell;

   OMEGA_SCOPE(LocH0, H0);
   OMEGA_SCOPE(LocEta0, Eta0);
   OMEGA_SCOPE(LocKx, Kx);
   OMEGA_SCOPE(LocKy, Ky);
   OMEGA_SCOPE(LocAngFreq, AngFreq);

   parallelFor(
       {Mesh->NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int ICell, int KLevel) {
          R8 X     = XCell(ICell);
          R8 Y     = YCell(ICell);
          R8 Phase = LocKx * X + LocKy * Y - LocAngFreq * ElapsedTimeSec;
          ThicknessTend(ICell, KLevel) +=
              LocEta0 *
              (-LocH0 * (LocKx + LocKy) * sin(Phase) - LocAngFreq * cos(Phase) +
               LocEta0 * (LocKx + LocKy) * cos(2.0_Real * Phase));
       });

} // end void ManufacturedThicknessTendency

//===--------------------------------------------------------------------===/
// Manufactured tendency term for the momentum equation
//===--------------------------------------------------------------------===/
void ManufacturedSolution::ManufacturedVelocityTendency::operator()(
    Array2DReal NormalVelTend, const OceanState *State,
    const AuxiliaryState *AuxState, int ThickTimeLevel, int VelTimeLevel,
    TimeInstant Time) const {

   // Get elapsed time since reference time
   R8 ElapsedTimeSec;
   TimeInterval ElapsedTimeInterval = Time - ReferenceTime;
   ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

   auto *Mesh       = HorzMesh::getDefault();
   auto NVertLayers = NormalVelTend.extent_int(1);

   Array1DReal FEdge     = Mesh->FEdge;
   Array1DReal XEdge     = Mesh->XEdge;
   Array1DReal YEdge     = Mesh->YEdge;
   Array1DReal AngleEdge = Mesh->AngleEdge;

   OMEGA_SCOPE(LocGrav, Gravity);
   OMEGA_SCOPE(LocEta0, Eta0);
   OMEGA_SCOPE(LocKx, Kx);
   OMEGA_SCOPE(LocKy, Ky);
   OMEGA_SCOPE(LocAngFreq, AngFreq);
   OMEGA_SCOPE(LocViscDel2, ViscDel2);
   OMEGA_SCOPE(LocViscDel4, ViscDel4);
   OMEGA_SCOPE(LocVelDiffTendencyEnable, VelDiffTendencyEnable);
   OMEGA_SCOPE(LocVelHyperDiffTendencyEnable, VelHyperDiffTendencyEnable);

   R8 LocKx2 = LocKx * LocKx;
   R8 LocKy2 = LocKy * LocKy;
   R8 LocKx4 = LocKx2 * LocKx2;
   R8 LocKy4 = LocKy2 * LocKy2;

   parallelFor(
       {Mesh->NEdgesAll, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int KLevel) {
          R8 X = XEdge(IEdge);
          R8 Y = YEdge(IEdge);

          R8 Phase       = LocKx * X + LocKy * Y - LocAngFreq * ElapsedTimeSec;
          R8 SourceTerm0 = LocAngFreq * sin(Phase) - 0.5_Real * LocEta0 *
                                                         (LocKx + LocKy) *
                                                         sin(2.0_Real * Phase);

          R8 U = LocEta0 *
                 ((-FEdge(IEdge) + LocGrav * LocKx) * cos(Phase) + SourceTerm0);
          R8 V = LocEta0 *
                 ((FEdge(IEdge) + LocGrav * LocKy) * cos(Phase) + SourceTerm0);

          // Del2 and del4 source terms
          if (LocVelDiffTendencyEnable) {
             U += LocViscDel2 * LocEta0 * (LocKx2 + LocKy2) * cos(Phase);
             V += LocViscDel2 * LocEta0 * (LocKx2 + LocKy2) * cos(Phase);
          }
          if (LocVelHyperDiffTendencyEnable) {
             U -= LocViscDel4 * LocEta0 *
                  ((LocKx4 + LocKy4 + LocKx2 * LocKy2) * cos(Phase));
             V -= LocViscDel4 * LocEta0 *
                  ((LocKx4 + LocKy4 + LocKx2 * LocKy2) * cos(Phase));
          }

          R8 NormalCompSourceTerm =
              cos(AngleEdge(IEdge)) * U + sin(AngleEdge(IEdge)) * V;
          NormalVelTend(IEdge, KLevel) += NormalCompSourceTerm;
       });

} // end void ManufacturedVelocityTendency

namespace {
constexpr Real kTwelveDays = 12.0 * CDay;
constexpr Real kOmegaScaler = 0.2 * 0.2 + 0.7 * 0.7 + 1.0 * 1.0;
const Real kOmegaMag = sqrt(kOmegaScaler);
constexpr Real kOmegaBase[3] = {0.2, 0.7, 1.0};

KOKKOS_INLINE_FUNCTION Real vectorNorm(const Real v[3]) {
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

KOKKOS_INLINE_FUNCTION void crossProduct(const Real a[3], const Real b[3],
                                                      Real c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

KOKKOS_INLINE_FUNCTION void lonlat2xyz(Real lon, Real lat, Real &x, Real &y,
                                                     Real &z) {
    Real cosLat = cos(lat);
    x = cos(lon) * cosLat;
    y = sin(lon) * cosLat;
    z = sin(lat);
}

KOKKOS_INLINE_FUNCTION void calcLocalEastNorth(const Real r[3], Real east[3],
                                                              Real north[3]) {
    constexpr Real kAxis[3] = {0.0, 0.0, 1.0};
    Real temp[3];
    crossProduct(kAxis, r, temp);
    Real normEast = vectorNorm(temp);
    if (normEast > 0) {
        east[0] = temp[0] / normEast;
        east[1] = temp[1] / normEast;
        east[2] = temp[2] / normEast;
    } else {
        east[0] = 1.0;
        east[1] = 0.0;
        east[2] = 0.0;
    }
    crossProduct(r, east, north);
    Real normNorth = vectorNorm(north);
    if (normNorth > 0) {
        north[0] /= normNorth;
        north[1] /= normNorth;
        north[2] /= normNorth;
    }
}

KOKKOS_INLINE_FUNCTION void flowRotation(Real lon, Real lat, Real &u, Real &v) {
    Real r[3];
    lonlat2xyz(lon, lat, r[0], r[1], r[2]);
    Real omegaScaled[3];
    Real scaleFactor = (TwoPi / kTwelveDays) / kOmegaMag;
    omegaScaled[0] = kOmegaBase[0] * scaleFactor;
    omegaScaled[1] = kOmegaBase[1] * scaleFactor;
    omegaScaled[2] = kOmegaBase[2] * scaleFactor;
    Real vel[3];
    crossProduct(omegaScaled, r, vel);
    Real east[3];
    Real north[3];
    calcLocalEastNorth(r, east, north);
    u = vel[0] * east[0] + vel[1] * east[1] + vel[2] * east[2];
    v = vel[0] * north[0] + vel[1] * north[1] + vel[2] * north[2];
}

KOKKOS_INLINE_FUNCTION void flowNondivergent(Real timeSec, Real lon, Real lat,
                                                             Real &u, Real &v) {
    Real lonP = lon - TwoPi * timeSec / kTwelveDays;
    Real cosLat = cos(lat);
    Real cost = cos(Pi * timeSec / kTwelveDays);
    Real sinLonP = sin(lonP);
    Real sinLonPSq = sinLonP * sinLonP;
    Real sin2LonP = sin(2.0_Real * lonP);
    Real sin2Lat = sin(2.0_Real * lat);
    Real factor = 1.0 / kTwelveDays;
    u = factor * (10.0_Real * sinLonPSq * sin2Lat * cost + TwoPi * cosLat);
    v = (10.0_Real * factor) * sin2LonP * cosLat * cost;
}

KOKKOS_INLINE_FUNCTION void flowDivergent(Real timeSec, Real lon, Real lat,
                                                         Real &u, Real &v) {
    Real lonP = lon - TwoPi * timeSec / kTwelveDays;
    Real cosLat = cos(lat);
    Real cost = cos(Pi * timeSec / kTwelveDays);
    Real sinLonP = sin(lonP);
    Real sinLonPHalf = sin(lonP * 0.5_Real);
    Real sinLonPHalfSq = sinLonPHalf * sinLonPHalf;
    Real sin2Lat = sin(2.0_Real * lat);
    Real cosLatSq = cosLat * cosLat;
    Real cosLatCubed = cosLatSq * cosLat;
    Real factor = 1.0 / kTwelveDays;
    u = factor * (-5.0_Real * sinLonPHalfSq * sin2Lat * cosLatSq * cost +
                      TwoPi * cosLat);
    v = (2.5_Real * factor) * sinLonP * cosLatCubed * cost;
}

KOKKOS_INLINE_FUNCTION void evaluateTransportFlow(int flowID, Real timeSec,
                                                                    Real lon, Real lat, Real &u,
                                                                    Real &v) {
    switch (flowID) {
    case 1:
        flowRotation(lon, lat, u, v);
        break;
    case 2:
    case 4:
        flowNondivergent(timeSec, lon, lat, u, v);
        break;
    case 3:
        flowDivergent(timeSec, lon, lat, u, v);
        break;
    default:
        u = 0.0;
        v = 0.0;
        break;
    }
}

KOKKOS_INLINE_FUNCTION Real computeTransportNormal(int flowID, Real timeSec,
                                                                    Real lon, Real lat,
                                                                    Real angle,
                                                                    Real sphereRadius) {
    Real u = 0.0;
    Real v = 0.0;
    evaluateTransportFlow(flowID, timeSec, lon, lat, u, v);
    return sphereRadius * (cos(angle) * u + sin(angle) * v);
}
} // namespace

TransportTestVelocityTendency::TransportTestVelocityTendency()
     : FlowID(0), ReferenceTime() {}

bool TransportTestVelocityTendency::init(Config *OmegaConfig) {
    Config TransportConfig("TransportTests");
    Error Err = OmegaConfig->get(TransportConfig);
    if (Err.isFail()) {
        FlowID = 0;
        return false;
    }

    int FlowIdValue = 0;
    Err = TransportConfig.get("FlowID", FlowIdValue);
    if (Err.isFail() || FlowIdValue <= 0) {
        FlowID = 0;
        return false;
    }

    FlowID = FlowIdValue;
    TimeStepper *DefStepper = TimeStepper::getDefault();
    Clock *ModelClock = DefStepper->getClock();
    ReferenceTime = ModelClock->getCurrentTime();
    return true;
}

bool TransportTestVelocityTendency::isEnabled() const { return FlowID > 0; }

void TransportTestVelocityTendency::operator()(
     Array2DReal NormalVelTend, const OceanState *State,
     const AuxiliaryState *AuxState, int ThickTimeLevel, int VelTimeLevel,
     TimeInstant Time) const {

    if (FlowID <= 0) {
        return;
    }

    R8 ElapsedTimeSec;
    TimeInterval ElapsedTimeInterval = Time - ReferenceTime;
    ElapsedTimeInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

    auto *Mesh       = HorzMesh::getDefault();
    auto NVertLayers = NormalVelTend.extent_int(1);

    Array1DReal LonEdge = Mesh->LonEdge;
    Array1DReal LatEdge = Mesh->LatEdge;
    Array1DReal AngleEdge = Mesh->AngleEdge;
    const int LocFlowID = FlowID;
    const R8 LocTimeSec = ElapsedTimeSec;
    const R8 SphereRadius = REarth;

    State->copyToHost(VelTimeLevel);
    HostArray2DReal NormalVelHost = State->getNormalVelocityH(VelTimeLevel);
    HostArray1DReal LonEdgeH = Mesh->LonEdgeH;
    HostArray1DReal LatEdgeH = Mesh->LatEdgeH;
    HostArray1DReal AngleEdgeH = Mesh->AngleEdgeH;
    auto NormalVelTendH = createHostMirrorCopy(NormalVelTend);
    deepCopy(NormalVelTendH, NormalVelTend);

    for (int IEdge = 0; IEdge < Mesh->NEdgesAll; ++IEdge) {
        R8 Lon = LonEdgeH(IEdge);
        R8 Lat = LatEdgeH(IEdge);
        R8 Angle = AngleEdgeH(IEdge);
        for (int KLevel = 0; KLevel < NVertLayers; ++KLevel) {
            R8 Source = computeTransportNormal(LocFlowID, LocTimeSec, Lon,
                                              Lat, Angle, SphereRadius);
            NormalVelTendH(IEdge, KLevel) +=
                Source - NormalVelHost(IEdge, KLevel);
        }
    }

    deepCopy(NormalVelTend, NormalVelTendH);
}

} // end namespace OMEGA

//=-------------------------------------------------------------------------===/
