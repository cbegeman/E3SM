//===-- Test driver for OMEGA time steppers -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA time steppers
///
/// This driver tests the OMEGA time stepping module. It checks that every
/// implemented time stepping scheme converges to an exact solution at
/// the expected theoretical rate. To only measure time convergence,
/// the initial conditions and tendency terms are set up such that for
/// every mesh element an independent ODE is solved.
///
//
//===-----------------------------------------------------------------------===/

#include "TimeStepper.h"
#include "../ocn/OceanTestCommon.h"
#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "CustomTendencyTerms.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

// Only one vertical layer is needed
constexpr int NVertLayers = 1;

namespace {
ManufacturedSolution gManufacturedSolution;

// Produce the manufactured solution normal velocity at a given time
Array2DReal manufacturedNormalVelocity(
    const ManufacturedSolution::ManufacturedVelocityTendency &VelTend,
    const HorzMesh *Mesh, TimeInstant Time) {
   Array2DReal AnalyticVel("ManufacturedNormalVel", Mesh->NEdgesSize,
                            Mesh->NVertLayers);

   TimeInterval ElapsedInterval = Time - VelTend.ReferenceTime;
   R8 ElapsedTimeSec;
   ElapsedInterval.get(ElapsedTimeSec, TimeUnits::Seconds);

   const auto FEdge     = Mesh->FEdge;
   const auto AngleEdge = Mesh->AngleEdge;
   const auto XEdge     = Mesh->XEdge;
   const auto YEdge     = Mesh->YEdge;
   const R8 LocGrav     = VelTend.Grav;
   const R8 LocEta0     = VelTend.Eta0;
   const R8 LocKx       = VelTend.Kx;
   const R8 LocKy       = VelTend.Ky;
   const R8 LocAngFreq  = VelTend.AngFreq;
   const int NLevels     = Mesh->NVertLayers;

   Kokkos::parallel_for(
       "ManufacturedNormalVel", {Mesh->NEdgesAll},
       KOKKOS_LAMBDA(int IEdge) {
          R8 X     = XEdge(IEdge);
          R8 Y     = YEdge(IEdge);
          R8 Phase = LocKx * X + LocKy * Y - LocAngFreq * ElapsedTimeSec;
          R8 SourceTerm0 = LocAngFreq * sin(Phase) -
                           0.5_Real * LocEta0 * (LocKx + LocKy) *
                               sin(2.0_Real * Phase);
          R8 U = LocEta0 *
                 ((-FEdge(IEdge) + LocGrav * LocKx) * cos(Phase) +
                  SourceTerm0);
          R8 V = LocEta0 *
                 ((FEdge(IEdge) + LocGrav * LocKy) * cos(Phase) +
                  SourceTerm0);
          R8 NormalComp = cos(AngleEdge(IEdge)) * U +
                          sin(AngleEdge(IEdge)) * V;
          for (int K = 0; K < NLevels; ++K) {
             AnalyticVel(IEdge, K) = NormalComp;
          }
       });

   return AnalyticVel;
}
} // namespace

int initState(TimeInstant Time,
              const ManufacturedSolution::ManufacturedVelocityTendency &VelTend,
              const HorzMesh *Mesh) {
   int Err = 0;
   auto *State             = OceanState::get("TestState");
   Array3DReal TracerArray = Tracers::getAll(0);

   Array2DReal LayerThickCell = State->getLayerThickness(0);
   Array2DReal NormalVelEdge  = State->getNormalVelocity(0);

   // Initially set thickness and velocity and tracers to 1
   deepCopy(LayerThickCell, 1);
   Array2DReal AnalyticVel = manufacturedNormalVelocity(VelTend, Mesh, Time);
   deepCopy(NormalVelEdge, AnalyticVel);
   deepCopy(TracerArray, 1);

   return Err;
}

//------------------------------------------------------------------------------
// The initialization routine for time stepper testing
int initTimeStepperTest(const std::string &mesh) {
   int Err = 0;
   Error Err1;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Default init

   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   Config *OmegaConfig = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err += OmegaConfig->get(TendConfig);
   TendConfig.set("UseCustomTendency", true);
   TendConfig.set("ManufacturedSolutionTendency", true);

   // Note that the default time stepper is not used in subsequent tests
   // but is initialized here because the number of time levels is needed
   // to initialize the Tracers. If a later timestepper test uses more time
   // levels than the default, this unit test may fail.
   TimeStepper::init1();

   IO::init(DefComm);
   Decomp::init(mesh);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default halo");
   }

   HorzMesh::init();

   // initialize vertical coordinate, do not read stream and use local
   // NVertLayers value
   VertCoord::init(false, NVertLayers);
   auto *DefVertCoord = VertCoord::getDefault();

   Tracers::init();
   AuxiliaryState::init();
   Tendencies::init();
   gManufacturedSolution.init();

   // finish initializing default time stepper
   TimeStepper::init2();

   // Default time stepper never used and time stepper tests require
   // the no calendar option, so we reset calendar here
   Calendar::reset();
   Calendar::init("No Calendar");

   // Non-default init
   // Creating non-default state and auxiliary state to use only one vertical
   // layer

   auto *DefMesh = HorzMesh::getDefault();
   auto *DefHalo = Halo::getDefault();

   int NTracers          = Tracers::getNumTracers();
   const int NTimeLevels = 2;
   auto *TestOceanState  = OceanState::create("TestState", DefMesh, DefHalo,
                                              NVertLayers, NTimeLevels);
   if (!TestOceanState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test state");
   }

   auto *TestAuxState = AuxiliaryState::create("TestAuxState", DefMesh, DefHalo,
                                               DefVertCoord, NTracers);

   Config *OmegaConfig = Config::getOmegaConfig();
   TestAuxState->readConfigOptions(OmegaConfig);

   if (!TestAuxState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test auxiliary state");
   }

   Config Options;

   // Creating non-default tendencies with custom velocity tendencies
   auto *TestTendencies = Tendencies::create(
      "TestTendencies", DefMesh, DefVertCoord, NTracers, &Options);
   if (!TestTendencies) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test tendencies");
   }

   // Disable all other tendencies
   TestTendencies->ThicknessFluxDiv.Enabled   = false;
   TestTendencies->PotientialVortHAdv.Enabled = false;
   TestTendencies->KEGrad.Enabled             = false;
   TestTendencies->SSHGrad.Enabled            = false;
   TestTendencies->VelocityDiffusion.Enabled  = false;
   TestTendencies->VelocityHyperDiff.Enabled  = false;
   TestTendencies->TracerHorzAdv.Enabled      = false;
   TestTendencies->TracerDiffusion.Enabled    = false;
   TestTendencies->TracerHyperDiff.Enabled    = false;
   TestTendencies->WindForcing.Enabled        = false;
   TestTendencies->BottomDrag.Enabled         = false;

   return Err;
}

// Slightly adjust time step so that it evenly divides TimeEnd and return number
// of steps
int adjustTimeStep(TimeStepper *Stepper, Real TimeEnd) {
   TimeInterval TimeStep = Stepper->getTimeStep();
   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);

   const int NSteps = std::ceil(TimeEnd / TimeStepSeconds);

   TimeStepSeconds = TimeEnd / NSteps;
   TimeStep.set(TimeStepSeconds, TimeUnits::Seconds);
   Stepper->changeTimeStep(TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

   return NSteps;
}

void timeLoop(TimeInstant TimeStart, Real TimeEnd) {
   auto *Stepper = TimeStepper::get("TestTimeStepper");
   auto *State   = OceanState::get("TestState");

   const int NSteps            = adjustTimeStep(Stepper, TimeEnd);
   const TimeInterval TimeStep = Stepper->getTimeStep();
   TimeInstant CurTime         = Stepper->getStartTime();

   // Time loop
   for (int Step = 0; Step < NSteps; ++Step) {
      Stepper->doStep(State, CurTime);
   }
}

void finalizeTimeStepperTest() {
   Tracers::clear();
   TimeStepper::clear();
   Tendencies::clear();
   AuxiliaryState::clear();
   OceanState::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Dimension::clear();
   Field::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int testTimeStepper(const std::string &Name, TimeStepperType Type,
                    Real ExpectedOrder, Real ATol) {
   int Err = 0;

   // Set pointers to data
   auto *DefMesh        = HorzMesh::getDefault();
   auto *DefVCoord      = VertCoord::getDefault();
   auto *DefHalo        = Halo::getDefault();
   auto *TestAuxState   = AuxiliaryState::get("TestAuxState");
   auto *TestTendencies = Tendencies::get("TestTendencies");

   // Set time information
   const TimeInstant TimeStart(0, 0, 0, 0, 0, 0);

   const Real TimeEnd = 1;
   TimeInstant TimeEndTI(0, 0, 0, 0, 0, 1);

   const Real BaseTimeStepSeconds = 0.2;
   TimeInterval TimeStepTI(BaseTimeStepSeconds, TimeUnits::Seconds);

   auto *TestTimeStepper = TimeStepper::create(
       "TestTimeStepper", Type, TimeStart, TimeEndTI, TimeStepTI,
       TestTendencies, TestAuxState, DefMesh, DefVCoord, DefHalo);

   if (!TestTimeStepper) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test time stepper {}", Name);
   }

   int NRefinements = 2;
   std::vector<ErrorMeasures> Errors(NRefinements);

   // This creates global exact solution and needs to be done only once
   R8 TimeStepSeconds = BaseTimeStepSeconds;

   // Convergence loop
   for (int RefLevel = 0; RefLevel < NRefinements; ++RefLevel) {
      TestTimeStepper->changeTimeStep(
          TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

      Err += initState(TimeStart, gManufacturedSolution.ManufacturedVelTend,
                       DefMesh);

      timeLoop(TimeStart, TimeEnd);

      auto *State = OceanState::get("TestState");
      Array2DReal NormalVelEdge = State->getNormalVelocity(0);
      Array2DReal ExactNormalVelEdge = manufacturedNormalVelocity(
          gManufacturedSolution.ManufacturedVelTend, DefMesh, TimeEndTI);
      ErrorMeasures VelErrors;
      computeErrors(VelErrors, NormalVelEdge, ExactNormalVelEdge, DefMesh,
                    OnEdge);
      Errors[RefLevel] = VelErrors;

      TimeStepSeconds /= 2;
   }

   std::vector<Real> ConvRates(NRefinements - 1);
   for (int RefLevel = 0; RefLevel < NRefinements - 1; ++RefLevel) {
      ConvRates[RefLevel] =
          std::log2(Errors[RefLevel].LInf / Errors[RefLevel + 1].LInf);
   }

   if (std::abs(ConvRates.back() - ExpectedOrder) > ATol) {
      Err++;
      LOG_ERROR(
          "Wrong convergence rate for time stepper {}, got {:.3f}, expected {}",
          Name, ConvRates.back(), ExpectedOrder);
   }

   TimeStepper::erase("TestTimeStepper");

   return Err;
}

int timeStepperTest(const std::string &MeshFile = "OmegaMesh.nc") {

   int Err = initTimeStepperTest(MeshFile);

   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Error initializing");
   }

   // Test convergence rate of different time steppers

   Real ExpectedOrder = 4;
   Real ATol          = 0.1;
   Err += testTimeStepper("RungeKutta4", TimeStepperType::RungeKutta4,
                          ExpectedOrder, ATol);

   ExpectedOrder = 1;
   ATol          = 0.1;
   Err += testTimeStepper("ForwardBackward", TimeStepperType::ForwardBackward,
                          ExpectedOrder, ATol);

   ExpectedOrder = 2;
   ATol          = 0.1;
   Err += testTimeStepper("RungeKutta2", TimeStepperType::RungeKutta2,
                          ExpectedOrder, ATol);

   if (Err == 0) {
      LOG_INFO("TimeStepperTest: Successful completion");
   }

   finalizeTimeStepperTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   LOG_INFO("----- Time Stepper Unit Test -----");

   RetVal += timeStepperTest("OmegaSphereMesh.nc");

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
