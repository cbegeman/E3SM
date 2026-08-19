//===-- ocn/MomBudgetDiag.cpp - Momentum budget diagnostics -----*- C++ -*-===//
//
// TEMPORARY debug utility. See MomBudgetDiag.h for details.
//
//===----------------------------------------------------------------------===//

#include "MomBudgetDiag.h"
#include "Decomp.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

#include <cmath>
#include <cstdio>

namespace OMEGA {

// Sampling controls. Every tendency evaluation is recorded for the first
// FullDetailCalls calls (the initial transient, where the flow starts from
// rest and each term can be checked in isolation); after that only every
// StrideCalls-th call is recorded so the file stays a manageable size over a
// multi-day run.
static constexpr I4 FullDetailCalls = 12;
static constexpr I4 StrideCalls     = 400;
static constexpr I4 MaxLevelsOut    = 24;

static const char *BudgetFileName = "mom_budget.csv";
static std::FILE *BudgetFile      = nullptr;

bool MomBudgetDiag::Initialized   = false;
bool MomBudgetDiag::Enabled       = false;
bool MomBudgetDiag::WriteThisCall = false;
I4 MomBudgetDiag::CallCount       = 0;
I4 MomBudgetDiag::StageCount      = 0;
I4 MomBudgetDiag::StepCount       = 0;
R8 MomBudgetDiag::CurrentTime     = 0.0;
I4 MomBudgetDiag::NLevelsOut      = 0;
std::vector<I4> MomBudgetDiag::SampleEdges;
std::vector<I4> MomBudgetDiag::SampleEdgeIDs;
std::vector<R8> MomBudgetDiag::SampleAngles;
HostArray2DReal MomBudgetDiag::Baseline;
HostArray2DReal MomBudgetDiag::ImplicitBefore;
HostArray2DReal MomBudgetDiag::Scratch;

//------------------------------------------------------------------------------
// Select one owned edge per distinct AngleEdge value and open the CSV

void MomBudgetDiag::init(const HorzMesh *Mesh, I4 NVertLayers) {

   if (Initialized)
      return;
   Initialized = true;

   // Only the master task writes the file. On the single-column tests this
   // is a single-rank run so all sample edges are local anyway.
   const MachEnv *DefEnv = MachEnv::getDefault();
   if (DefEnv->getMyTask() != DefEnv->getMasterTask())
      return;

   NLevelsOut = std::min(NVertLayers, MaxLevelsOut);

   // Pick the first owned edge for each distinct AngleEdge value. Angles are
   // compared with a tolerance since they are stored as floating point.
   const HostArray1DReal &AngleEdgeH = Mesh->AngleEdgeH;
   const HostArray1DI4 &EdgeIDH      = Decomp::getDefault()->EdgeIDH;
   const R8 AngleTol                 = 1.0e-6;

   for (I4 IEdge = 0; IEdge < Mesh->NEdgesOwned; ++IEdge) {
      const R8 Angle = AngleEdgeH(IEdge);
      bool Seen      = false;
      for (std::size_t I = 0; I < SampleAngles.size(); ++I) {
         if (std::fabs(SampleAngles[I] - Angle) < AngleTol) {
            Seen = true;
            break;
         }
      }
      if (!Seen) {
         SampleAngles.push_back(Angle);
         SampleEdges.push_back(IEdge);
         SampleEdgeIDs.push_back(EdgeIDH(IEdge));
      }
   }

   Baseline = HostArray2DReal("MomBudgetBaseline", Mesh->NEdgesSize,
                              NVertLayers);
   ImplicitBefore =
       HostArray2DReal("MomBudgetImplicitBefore", Mesh->NEdgesSize,
                       NVertLayers);
   Scratch =
       HostArray2DReal("MomBudgetScratch", Mesh->NEdgesSize, NVertLayers);

   BudgetFile = std::fopen(BudgetFileName, "w");
   if (BudgetFile == nullptr) {
      LOG_WARN("MomBudgetDiag: could not open {} for writing", BudgetFileName);
      return;
   }
   std::fprintf(BudgetFile,
                "step,stage,time,iEdge,globalEdge,angleEdge,k,term,tendency\n");
   Enabled = true;

   LOG_INFO("MomBudgetDiag: writing {} for {} sample edges, {} levels",
            BudgetFileName, SampleEdges.size(), NLevelsOut);
   for (std::size_t I = 0; I < SampleEdges.size(); ++I) {
      LOG_INFO("MomBudgetDiag:   IEdge={} (glob {}) AngleEdge={}",
               SampleEdges[I], SampleEdgeIDs[I], SampleAngles[I]);
   }
}

//------------------------------------------------------------------------------

void MomBudgetDiag::finalize() {
   if (BudgetFile != nullptr) {
      std::fclose(BudgetFile);
      BudgetFile = nullptr;
   }
   Enabled = false;
}

//------------------------------------------------------------------------------

bool MomBudgetDiag::active() { return Enabled && WriteThisCall; }

//------------------------------------------------------------------------------
// Start a new tendency evaluation: decide whether to record it, advance the
// counters and zero the differencing baseline.

void MomBudgetDiag::beginTendencies(R8 SimTime) {

   if (!Enabled)
      return;

   // A new time step is detected by a change in simulation time. Within a
   // step the multi-stage time integrators call the tendency routines
   // several times, tracked by StageCount.
   if (CallCount == 0 || SimTime != CurrentTime) {
      if (CallCount > 0)
         ++StepCount;
      CurrentTime = SimTime;
      StageCount  = 0;
   } else {
      ++StageCount;
   }

   WriteThisCall =
       (CallCount < FullDetailCalls) || (CallCount % StrideCalls == 0);
   ++CallCount;

   if (WriteThisCall)
      deepCopy(Baseline, 0._Real);
}

//------------------------------------------------------------------------------
// Write one term's rows for the sample edges

void MomBudgetDiag::writeRows(const std::string &TermName,
                              const HostArray2DReal &Values) {

   for (std::size_t I = 0; I < SampleEdges.size(); ++I) {
      const I4 IEdge = SampleEdges[I];
      for (I4 K = 0; K < NLevelsOut; ++K) {
         std::fprintf(BudgetFile, "%d,%d,%.6f,%d,%d,%.6f,%d,%s,%.17g\n",
                      StepCount, StageCount, CurrentTime, IEdge,
                      SampleEdgeIDs[I], SampleAngles[I], K, TermName.c_str(),
                      static_cast<double>(Values(IEdge, K)));
      }
   }
}

//------------------------------------------------------------------------------
// Explicit term: difference the accumulated tendency against the baseline

void MomBudgetDiag::recordTerm(const std::string &TermName,
                               const Array2DReal &Tend) {

   if (!active())
      return;

   deepCopy(Scratch, Tend);

   // Contribution of this term is the change since the previous term
   const I4 NRows = Scratch.extent(0);
   const I4 NCols = Scratch.extent(1);
   HostArray2DReal Delta("MomBudgetDelta", NRows, NCols);
   for (I4 I = 0; I < NRows; ++I) {
      for (I4 K = 0; K < NCols; ++K) {
         Delta(I, K) = Scratch(I, K) - Baseline(I, K);
      }
   }

   writeRows(TermName, Delta);

   // Advance the baseline so the next term reports only its own contribution
   deepCopy(Baseline, Scratch);
}

//------------------------------------------------------------------------------

void MomBudgetDiag::recordTotal(const std::string &TermName,
                                const Array2DReal &Tend) {

   if (!active())
      return;

   deepCopy(Scratch, Tend);
   writeRows(TermName, Scratch);
}

//------------------------------------------------------------------------------

void MomBudgetDiag::saveImplicitBefore(const Array2DReal &Field) {

   if (!active())
      return;

   deepCopy(ImplicitBefore, Field);
}

//------------------------------------------------------------------------------
// Implicit term: recover an effective tendency from the change in the state

void MomBudgetDiag::recordImplicit(const std::string &TermName,
                                   const Array2DReal &Field, R8 DT) {

   if (!active())
      return;
   if (DT == 0.0)
      return;

   deepCopy(Scratch, Field);

   const I4 NRows = Scratch.extent(0);
   const I4 NCols = Scratch.extent(1);
   HostArray2DReal Delta("MomBudgetImplicitDelta", NRows, NCols);
   for (I4 I = 0; I < NRows; ++I) {
      for (I4 K = 0; K < NCols; ++K) {
         Delta(I, K) = (Scratch(I, K) - ImplicitBefore(I, K)) / DT;
      }
   }

   writeRows(TermName, Delta);
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
