#ifndef OMEGA_MOMBUDGETDIAG_H
#define OMEGA_MOMBUDGETDIAG_H
//===-- ocn/MomBudgetDiag.h - Momentum budget diagnostics -------*- C++ -*-===//
//
/// \file
/// \brief TEMPORARY debug utility: per-term momentum tendency budget
///
/// Writes a CSV file containing the contribution of each momentum tendency
/// term at a small set of sample edges. The explicit terms are recovered by
/// differencing the accumulated NormalVelocityTend array between successive
/// term evaluations. The implicit vertical mixing term does not pass through
/// NormalVelocityTend at all (it overwrites NormalVelocity in place via a
/// tridiagonal solve) so it is recovered separately as
/// (u_after - u_before)/DT.
///
/// Sample edges are chosen as one owned edge per distinct AngleEdge value so
/// that every edge orientation on the mesh is represented exactly once.
///
/// THIS FILE IS DEBUG SCAFFOLDING AND IS INTENDED TO BE DELETED.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"

#include <string>
#include <vector>

namespace OMEGA {

class HorzMesh;

/// Collects and writes a per-term momentum tendency budget
class MomBudgetDiag {
 public:
   /// Select sample edges and open the CSV file. Safe to call more than once;
   /// only the first call has an effect.
   static void init(const HorzMesh *Mesh, I4 NVertLayers);

   /// Close the CSV file
   static void finalize();

   /// True when a budget record should be written for the current call. All
   /// of the record* routines below check this internally, this is exposed
   /// so callers can skip the device-to-host copies entirely.
   static bool active();

   /// Begin a new tendency evaluation. Resets the running baseline used to
   /// difference the explicit terms and advances the call/stage counters.
   static void beginTendencies(R8 SimTime);

   /// Record the contribution of one explicit term. Tend is the accumulated
   /// tendency array *after* the term has been applied; the contribution is
   /// Tend minus the previously recorded baseline. The baseline is then
   /// updated so the next call reports only the next term.
   static void recordTerm(const std::string &TermName, const Array2DReal &Tend);

   /// Record the total accumulated explicit tendency. Does not modify the
   /// baseline.
   static void recordTotal(const std::string &TermName,
                           const Array2DReal &Tend);

   /// Stash a copy of Field for a later call to recordImplicit
   static void saveImplicitBefore(const Array2DReal &Field);

   /// Record (Field - saved)/DT as the effective tendency of an implicit term
   static void recordImplicit(const std::string &TermName,
                              const Array2DReal &Field, R8 DT);

 private:
   static bool Initialized;      ///< init has run
   static bool Enabled;          ///< file opened successfully and rank is 0
   static bool WriteThisCall;    ///< current call is being recorded
   static I4 CallCount;          ///< number of tendency evaluations so far
   static I4 StageCount;         ///< RK stage within the current step
   static I4 StepCount;          ///< number of time steps so far
   static R8 CurrentTime;        ///< simulation time of the current call
   static I4 NLevelsOut;         ///< number of vertical levels written
   static std::vector<I4> SampleEdges;   ///< local indices of sample edges
   static std::vector<I4> SampleEdgeIDs; ///< global IDs of sample edges
   static std::vector<R8> SampleAngles;  ///< AngleEdge of sample edges
   static HostArray2DReal Baseline;      ///< running explicit-tendency
                                         ///< baseline
   static HostArray2DReal ImplicitBefore; ///< pre-solve velocity snapshot
   static HostArray2DReal Scratch;        ///< host mirror scratch space

   /// Write one term's rows for all sample edges and levels
   static void writeRows(const std::string &TermName,
                         const HostArray2DReal &Values);
};

} // namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_MOMBUDGETDIAG_H
