// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_TMOP_TOOLS_HPP
#define MFEM_TMOP_TOOLS_HPP

#include "bilinearform.hpp"
#include "pbilinearform.hpp"
#include "tmop.hpp"
#include "gslib.hpp"

namespace mfem
{

// Performs the full remap advection loop.
class AdvectorCG : public AdaptivityEvaluator
{
private:
   RK4Solver ode_solver;
   Vector nodes0;
   Vector field0;
   const double dt_scale;
   const AssemblyLevel al;

   void ComputeAtNewPositionScalar(const Vector &new_nodes, Vector &new_field);
public:
   AdvectorCG(AssemblyLevel al = AssemblyLevel::FULL,
              double timestep_scale = 0.5)
      : AdaptivityEvaluator(),
        ode_solver(), nodes0(), field0(), dt_scale(timestep_scale), al(al) { }

   virtual void SetInitialField(const Vector &init_nodes,
                                const Vector &init_field);

   virtual void ComputeAtNewPosition(const Vector &new_nodes,
                                     Vector &new_field);
};

#ifdef MFEM_USE_GSLIB
class InterpolatorFP : public AdaptivityEvaluator
{
private:
   Vector nodes0;
   GridFunction field0_gf;
   FindPointsGSLIB *finder;
   Array<uint> el_id_out, code_out, task_id_out;
   Vector pos_r_out, dist_p_out;
   int dim;
public:
   InterpolatorFP() : finder(NULL) { }

   virtual void SetInitialField(const Vector &init_nodes,
                                const Vector &init_field);

   virtual void ComputeAtNewPosition(const Vector &new_nodes,
                                     Vector &new_field);

   ~InterpolatorFP()
   {
      finder->FreeData();
      delete finder;
   }
};
#endif

/// Performs a single remap advection step in serial.
class SerialAdvectorCGOper : public TimeDependentOperator
{
protected:
   const Vector &x0;
   Vector &x_now;
   GridFunction &u;
   VectorGridFunctionCoefficient u_coeff;
   mutable BilinearForm M, K;
   const AssemblyLevel al;

public:
   /** Here @a fes is the FESpace of the function that will be moved. Note
       that Mult() moves the nodes of the mesh corresponding to @a fes. */
   SerialAdvectorCGOper(const Vector &x_start, GridFunction &vel,
                        FiniteElementSpace &fes,
                        AssemblyLevel al = AssemblyLevel::FULL);

   virtual void Mult(const Vector &ind, Vector &di_dt) const;
};

#ifdef MFEM_USE_MPI
/// Performs a single remap advection step in parallel.
class ParAdvectorCGOper : public TimeDependentOperator
{
protected:
   const Vector &x0;
   Vector &x_now;
   GridFunction &u;
   VectorGridFunctionCoefficient u_coeff;
   mutable ParBilinearForm M, K;
   const AssemblyLevel al;

public:
   /** Here @a pfes is the ParFESpace of the function that will be moved. Note
       that Mult() moves the nodes of the mesh corresponding to @a pfes. */
   ParAdvectorCGOper(const Vector &x_start, GridFunction &vel,
                     ParFiniteElementSpace &pfes,
                     AssemblyLevel al = AssemblyLevel::FULL);

   virtual void Mult(const Vector &ind, Vector &di_dt) const;
};
#endif

class TMOPNewtonSolver : public LBFGSSolver
{
protected:
   // 0 - Newton, 1 - LBFGS.
   int solver_type;
   bool parallel;

   // Quadrature points that are checked for negative Jacobians etc.
   const IntegrationRule &ir;

   void UpdateDiscreteTC(const TMOP_Integrator &ti, const Vector &x_new) const;

public:
#ifdef MFEM_USE_MPI
   TMOPNewtonSolver(MPI_Comm comm, const IntegrationRule &irule, int type = 0)
      : LBFGSSolver(comm), solver_type(type), parallel(true), ir(irule) { }
#endif
   TMOPNewtonSolver(const IntegrationRule &irule, int type = 0)
      : LBFGSSolver(), solver_type(type), parallel(false), ir(irule) { }

   virtual double ComputeScalingFactor(const Vector &x, const Vector &b) const;

   virtual void ProcessNewState(const Vector &x) const;

   virtual void Mult(const Vector &b, Vector &x) const
   {
      if (solver_type == 0)
      {
         NewtonSolver::Mult(b, x);
      }
      else if (solver_type == 1)
      {
         LBFGSSolver::Mult(b, x);
      }
      else { MFEM_ABORT("Invalid type"); }
   }

   virtual void SetSolver(Solver &solver)
   {
      if (solver_type == 0)
      {
         NewtonSolver::SetSolver(solver);
      }
      else if (solver_type == 1)
      {
         LBFGSSolver::SetSolver(solver);
      }
      else { MFEM_ABORT("Invalid type"); }
   }
   virtual void SetPreconditioner(Solver &pr) { SetSolver(pr); }
   int CheckDetJpr_2D(const FiniteElementSpace*, const Vector&) const;
   int CheckDetJpr_3D(const FiniteElementSpace*, const Vector&) const;
   double MinDetJpr_2D(const FiniteElementSpace*, const Vector&) const;
   double MinDetJpr_3D(const FiniteElementSpace*, const Vector&) const;
};

void vis_tmop_metric_s(int order, TMOP_QualityMetric &qm,
                       const TargetConstructor &tc, Mesh &pmesh,
                       char *title, int position);
#ifdef MFEM_USE_MPI
void vis_tmop_metric_p(int order, TMOP_QualityMetric &qm,
                       const TargetConstructor &tc, ParMesh &pmesh,
                       char *title, int position);
#endif

}

#endif
