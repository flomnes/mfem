#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

void velocity_function(const Vector &x, Vector &v)
{
   int dim = x.Size();
   switch (dim)
   {
      case 1: v(0) = 1.0; break;
      case 2: v(0) = x(1); v(1) = -x(0); break;
      case 3: v(0) = x(1); v(1) = -x(0); v(2) = x(0); break;
   }
}

void AddConvectionIntegrators(BilinearForm &k, VectorCoefficient &velocity,
                              bool dg)
{
   k.AddDomainIntegrator(new ConvectionIntegrator(velocity, -1.0));

   if (dg)
   {
      k.AddInteriorFaceIntegrator(
         new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
      k.AddBdrFaceIntegrator(
         new TransposeIntegrator(new DGTraceIntegrator(velocity, 1.0, -0.5)));
   }
}

void test_assembly_level(Mesh &mesh, int order, bool dg, const int pb,
                         const AssemblyLevel assembly)
{
   mesh.EnsureNodes();
   mesh.SetCurvature(mesh.GetNodalFESpace()->GetOrder(0));
   int dim = mesh.Dimension();

   FiniteElementCollection *fec;
   if (dg)
   {
      fec = new L2_FECollection(order, dim, BasisType::GaussLobatto);
   }
   else
   {
      fec = new H1_FECollection(order, dim);
   }

   FiniteElementSpace fespace(&mesh, fec);

   BilinearForm k_test(&fespace);

   ConstantCoefficient one(1.0);
   VectorFunctionCoefficient vel_coeff(dim, velocity_function);

   if (pb==0) // Mass
   {
      k_test.AddDomainIntegrator(new MassIntegrator(one));
   }
   else if (pb==1) // Convection
   {
      AddConvectionIntegrators(k_test, vel_coeff, dg);
   }
   else if (pb==2) // Diffusion
   {
      k_test.AddDomainIntegrator(new DiffusionIntegrator(one));
   }

   tic_toc.Clear();
   tic_toc.Start();
   k_test.SetAssemblyLevel(assembly);
   k_test.Assemble();
   k_test.Finalize();
   tic_toc.Stop();

   cout << "      Assembly time: " <<  tic_toc.RealTime() << "sec." << endl;


   // GridFunction x(&fespace), y_test(&fespace);

   // x.Randomize(1);

   // tic_toc.Clear();
   // tic_toc.Start();
   // k_test.Mult(x,y_test);
   // tic_toc.Stop();
  
   // cout << "      Apply time: " <<  tic_toc.RealTime() << "sec." << endl;

   delete fec;
}

void test_assembly_levels(Mesh &&mesh, int order, bool dg, const int pb)
{
   const int dim = mesh.Dimension();
   int ref_levels = (int)floor(log(50000./mesh.GetNE())/log(2.)/dim);
   for (int l = 0; l < ref_levels; l++)
   {
      mesh.UniformRefinement();
   }
   cout << "    The number of elements is: " << mesh.GetNE() << endl;
   for (AssemblyLevel assembly : {AssemblyLevel::PARTIAL,AssemblyLevel::ELEMENT,AssemblyLevel::FULL,AssemblyLevel::LEGACYFULL})
   {
      cout << "    Assembly level: ";
      switch (assembly)
      {
      case AssemblyLevel::PARTIAL:
         cout << "PARTIAL" << endl;
         break;
      case AssemblyLevel::ELEMENT:
         cout << "ELEMENT" << endl;
         break;
      case AssemblyLevel::FULL:
         cout << "FULL" << endl;
         break;
      case AssemblyLevel::LEGACYFULL:
         cout << "LEGACYFULL" << endl;
         break;
      default:
         break;
      }
      test_assembly_level(mesh, order, dg, pb, assembly);
   }
}

int main(int argc, char *argv[])
{
   for (int pb : {0, 1, 2})
   {
      switch (pb)
      {
      case 0:
         cout << "Mass" << endl;
         break;
      case 1:
         cout << "Convection" << endl;
         break;
      case 2:
         cout << "Diffusion" << endl;
         break;      
      default:
         break;
      }
      for (bool dg : {true, false})
      {
         cout << (dg ? "DG" : "CG") << endl;
         for (int order : {1, 2, 3, 4, 5, 6, 7, 8})
         {
            cout << "  order " << order << endl;
            test_assembly_levels(Mesh("../data/periodic-square.mesh", 1, 1), order, dg, pb);
            cout << endl;
            test_assembly_levels(Mesh("../data/periodic-hexagon.mesh", 1, 1), order, dg, pb);
            cout << endl;
            test_assembly_levels(Mesh("../data/star-q3.mesh", 1, 1), order, dg, pb);
            cout << endl;
            test_assembly_levels(Mesh("../data/periodic-cube.mesh", 1, 1), order, dg, pb);
            cout << endl;
            test_assembly_levels(Mesh("../data/fichera-q3.mesh", 1, 1), order, dg, pb);
            cout << endl;
         }
      }

      cout << "  AMR" << endl;
      // Test AMR cases (DG not implemented)
      for (int order : {2, 3, 4})
      {
         cout << "  order " << order << endl;
         test_assembly_levels(Mesh("../data/amr-quad.mesh", 1, 1), order, false, 0);
         cout << endl;
      }
      int order = 2;
      cout << "  order " << order << endl;
      test_assembly_levels(Mesh("../data/fichera-amr.mesh", 1, 1), order, false, 0);
      cout << endl;
   }
}