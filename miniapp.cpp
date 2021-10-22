// Charles Goodman
// 18 October 2021
// Deterministic Transport Solver
// monoenergetic, 1,2,3-D, implicit time, discrete ordinates

#include <mfem.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>

double inflow_function(const mfem::Vector &x);
double u0_function(const mfem::Vector &x);


int main(int argc, char *argv[]) {
	

        // Enable GPU (when executed with option -d cuda)
        const char *device_config = "cpu";

        mfem::OptionsParser args(argc, argv);
        args.AddOption(&device_config, "-d", "--device",
            "Device configuration string, see Device::Configure().");

        args.Parse();
        if (!args.Good())
        {
          args.PrintUsage(std::cout);
          return 1;
        }
        args.PrintOptions(std::cout);

        mfem::Device device(device_config);
        device.Print();

	// Read Input Files
	// Geometry
	mfem::Mesh mesh("1_10x1_10.mesh", 1, 1);
	int dim = mesh.Dimension();

	// Problem Data
	std::ifstream matf("B1.mat");
	
	// Boundary Conditions (unused)
	int lboundtype;
	double lboundval;
	int rboundtype;
	double rboundval;
	
	matf >> lboundtype;
	matf >> lboundval;
	matf >> rboundtype;
	matf >> rboundval;

	//Materials
	int Nmats;
	matf >> Nmats;

	double tmp;
	double tmp2;
	double tmp3;
	
	mfem::Vector SigmaT(mesh.attributes.Max());
	mfem::Vector SigmaS(mesh.attributes.Max());
	mfem::Vector Source(mesh.attributes.Max());

	for (int i = 0; i < Nmats; i++) {
		matf >> tmp;
		SigmaT(i) = tmp;
		matf >> tmp;
		matf >> tmp2;
		matf >> tmp3;
		SigmaS(i) = (tmp + tmp2*tmp3)/8; //2 to power of quad dim;
		matf >> tmp;
		Source(i) = tmp/4/3.141592;
	}
	matf.close();

	mfem::PWConstCoefficient sigmaT(SigmaT);
	mfem::PWConstCoefficient sigmaS(SigmaS);
	mfem::PWConstCoefficient source(Source);


	// Quadrature Set
	// Write Quadrature set angles from file to a vector of unit vectors
	std::ifstream quadf("LS16.quad");
	int quadDim;
	quadf >> quadDim;
	int M;
	quadf >> M;
	std::vector<double> QuadW;
	std::vector<mfem::Vector> Quad;
	for (int i = 0; i < M; i++) {
		
		double angle[quadDim];
		double tmp;
		for (int j = 0; j < quadDim; j++) {
			quadf >> angle[j];
		}
		mfem::Vector Omega;
		Omega.NewDataAndSize(angle,quadDim);
		Quad.push_back(Omega);
		quadf >> tmp;
		QuadW.push_back(tmp);
	}
	quadf.close();

	int kend = 50;
	int k = 0;
	double dt = .02;
	double v = 30; //52; // cm/shake for 14 MeV neutron
	
	//Define the discontinous DG finite element space of given polynomial order on the mesh
	int order = 1; //This is Bilinear Discontinuous for 2D problems
	mfem::DG_FECollection fec(order, dim, mfem::BasisType::GaussLobatto);
	mfem::FiniteElementSpace fes(&mesh, &fec);
	mfem::FiniteElementSpace vfes(&mesh, &fec, dim);
				

	mfem::Array<int> ess_tdof_list, ess_bdr(mesh.bdr_attributes.Max());
	ess_bdr = 0;
	ess_bdr[0] = 1;
  fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

	std::vector<mfem::BilinearForm*> A;
	for ( int m = 0; m < M; m++) { //For each discrete ordinate
		mfem::VectorConstantCoefficient Omega(Quad[m]);
		
		mfem::BilinearForm *S = new mfem::BilinearForm(&fes);
		mfem::ConstantCoefficient timeconst (1.0/v/dt);
		S->AddDomainIntegrator(new mfem::MassIntegrator(timeconst));
		S->AddDomainIntegrator(new mfem::ConvectionIntegrator(Omega, 1.0));
		S->AddDomainIntegrator(new mfem::MassIntegrator(sigmaT));
		S->AddInteriorFaceIntegrator(new mfem::TransposeIntegrator(new mfem::DGTraceIntegrator(Omega, -1.0, 0.5)));
		S->AddBdrFaceIntegrator(new mfem::TransposeIntegrator(new mfem::DGTraceIntegrator(Omega, -1.0, 0.5)));
		
		S->Assemble();
		S->Finalize();
		
		//A is a vector of Matrix such that Ax = b
		A.push_back(S);

	}
	
	//Set initial condition
	mfem::FunctionCoefficient u0(u0_function);
	std::vector<mfem::GridFunction> u;
	std::vector<mfem::GridFunction> uprev;
	for (int m = 0; m < M; m++) {
		mfem::GridFunction ut(&fes);
		ut.ProjectCoefficient(u0);
		u.push_back(ut);
		ut /= v*dt;
		uprev.push_back(ut);
	}
		
	
	mfem::GridFunction scalarflux(&fes); 
	mfem::GridFunction current(&vfes);
	std::vector<mfem::GridFunction*> currentcom (dim);
	for (int d = 0; d < dim; d++) {
		currentcom[d] = new mfem::GridFunction(&fes);
	}
	 
	scalarflux = 0.0;
		
	// Time Loop
	while (k < kend) {
		k++;

		//Initialize source from previous time step
		
		int inneriterations = 0;
		while (true) {
			inneriterations++;
			//converge source iterations
			mfem::GridFunction oldscalarflux = scalarflux;
			mfem::GridFunctionCoefficient oldscalarfluxc(&oldscalarflux);
			mfem::ProductCoefficient scatsource(sigmaS, oldscalarfluxc); 
			scalarflux = 0.0;
			current = 0.0;
			for (int d = 0; d < dim; d++) {
				*currentcom[d] = 0;
			}
				
			for (int m = 0; m < M; m++) {
				//run mfem over domain for each ordinate
				
				//make these two function coefficient functions of m and t
				mfem::VectorConstantCoefficient Omega(Quad[m]);
				mfem::FunctionCoefficient inflow(inflow_function);

				//previous time step angular flux
				mfem::GridFunctionCoefficient tprev(&uprev[m]);

				//Here is the iterated source
				mfem::LinearForm b(&fes);
				//Scattering Source
				b.AddDomainIntegrator(new mfem::DomainLFIntegrator(scatsource));

				//Time Source
				b.AddDomainIntegrator(new mfem::DomainLFIntegrator(tprev));
			
				//External source
				b.AddDomainIntegrator(new mfem::DomainLFIntegrator(source));

				//Boundary Conditions
				b.AddBdrFaceIntegrator(new mfem::BoundaryFlowIntegrator(inflow, Omega, -1.0, -0.5));
				b.Assemble();
		
				//A[m] is a pointer
				//Need to perform u = u + Ainv*b	
				mfem::BlockILU ilu(A[m]->SpMat());

				mfem::GMRESSolver solver;
				solver.iterative_mode=false;
				solver.SetRelTol(1e-14);
				solver.SetAbsTol(0.0);
				solver.SetMaxIter(100);
				solver.SetPrintLevel(0);
				solver.SetPreconditioner(ilu);
				solver.SetOperator(A[m]->SpMat());
				
				u[m] = 0;
				solver.Mult(b, u[m]);
			
				
				//need to use quadrature weights
				mfem::GridFunction tmpgf = u[m];
				tmpgf *= QuadW[m];
				//calculate scalar flux
				scalarflux += tmpgf;
				//calculate current
				//can be done componentwise and recombined in vector at end
				for (int d = 0; d < dim; d++) {
					mfem::GridFunction tmpgf2 = tmpgf;
					tmpgf2 *= Quad[m].Elem(d);
					*currentcom[d] += tmpgf2;
				}
				//tries to perform inner product
				//current += Quad[m]*u[m];
				
			}

			//compute residual in angular flux calc

			int order_to_integrate = 4;
			const mfem::IntegrationRule *irs[mfem::Geometry::NumGeom];
			for (int i=0; i < mfem::Geometry::NumGeom; ++i) {
    		irs[i] = &(mfem::IntRules.Get(i, order_to_integrate));
			}			
			
			std::cout << scalarflux.ComputeMaxError(oldscalarfluxc,irs) << std::endl;

			//end the loop when flux is converged
			if (scalarflux.ComputeMaxError(oldscalarfluxc,irs) < 1E-12) {
				std::cout << inneriterations << std::endl;
				std::cout << scalarflux.ComputeMaxError(oldscalarfluxc,irs) << std::endl;
				break;
			}

		}
	
		//save previous time step angular flux
		for (int m = 0; m < M; m++) {
			uprev[m] = u[m];
			uprev[m] /= v*dt;
		}
		
	
		//compute local residuals
		//mfem::GridFunction absorption(&fes);
		//absorption = scalarflux;
		//absorption *= (SigmaT-SigmaS);
		//mfem::DivergenceGridFunctionCoefficient divcurrent(&current);
	
                //write solution at step k to output file
                system(("mkdir -p /gpfs/wolf/gen167/scratch/reynolaa/output/"+std::to_string(k)).c_str());
                std::ofstream osol("/gpfs/wolf/gen167/scratch/reynolaa/output/"+std::to_string(k)+"/sf.gf");
		osol.precision(16);
		scalarflux.Save(osol);
		osol.close();
		for (int d = 0; d < dim; d++) {
                        std::ofstream osol("/gpfs/wolf/gen167/scratch/reynolaa/output/"+std::to_string(k)+"/curr"+std::to_string(d)+".gf");
			osol.precision(16);
			currentcom[d]->Save(osol);
			osol.close();
		}
		
	}	

}

// Inflow boundary condition
double inflow_function(const mfem::Vector &x) {
	if (x[0] == 0) {
		return 100;
	} else {
		return 0;
	}
}

// Initial Condition
double u0_function(const mfem::Vector &x) {
	return .001;
}

