#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <iomanip>
#include <mpi.h>

#include "DataStructs.h"
#include "rk4.h"
#include "FluxFunctions.h"
#include "RHSoperator.h"

#ifdef _DOUBLE_
#define FLOATTYPE double
#else
#define FLOATTYPE float
#endif

void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name);
FLOATTYPE calcL2norm(DataStruct<FLOATTYPE> &u, DataStruct<FLOATTYPE> &uinit);

int main(int narg, char **argv)
{
  // --- INICIO DE MPI ---
  int rank, size;
  MPI_Init(&narg, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  // ---------------------

  int numPoints = 80;
  FLOATTYPE k = 2.; 

  if(narg != 3)
  {
    if(rank == 0) {
        std::cout<< "Wrong number of arguments. You should include:" << std::endl;
        std::cout<< "    Num points" << std::endl;
        std::cout<< "    Wave number" << std::endl;
    }
    MPI_Finalize();
    return 1;
  }
  else
  {
    numPoints = std::stoi(argv[1]);
    k         = std::stod(argv[2]);
  }

  DataStruct<FLOATTYPE> u(numPoints), f(numPoints), xj(numPoints);
  LinearFlux<FLOATTYPE> lf;
  RungeKutta4<FLOATTYPE> rk(u);

  FLOATTYPE *datax = xj.getData();
  FLOATTYPE *dataU = u.getData();
  for(int j = 0; j < numPoints; j++)
  {
    datax[j] = FLOATTYPE(j)/FLOATTYPE(numPoints-1);
    dataU[j] = sin(k*2. * M_PI * datax[j]);
  }

  DataStruct<FLOATTYPE> Uinit;
  Uinit = u;

  Central1D<FLOATTYPE> rhs(u,xj,lf);

  FLOATTYPE CFL = 2.4;
  FLOATTYPE dt = CFL*datax[1];

  // Solo el rank 0 escribe la condición inicial
  if(rank == 0) {
      write2File(xj, u, "initialCondition.csv");
  }

  FLOATTYPE t_final = 1.;
  FLOATTYPE time = 0.;
  DataStruct<FLOATTYPE> Ui(u.getSize()); 

  double compTime = MPI_Wtime();

  while(time < t_final)
  {
    if(time+dt >= t_final) dt = t_final - time;

    rk.initRK();
    for(int s = 0; s < rk.getNumSteps(); s++)
    {
      rk.stepUi(dt);
      Ui = *rk.currentU();
      rhs.eval(Ui);
      rk.setFi(rhs.ref2RHS());
    }
    rk.finalizeRK(dt);
    time += dt;
  }

  compTime = MPI_Wtime() - compTime;

  // Solo el rank 0 escribe el archivo final
  if(rank == 0) {
      write2File(xj, u, "final.csv");
  }

  FLOATTYPE err = calcL2norm(Uinit, u);
  
  if(rank == 0) {
      std::cout << std::setprecision(4) << "Comp. time: " << compTime;
      std::cout << " sec. Error: " << err/k;
      std::cout << " kdx: " << k*datax[1]*2.*M_PI;
      std::cout << std::endl;
  }

  // Sincronizar y cerrar MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  return 0;
}

void write2File(DataStruct<FLOATTYPE> &X, DataStruct<FLOATTYPE> &U, std::string name)
{
  std::ofstream file;
  file.open(name, std::ios_base::trunc);
  if(!file.is_open()) 
  {
    std::cout << "Couldn't open file for " << name << std::endl;
    exit(1);
  }
  
  // Evitamos imprimir las celdas fantasma (empezamos en 1 y terminamos en size-1)
  for(int j = 1; j < U.getSize() - 1; j++)
  {
    file << X.getData()[j] << " ," << U.getData()[j] << std::endl;
  }

  file.close();
}

FLOATTYPE calcL2norm(DataStruct<FLOATTYPE> &u, DataStruct<FLOATTYPE> &uinit)
{
  FLOATTYPE err = 0.;
  const FLOATTYPE *dataU = u.getData();
  const FLOATTYPE *dataInit = uinit.getData();

  for(int n = 0; n < u.getSize(); n++)
  {
    err += (dataU[n] - dataInit[n])*(dataU[n] - dataInit[n]);
  }

  return sqrt( err );
}