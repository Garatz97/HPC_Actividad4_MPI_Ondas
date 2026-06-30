#include "RHSoperator.h"
#include <mpi.h> // Imprescindible para que los procesos se comuniquen

template<class T>
RHSOperator<T>::RHSOperator()
{

}

template<class T>
RHSOperator<T>::~RHSOperator()
{

}

template<class T>
Central1D<T>::Central1D(DataStruct<T> &_U, 
                     DataStruct<T> &_mesh, 
                     FluxFunction<T> &_F):
U(_U), mesh(_mesh), F(_F)
{
  RHS.setSize(_U.getSize());
}

template<class T>
Central1D<T>::~Central1D()
{

}

template<class T>
void Central1D<T>::evalRHS(DataStruct<T> &Uin)
{
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  T *dataRHS = RHS.getData();
  T *dataU = Uin.getData(); // Sin 'const' para poder actualizar las celdas fantasma
  const T *dataMesh = mesh.getData();
  const int len = U.getSize();
  
  // Adaptamos el tipo de dato de MPI dinámicamente según el template (float o double)
  MPI_Datatype mpi_type = (sizeof(T) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT;

  // 1. COMUNICACIÓN MPI, Celdas Fantasma
  int left  = (rank - 1 + size) % size;
  int right = (rank + 1) % size;

  // Enviamos nuestro primer punto físico (1) al vecino izquierdo 
  // y recibimos su último punto para guardarlo en nuestra celda derecha (len-1)
  MPI_Sendrecv(&dataU[1], 1, mpi_type, left, 0, 
               &dataU[len-1], 1, mpi_type, right, 0, 
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Enviamos nuestro último punto físico (len-2) al vecino derecho 
  // y recibimos su primer punto para guardarlo en nuestra celda izquierda (0)
  MPI_Sendrecv(&dataU[len-2], 1, mpi_type, right, 1, 
               &dataU[0], 1, mpi_type, left, 1, 
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);


  // 2. CÁLCULO NUMÉRICO
  // Iteramos EXCLUSIVAMENTE sobre los puntos físicos interiores (de 1 a len-2)
  for(int j = 1; j < len - 1; j++)
  {
    T dx = dataMesh[j+1] - dataMesh[j-1];
    dataRHS[j] = -(F.computeFlux(dataU[j+1]) - F.computeFlux(dataU[j-1]))/dx;
  }

  // Las fronteras locales (fantasmas) toman el valor adyacente para mantener consistencia
  dataRHS[0] = dataRHS[1];
  dataRHS[len-1] = dataRHS[len-2];
}

template<class T>
void Central1D<T>::eval()
{
  evalRHS(U);
}

template<class T>
void Central1D<T>::eval(DataStruct<T> &Uin)
{
  evalRHS(Uin);
}

template<class T>
DataStruct<T>& Central1D<T>::ref2RHS()
{
  return RHS;
}

template class RHSOperator<float>;
template class RHSOperator<double>;

template class Central1D<float>;
template class Central1D<double>;