#include "mpi.h"
#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Dense>
void IdentityMatrix(int reps, int pid, int np);

int main(int argc, char **argv)
{
  /* MPI Variables */
  int np, pid;

  /* MPI setup */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  const int N = std::atoi(argv[1]);  
  IdentityMatrix(N, pid, np);

  /* finish */
  MPI_Finalize();

  return 0;
}

void IdentityMatrix(int N, int pid, int np)
{
  int Nlocal=N/np;
  Eigen::MatrixXd mat(Nlocal,N);
;
  mat.setZero();
  for(int i=0;i<Nlocal;i++){   
      mat(i,(Nlocal*pid)+i)=1;
      
  }

  
  MPI_Status status;
  int tag = 0;
  if(0==pid){//master
    std::cout<<mat<<"\n";
                                //printmat(mat,N,Nlocal);
    for(int ipid=1;ipid<np;++ipid){
      MPI_Recv(&mat(0,0), Nlocal*N, MPI_DOUBLE, ipid, tag, MPI_COMM_WORLD, &status );
  std::cout<<mat<<"\n";
                          	//printmat(mat,N,Nlocal);
    }
  }else{
      MPI_Send(&mat(0,0), Nlocal*N, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }

  
}

