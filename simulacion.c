/*
========================================================================
simulacion.c
========================================================================
Programa principal

    Copyright (C) 2008  Edison Montoya, eamonto@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


Up to date: 28 Feb 2012					
*/


#include <header.h>


int main (int argc,char **argv)
{
  extern char *inputfile,*initfile,*outputfile;
  int input_array_int[10];
  double input_array[10];

  int error;
 
  error=MPI_Init (&argc, &argv);  //MPI INITIALIZE
  com=MPI_COMM_WORLD;
   
  if(error!=MPI_SUCCESS){
    printf("MPI could not be initialized.\n");
    exit(1);
  }

  MPI_Comm_size (com, &size);
  MPI_Comm_rank (com, &rank);

  if(argc != 4) usage();     //VERIFICATION OF INPUT FILES

  if (rank == 0)
    {
      initfile=argv[1];      //ASIGNATION
      inputfile=argv[2];     //OF FILE'S 
      outputfile=argv[3];    //NAMES
      
      create_remove_dir();   //REMOVE OLD outputfile AND CREATE A NEW ONE

      read_in(input_array_int,input_array);         //READ INITIAL CONDITIONS
    }

  MPI_Bcast(input_array_int, 2 , MPI_INT, 0 , com);  
  MPI_Bcast(input_array, 6 , MPI_DOUBLE, 0 , com);  

  Num_par = input_array_int[0];      //NUMBER OF PARTICLES IN THE SIMULATION
  INTEGRATOR = input_array_int[1];   //INTEGRATOR ELECTION
  
  t = input_array[0];                //TIMESTEP
  G = input_array[1];                //GRAVITATION CONSTANT
  w = input_array[2];                //STAR'S ROTATION 
  Num_archive   = input_array[3];    //NUMBER OUTPUT ARCHIVES
  Num_simulation= input_array[4];    //NUMBER OF INTEGRATION FOR ARCHIVE
  EPS = input_array[5];              //SOFTENING
  
  allocate_memory();                 //ALLOCATE THE MEMORY
  
  if (rank == 0) read_input();       //READ INITIAL POSITIONS AND VELOCITIES 

  Bcast_part_struct(0);  //MASTER MAKE A BROADCAST OF INITIAL POSITIONS AND VELOCITIES

  emulator();            //SIMULATION 

  free_memory();         //FREE THE MEMORY

  error=MPI_Finalize();  //END OF MPI

  return error;
}
