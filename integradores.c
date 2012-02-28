/*
=========================================================================
integradores.c
=========================================================================
Esta libreria contiene la implementacion de un integrador simplectico 
de segundo orden (leapfrog) para N-cuerpos, y un integrador Runge-Kutta 4.


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

//////////////////////////RUNGE-KUTTA4/////////////////////////////
int RK4(int *par_proc)
{
  int i,j,min,max,num_par_proc[size+1];
  double acele[3],auxF[MAXDIM][Num_par],auxG[MAXDIM][Num_par];

  min=par_proc[rank];
  max=par_proc[rank+1];

  for(i=0;i<size;i++) 
    num_par_proc[i]=par_proc[i+1]-par_proc[i];  //PARTICLES IN THE PROSESSORS


//////////////CALCULA K1/////////////////
  for(j=min;j<max;j++)
    {
      aceleration(j,acele);

      for(i=0;i<MAXDIM;i++)
	{
	  auxF[i][j-min]=t*acele[i];
	  auxG[i][j-min]=t*par[j].vel[i];
	}
    }

///////ALL PROCESS KNOW THE FUNCTIONS F1 AND G1///////
  for(i=0;i<MAXDIM;i++)
    {
      MPI_Allgatherv(auxF[i],max-min,MPI_DOUBLE,F1[i],num_par_proc,par_proc,MPI_DOUBLE,com);
      MPI_Allgatherv(auxG[i],max-min,MPI_DOUBLE,G1[i],num_par_proc,par_proc,MPI_DOUBLE,com);
    }

//////////////CALCULA K1/////////////////



//////////////CALCULA K2/////////////////
  for(j=0;j<Num_par;j++)
    {
      for(i=0;i<MAXDIM;i++) //NEW POSITION FOR CALCULATE K2
	{
	  par[j].pos[i]=Auxiliar[i][j]+0.5*G1[i][j];
	}
    }


  for(j=min;j<max;j++)
    {
      aceleration(j,acele);

      for(i=0;i<MAXDIM;i++)
	{
	  auxF[i][j-min]=t*acele[i];
	  auxG[i][j-min]=t*(par[j].vel[i]+0.5*F1[i][j]);
	}
    }

//////////ALL PROCESS KNOW THE FUNCTIONS F2 AND G2////////// 
  for(i=0;i<MAXDIM;i++)
    {
      MPI_Allgatherv(auxF[i],max-min,MPI_DOUBLE,F2[i],num_par_proc,par_proc,MPI_DOUBLE,com);
      MPI_Allgatherv(auxG[i],max-min,MPI_DOUBLE,G2[i],num_par_proc,par_proc,MPI_DOUBLE,com);
    }

//////////////CALCULA K2/////////////////



//////////////CALCULA K3/////////////////

  for(j=0;j<Num_par;j++)
    {
      for(i=0;i<MAXDIM;i++)  //NEW POSITION FOR CALCULATE K3
	{
	  par[j].pos[i]=Auxiliar[i][j]+0.5*G2[i][j];
	}
    }

  for(j=min;j<max;j++)
    {
      aceleration(j,acele);

      for(i=0;i<MAXDIM;i++)
	{
	  auxF[i][j-min]=t*acele[i];
	  auxG[i][j-min]=t*(par[j].vel[i]+0.5*F2[i][j]);
	}
    }

//////////ALL PROCESS KNOW THE FUNCTIONS F3 AND G3////////// 
  for(i=0;i<MAXDIM;i++)
    {
      MPI_Allgatherv(auxF[i],max-min,MPI_DOUBLE,F3[i],num_par_proc,par_proc,MPI_DOUBLE,com);
      MPI_Allgatherv(auxG[i],max-min,MPI_DOUBLE,G3[i],num_par_proc,par_proc,MPI_DOUBLE,com);
    }

//////////////CALCULA K3/////////////////



//////////////CALCULA K4/////////////////

  for(j=0;j<Num_par;j++)
    {
      for(i=0;i<MAXDIM;i++)  //NEW POSITION FOR CALCULATE K4
	{
	  par[j].pos[i]=Auxiliar[i][j]+G3[i][j];
	}
    }


  for(j=min;j<max;j++)
    {
      aceleration(j,acele);

      for(i=0;i<MAXDIM;i++)
	{
	  auxF[i][j-min]=t*acele[i];
	  auxG[i][j-min]=t*(par[j].vel[i]+F3[i][j]);
	}
    }

//////////ALL PROCESS KNOW THE FUNCTIONS F4 AND G4////////// 
  for(i=0;i<MAXDIM;i++)
    {
      MPI_Allgatherv(auxF[i],max-min,MPI_DOUBLE,F4[i],num_par_proc,par_proc,MPI_DOUBLE,com);
      MPI_Allgatherv(auxG[i],max-min,MPI_DOUBLE,G4[i],num_par_proc,par_proc,MPI_DOUBLE,com);
    }

  return 0;
}


//////////////////////RK4 INTEGRATOR//////////////////////
int RK4_integration(int *par_proc)
{
  int i,j;
  double C,D,min,max;

  ////MAKE A COPY OF THE ORIGINAL POSITIONS
  for(j=0;j<Num_par;j++){
    for(i=0;i<MAXDIM;i++){
      Auxiliar[i][j]=par[j].pos[i];
    }
  }

  ////CALCULATE OF F's AND G's FOR ALL PARTICLE IN ALL DIRECTION
  RK4(par_proc);  

  min=par_proc[rank];
  max=par_proc[rank+1];

  //////////// PARTICLE j; VELOCITY AND POSITION IN DIRECTION i ////////
  for(j=min;j<max;j++)
    {
      for(i=0;i<MAXDIM;i++)
	{
	  C = F1[i][j]+2.0*F2[i][j]+2.0*F3[i][j]+F4[i][j];
	  par[j].vel[i] = par[j].vel[i]+C/6.0;

	  D = G1[i][j]+2.0*G2[i][j]+2.0*G3[i][j]+G4[i][j];
	  par[j].pos[i] = Auxiliar[i][j]+D/6.0;
	}
    }

  //POSITION'S BROADCAST 
  for(i=0;i<size;i++)
    Bcast_part_pos(i,par_proc[i],par_proc[i+1]);

  //VELOCITIES BROADCAST
  for(i=0;i<size;i++)
    Bcast_part_vel(i,par_proc[i],par_proc[i+1]);

  return 0;
}


///////////////////SYMPLECTIC INTEGRATOR//////////////////
int symplectic_integration2(int *par_proc)
{
  int i,j,min,max;
  double acele[3];

  min=par_proc[rank];
  max=par_proc[rank+1];

  //FIRST INTEGRATION OF POSITIONS
  for(j=0;j<Num_par;j++){
    for(i=0;i<MAXDIM;i++){
      par[j].pos[i] = par[j].pos[i] + 0.5*t*par[j].vel[i];
    }
  }


  //INTEGRATION OF VELOCITIES
  for(j=min;j<max;j++){

    aceleration(j,acele);
    
    for(i=0;i<MAXDIM;i++){
      par[j].vel[i] = par[j].vel[i] + t*acele[i];
    }
  }

  //VELOCITIES BROADCAST
  for(i=0;i<size;i++)
    Bcast_part_vel(i,par_proc[i],par_proc[i+1]);


  //SECOND INTEGRATION OF POSITIONS
  for(j=0;j<Num_par;j++){
    for(i=0;i<MAXDIM;i++){
      par[j].pos[i]=  par[j].pos[i] + 0.5*t*par[j].vel[i];
    }
  }

  return 0;
}


///////////////////SYMPLECTIC INTEGRATOR//////////////////
int symplectic_integration(int *par_proc)
{
  int i,j,min,max;
  double acele[3];

  min=par_proc[rank];
  max=par_proc[rank+1];

  //FIRST INTEGRATION OF POSITIONS
  for(j=min;j<max;j++){
    for(i=0;i<MAXDIM;i++){
      par[j].pos[i] = par[j].pos[i] + 0.5*t*par[j].vel[i];
    }
  }

  //POSITION'S BROADCAST 
  for(i=0;i<size;i++)
    Bcast_part_pos(i,par_proc[i],par_proc[i+1]);


  //INTEGRATION OF VELOCITIES
  for(j=min;j<max;j++){

    aceleration(j,acele);
    
    for(i=0;i<MAXDIM;i++){
      par[j].vel[i] = par[j].vel[i] + t*acele[i];
    }
  }

  //VELOCITIES BROADCAST
  for(i=0;i<size;i++)
    Bcast_part_vel(i,par_proc[i],par_proc[i+1]);


  //SECOND INTEGRATION OF POSITIONS
  for(j=min;j<max;j++){
    for(i=0;i<MAXDIM;i++){
      par[j].pos[i]=  par[j].pos[i] + 0.5*t*par[j].vel[i];
    }
  }

  //POSITION'S BROADCAST   
  for(i=0;i<size;i++)
    Bcast_part_pos(i,par_proc[i],par_proc[i+1]);


  return 0;
}



//////////////CALCULATE ARRAY FOR NUMBER OF PARTICLES PER PROCESSOR//////////
int Num_par_proc_division(int *par_proc)
{
  int aux1,aux2,i;

  aux1=Num_par%size;
  aux2=(Num_par-aux1)/size;

  par_proc[0]=0;
  par_proc[1]=aux2+aux1;

  if(rank==0)
    {
      printf("particulas sobrantes: %d\n",aux1);
      printf("particulas por procesador: %d\n",aux2);
      printf("Numero de particulas del primer procesador %d\n",par_proc[1]);
    }

  for(i=2;i<=size;i++) 
    par_proc[i]=par_proc[i-1]+aux2;

  if(par_proc[size]!=Num_par) 
    {  
      printf("\n Problema en la reparticion de las particulas\n");
      exit(1);
    }

  return 0;
}


////////THE BOSS OF INTEGRATOR AND EXIT DATES//////////
int emulator(void)
{
  int i,l,par_proc[size+2];
  double Ttotal,time;
  
  Ttotal = t*Num_simulation*Num_archive;
  if(rank == 0)  printf("\nTotal time of simulation =%lf\n\n",Ttotal);

  time=0.0;

  Num_par_proc_division(par_proc); //PARTICLES FOR PROCESSOR
 

  ////THIS "for" GENERATE "Num_archive" ARCHIVES WHIT POSITION AND VELOCITY
  for(i=0 ; i<Num_archive ; i++)
    {
      for (l=0; l<Num_simulation; l++) ///NUMBER OF SIMULATION IN ONE ARCHIVE
	{ 
	  if(l%10==0)  energy(time,par_proc); //ENERGY CALCULATION

	  if(rank == 0) printf("Time of Simulation: %lf\n",time);

	  if(INTEGRATOR==0)
	    { ////SYMPLECTIC INTEGRATION//
	      symplectic_integration(par_proc);
	    }
	  else if (INTEGRATOR==1)
	    { ////RK4 INTEGRATION ////////
	      RK4_integration(par_proc); 
	    }
	  else
	    { ////SYMPLECTIC INTEGRATION//
	      symplectic_integration2(par_proc);
	    }

	  time=time + t;	
	}

      centermasspos(par_proc); //CENTER MASS POSITIONS
      centermassvel(par_proc); //CENTER MASS VELOCITIES

      if(rank == 0)  write_output(i); //WRITE OUTPUT ARCHIVE i-esimo

    }

  return 0;
}

