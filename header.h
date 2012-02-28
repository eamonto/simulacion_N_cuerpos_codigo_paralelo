/*
=================================================================
header.h
=================================================================
Es la cabecera del programa, en el se declaran las variables globales y se 
crea la estructura principal del programa "Particle" que contiene la 
informacion de las particulas: posicion, velocidad, masa, distancia al
centro de coordenadas.

El programa se desarrollo con dos integradores distintos un RK4 y metodo
simplectico (leapfrog), ambos metodos necesitan de las rutinas contenidas en 
lib_simulaction.c, los metodos estan implementados en el archivo
integradores.c

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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define MAXDIM 3

///////////////////////// GLOBAL VARIABLES///////////////////////

int Num_par;
int INTEGRATOR;
double t,totalenergy,potencial,kinetic,w,G,EPS;
double Num_archive,Num_simulation; 

char *inputfile,*initfile,*outputfile; //IN-OUT VARIABLES 

double *Auxiliar[MAXDIM]; //AUXILIAR FOR INTEGRATION
double *F1[MAXDIM],*F2[MAXDIM],*F3[MAXDIM],*F4[MAXDIM]; //FOR RK4
double *G1[MAXDIM],*G2[MAXDIM],*G3[MAXDIM],*G4[MAXDIM]; //FOR RK4


//////////PARALELO///////////////
int rank,size;
MPI_Comm com;


/////////////PARTICLE STRUCTURE ////////////////////////
struct Particle
{
  double radius;
  double pos[MAXDIM];
  double vel[MAXDIM];
  double mass; 
};

struct Particle *par;

///////////////LYBRARIES////////////////////////
#include <lib_simulacion.h>
#include <integradores.h>
