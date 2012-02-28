/////// ACELERATION(in variable acele) OF PARTICLE i///////////
int aceleration(int i,double *acele);

/////////////CALCULATE OF GRAVITATIONAL POTENCIAL///////////////
double potcalc(int i);

/////////////////////////////TOTAL ENERGY/////////////////////////
int energy(double time,int *par_proc);

///////////////////CENTER OF MASS POSITION//////////////////////
int centermasspos(int *par_proc);

////////////////////CENTER OF MASS VELOCITY/////////////////////
int centermassvel(int *par_proc);

/////////////IO////////////////

///////////VERIFICATION OF INPUT FILES/////////////
int usage(void);

//REMOVE OLD outputfile AND CREATE A NEW ONE
int create_remove_dir(void);

//////////////////READ PARAMETERS OF THE STAR///////////////////
int read_in(int *input_array_int,double *input_array);

////////////READ INITIAL POSITIONS AND VELOCITYS////////////////
int read_input(void);

/////////WRITE RADIUS, POSITIONS, VELOCITIES AND MASS//////////////
int write_output(int j);

////////ALLOCATE MEMORY////////////
int allocate_memory(void);

/////////FREE MEMORY THAT WAS ALLOCATE IN THE PROGRAM///////////
int free_memory(void);

///////PROCESSOR myrank MAKE A BROADCAST OF POSITIONS////////
int Bcast_part_pos (int myrank,int min, int max);

///////PROCESSOR myrank MAKE A BROADCAST OF VELOCITIES////////
int Bcast_part_vel (int myrank,int min, int max);

/////PROCESSOR myrank MAKE A BROADCAST OF PATICLE STRUCTURE/////
int Bcast_part_struct (int myrank);
