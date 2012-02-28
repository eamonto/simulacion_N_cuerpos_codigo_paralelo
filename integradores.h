//////////////////////////RUNGE-KUTTA4/////////////////////////////
int RK4(int *par_proc);

//////////////////////SIMULATION OF DYNAMIC OF THE STAR//////////////////////
int RK4_integration(int *par_proc);

///////////////////SYMPLECTIC INTEGRATOR//////////////////
int symplectic_integration2(int *par_proc);

///////////////////SYMPLECTIC INTEGRATOR//////////////////
int symplectic_integration(int *par_proc);

//////////CALCULATE ARRAY FOR NUMBER OF PARTICLES PER PROCESATOR////////
int Num_par_proc_division(int *par_proc);

////////THE BOSS OF INTEGRATOR AND EXIT DATES//////////
int emulator(void);


