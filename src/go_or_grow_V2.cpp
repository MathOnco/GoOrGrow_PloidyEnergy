/*
    go_or_grow.cpp
    Author: Gregory J. Kimmel
    Date: 06/05/2020

    This file solves the coupled advection-diffusion equation that governs cell
    movement with ploidy:

    E_t = d(eta)*nabla^2 E + h(E)
    u_t = -nabla * (u nabla f(E)) + g(E,k) u
    v_t = -nabla * (u nabla f(E)) + g(E,k) v

    USAGE
        ./go_or_grow -i input.in -o output -v
        ./go_or_grow -v
        ./go_or_grow

    INPUTS
        Read from terminal:
            -v          Verbosity toggle (-v prints details to terminal)
            -i          Params read from file (otherwise default is used)
        Read from input file (if -i is given):
            tFinal:     Max final time to run the simulation to.
            R:          Radius of the dish
            R0:         Initial seeding radius
            dt:         Time step size
            dr:         Spatial step size
            d(eta):     Diffusion coefficient
            a:          Consumption rate
            xi_u:       Tendency for u to move (high = low movement)
            phi_u:      Impact of energy on u (high = sensitive to levels)
            xi_v:       Tendency for v to move (high = low movement)
            phi_v:      Impact of energy level on v (high = sensitive to levels)
            dim:        Dimension of system (n = 2 or 3 is standard)

            growOutsideR0:  Are cells allowed to grow outside seeding density?
            randomCellMotion: Random cell motility 
            alpha:      (OPTIONAL) - Controls the numerical solver
                        (alpha = 1: Backward Euler, alpha = 1/2 Crank-Nicolson)
            
            **THIS FILE IS ALWAYS MADE!!!**
            statisticsFile: Filename for stats

            saveFiles:  Set to true to allow variable files to be made
            outFile:    Parent name of files (_<var>.dat)          

    OUTPUTS
        Written to terminal (if verbosity>0):
            Various diagnostics (e.g. laplacian, gradient, matrices, etc.)
        Written to terminal:
            The elapsed time
        Written to file:
            statisticsFile: Contains biological parameters and important values
                        obtained in the simulation (subject to change).
        Written to files (if saveFiles = true):
            t: The time vector  (Nt x 1)
            r: The space vector (Nr x 1)
            E: The energy       (Nr x Nt)
            u: The u-cell type  (Nr x Nt)
            v: The v-cell type  (Nr x Nt)

 */

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <tuple>
#include <armadillo>
#include <ctime>

using namespace std;
using namespace arma;

// Details of the functions are found after main()
mat make_laplacian(int N, double dx);       // Make laplacian matrix
mat make_gradient(int N, double dx, vec r, int dim); // Make gradient matrix
mat make_gradient(int N, double dx);
vec f(vec E, double xi);                    //Chemotactic term
vec g(vec E, vec u, vec v, double k, vec c);       // Growth term
vec h(vec u, vec v, double a);              // Consumption term
void writeToFile(string filename, vec v);   // Write vector to file
void writeToFile(string filename, mat m);   // Write matrix to file

// Set initial conditions and make solution matrices
tuple<mat, mat, mat> setInitialConditions(int nrows, int ncols, vec r, double R,double R0, double u0, double v0);

// Write output to terminal
void writeToTerminal(string name, string str, int verbNeeded, int verbosity);
void writeToTerminal(string name, int i, int verbNeeded, int verbosity);
void writeToTerminal(string name, double d, int verbNeeded, int verbosity);
void writeToTerminal(string name, vec v, int verbNeeded, int verbosity);
void writeToTerminal(string name, mat M, int verbNeeded, int verbosity);

// Get statistics for output (e.g. mean, median, variance, skewness)
void getStatistics(vec x, vec r, double R, int Nr, double x0,
double total, double &median, double &mean, double &variance,
double &skewness, double &waveEdge, double &waveMid);

// Check to see if file has been created
bool is_file_exist(const char *filename)
{
    return static_cast<bool>(ifstream(filename));
}
bool is_file_exist(string filename)
{
    return static_cast<bool>(ifstream(filename));
}

int main(int argc, char* argv[])
{

    // Start the clock
    clock_t start = clock();

    // The user provides the final time and step size
    // Check if correct usage is given
    if(argc > 5)
    {
        printf("Sample usage: ./go_or_grow -v num -i input.in\n");
        return -1;
    }

    // Check for verbosity. Defaulted to 0
    int verbosity = 0;
    for ( int ii=1; ii<argc; ii++ ) {
        if ( strcmp(argv[ii], "-v") == 0 ) {
            // Verbosity flag.
            ii++;
            verbosity = stoi(argv[ii]);
        }
    }

    // Default parameters
    double tFinal = 100.0;          // Max time to run simulation.
    double R = 10.0;                // Radius of the dish
    double R0 = 2.0;                // Initial seeding radius
    double dt = 0.01;               // time step size
    double dr = 0.1;                // spatial step size
    double u0 = 0.1;                // Initial u0 density
    double v0 = 0.1;                // Initial v0 density
    double eta = 100.0;             // diffusion coefficient
    double a = 3.5;                 // consumption rate
    double xi_u = 0.8;              // tendency for u to move
    double phi_u = 0.25;            // u-sensitivity to low E
    double xi_v = 0.6;              // tendency for v to move
    double phi_v = 0.4;             // v-sensitivity to low E
    double alpha = 0.5;             // (alpha = 0.5,1 for C-N/BE, respectively)
    int dim = 2;                    // Spatial dimension
    bool growOutsideR0 = true;      // If can cells grow outside initial radius?
    bool saveFiles = true;
    double randomCellMotion = 1e-2;

    // Get current date
    time_t rawtime;
    struct tm * timeinfo;
    char buffer [80];
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (buffer,80,"%F",timeinfo);
    string dateString = buffer;

    // default outName and statistics file
    string outFile = "out_" + dateString;
    string statisticsFile = "statistics_" + dateString + ".dat";

    // Process I/O information and read the input if provided.
    string filepath = "";
    for ( int ii=1; ii<argc; ii++ ) {
        if ( strcmp(argv[ii], "-i") == 0 ) {
            // Input file provided.

            // Increment ii and read the filename/path from the next input.
            ii++;
            filepath = argv[ii];
            writeToTerminal("Reading input file", filepath, 1, verbosity);

            // Open the input file.
            string line;
            string varName;
            string varValue;
            ifstream inputFile(filepath);
            if ( inputFile.is_open() ) {
                // Extract data from each line.
                size_t found;
                while ( getline(inputFile,line) )
                {
                    // Find the first tab and grab var name and value.
                    found = line.find_first_of("\t");
                    varName = line.substr(0, found);
                    varValue = line.substr(found+1);

                    // Figure out which param is getting set and set it.
                    if( "outFile" == varName)
                    {
                        outFile = varValue;
                        writeToTerminal("Output name", outFile,
                        1, verbosity);
                    }

                    if( "statisticsFile" == varName)
                    {
                        statisticsFile = varValue + ".dat";
                        writeToTerminal("Statistics file", 
                        statisticsFile,1, verbosity);
                    }

                    if( "tFinal" == varName)
                    {
                        tFinal = stod(varValue);
                        writeToTerminal("tfinal",tFinal,1, verbosity);
                    }

                    if( "R" == varName)
                    {   
                        R = stod(varValue);
                        writeToTerminal("R",R,1, verbosity);
                    }
                    
                    if( "R0" == varName)
                    {
                        R0 = stod(varValue);
                        writeToTerminal("R0",R0,1,verbosity);
                    }
                    
                    if( "dt" == varName)
                    {
                        dt = stod(varValue);
                        writeToTerminal("dt",dt,1, verbosity);
                    }
                    
                    if( "dr" == varName)
                    {
                        dr = stod(varValue);
                        writeToTerminal("dr", dr,1, verbosity);
                    }
                    
                    if( "u0" == varName)
                    {
                        u0 = stod(varValue);
                        writeToTerminal("u0", u0,1, verbosity);
                    }
                    
                    if( "v0" == varName)
                    { 
                        v0 = stod(varValue);
                        writeToTerminal("v0",v0,1, verbosity);
                    }

                    if( "eta" == varName)
                    {
                        eta = stod(varValue);
                        writeToTerminal("eta",eta,1, verbosity);
                    }

                    if( "a" == varName)
                    {
                        a = stod(varValue);
                        writeToTerminal("a",a,1, verbosity);
                    }

                    if( "xi_u" == varName)
                    {
                        xi_u = stod(varValue);
                        writeToTerminal("xi_u",xi_u,1, verbosity);
                    }

                    if( "phi_u" == varName)
                    { 
                        phi_u = stod(varValue);
                        writeToTerminal("phi_u",phi_u,1, verbosity);
                    }

                    if( "xi_v" == varName)
                    {
                        xi_v = stod(varValue);
                        writeToTerminal("xi_v",xi_v,1, verbosity);
                    }

                    if( "phi_v" == varName)
                    {
                        phi_v = stod(varValue);
                        writeToTerminal("phi_v", phi_v,1, verbosity);
                    }

                    if( "alpha" == varName)
                    {
                        alpha = stod(varValue);
                        writeToTerminal("alpha",alpha,1, verbosity);
                    }

                    if( "dim" == varName)
                    {
                        dim = stoi(varValue);
                        writeToTerminal("dim",dim,1, verbosity);
                    }

                    if( "growOutsideR0" == varName)
                    {
                        if(varValue == "true" || varValue == "True" ||
                        varValue == "t")
                            growOutsideR0 = true;
                        else
                            growOutsideR0 = false;
                        writeToTerminal("growOutsideR0",growOutsideR0,1,
                        verbosity);
                    }

                    if( "saveFiles" == varName)
                    {
                        if(varValue == "true" || varValue == "True" ||
                        varValue == "t")
                            saveFiles = true;
                        else
                            saveFiles = false;
                        writeToTerminal("saveFiles",saveFiles,1,
                        verbosity);
                    }

                    if( "randomCellMotion" == varName)
                    {
                        randomCellMotion = stod(varValue);
                        writeToTerminal("randomCellMotion",
                        randomCellMotion,1, verbosity);
                    }
                }

                // Close the input value.
                inputFile.close();
            }
            else
            {
                // File was not read, abort.
                cout << "NO INPUT FILE READ!" << endl;
                return -1;
            }
        }
    }

    // Create variable file names
    string tfile = outFile + "_t.dat";
    string rfile = outFile + "_x.dat";
    string Efile = outFile + "_E.dat";
    string ufile = outFile + "_u.dat";
    string vfile = outFile + "_v.dat";

    // Number of temporal and spatial steps.
    // NOTE: An extra 1 is added (the tFinal/dt gives total number of steps
    // taken but does not include the first step t = 0 or r = 0!!)
    int Nt = tFinal/dt + 1;
    int Nr = R/dr + 1;

    // Discretize time and space
    vec t = linspace(0.0,tFinal,Nt);
    vec r = linspace(0.0,R,Nr);

    // Output time and space vectors for error checking
    writeToTerminal("t", t, 10, verbosity);
    writeToTerminal("r", r, 10, verbosity);

    // Make laplacian, gradient and I matrix
    mat laplacian(Nr,Nr), gradLaplace(Nr,Nr), gradient(Nr,Nr), I = eye(Nr,Nr);
    laplacian = make_laplacian(Nr, dr);
    gradLaplace = make_gradient(Nr,dr,r,dim);
    gradient = make_gradient(Nr,dr);

    // Output Gradient and Laplacian to terminal for error checking
    writeToTerminal("Laplacian Matrix", laplacian, 11, verbosity);
    writeToTerminal("Gradient Matrix", gradient, 11, verbosity);

    // We now build our solver using Crank-Nicolson so we will create
    // LHS*x^(n+1) = RHS*x^(n) and then call a solve routine and iterate.
    // The LHS and RHS depend on the iteration so there is no benefit to
    // creating them before beginning our iteration
    mat LHS(Nr,Nr), RHS(Nr,Nr), z(Nr,Nr);

    // make solution vector and set initial conditions
    mat E,u,v;
    tie(E,u,v) = setInitialConditions(Nr,Nt,r,R,R0,u0,v0);

    // Test initial conditions
    writeToTerminal("E(r,0)",E,10,verbosity);
    writeToTerminal("u(r,0)",u,10,verbosity);
    writeToTerminal("v(r,0)",v,10,verbosity);

    // Make commonly used matrices
    mat radL, temp(Nr,Nr,fill::zeros);
    radL = laplacian + gradLaplace;

    // Used to generate time-dependent part of matrices
    vec f0,g0,h0, soln, growthAllowed(Nr,fill::ones);

    if (!growOutsideR0)
    {
        for(int i=0;i<Nr;i++)
        {
            if (r(i) > R0)
                growthAllowed(i) = 0.0;
        }
    }

    cout << growthAllowed << endl;

    // We generate

    // Exit if the sum(E) < tol (no energy left)
    double tol = 1e-3;
    bool energyLeft = true;
    bool steadyState = false;

    // Begin iteration
    for(int n=0;n<Nt-1;n++)
    {

        // Build matrix solve for E
        h0 = h(u.col(n), v.col(n), a);

        // z := eta*nabla^2 + h0
        z = eta*radL + diagmat(h0);
        LHS = I - dt*alpha*z; RHS = I + dt*(1.0-alpha)*z;

        // update E
        soln = solve(LHS, RHS*E.col(n));

        // Remove small negative numbers introduced by num. error.
        soln.transform( [](double val) { return (val < 0.0) ? 0.0 : val; } );
        // soln.clean(datum::eps);
        E.col(n+1) = soln;

        if(abs(soln(Nr-1)-soln(Nr-2))>dr*dr)
        {
            cout << "dE/dr large near edge. Increase cell movement." << endl;
            // randomCellMotion *= 2.0;
        }

        // Build matrix solve for u
        f0 = f(E.col(n),xi_u);
        g0 = g(E.col(n),u.col(n),v.col(n),phi_u,growthAllowed);

        // z = eps*nabla^2 - nabla^2(f(E)) - nabla(f(E))*nabla + g0
        z = randomCellMotion*radL + diagmat(g0 - radL*f0) -
        diagmat(gradient*f0)*gradient;
        
        LHS = I - dt*alpha*z; RHS = I + dt*(1.0-alpha)*z;

        // update u
        soln = solve(LHS, RHS*u.col(n));

        // Remove small negative numbers introduced by num. error.
        soln.transform( [](double val) { return (val < 0.0) ? 0.0 : val; } );
        // soln.clean(datum::eps);
        u.col(n+1) = soln;

        if(abs(soln(Nr-1)-soln(Nr-2))>dr*dr)
        {
            cout << "du/dr large near edge. Increase cell movement." << endl;
            // randomCellMotion *= 2.0;
        }

        // Build matrix solve for v
        f0 = f(E.col(n),xi_v);
        g0 = g(E.col(n),u.col(n),v.col(n),phi_v, growthAllowed);

        // z = eps*nabla^2 - nabla^2(f(E)) - nabla(f(E))*nabla + g0
        z = randomCellMotion*radL + diagmat(g0 - radL*f0) -
        diagmat(gradient*f0)*gradient;
        
        LHS = I - dt*alpha*z; RHS = I + dt*(1.0-alpha)*z;

        // update v
        soln = solve(LHS, RHS*v.col(n));

        // Remove small negative numbers introduced by num. error.
        soln.transform( [](double val) { return (val < 0.0) ? 0.0 : val; } );
        // soln.clean(datum::eps);
        v.col(n+1) = soln;

        if(abs(soln(Nr-1)-soln(Nr-2))>dr*dr)
        {
            cout << "dv/dr large near edge. Increase cell movement." << endl;
            // randomCellMotion *= 2.0;
        }

        // Check if total energy is small
        if (sum(E.col(n+1))/Nr < tol)
        {
            Nt = n+1; // New total time length vector
            energyLeft = false;
            break;
        }

        // // Check if system is in equilibrium dF/dt ~ 0
        // if(sum(abs(u.col(n+1) - u.col(n)))/Nr < dt*tol &&
        // sum(abs(v.col(n+1) - v.col(n)))/Nr < dt*tol)
        // {
        //     Nt = n+1; // New total time length vector
        //     steadyState = true;
        //     break;
        // }

    }

    // Resize data to remove extra elements (not run because energy was gone).
    E.resize(Nr,Nt);
    u.resize(Nr,Nt);
    v.resize(Nr,Nt);
    t.resize(Nt);

    // Output final time
    tFinal = t(Nt-1);
    writeToTerminal("Final time",tFinal,1,verbosity);

    // Output for error checking
    writeToTerminal("E",E,5,verbosity);
    writeToTerminal("u",u,5,verbosity);
    writeToTerminal("v",v,5,verbosity);

    // Save variables if desired
    if(saveFiles==1)
    {
        writeToFile(tfile,t);
        writeToFile(rfile,r);
        writeToFile(Efile,E);
        writeToFile(ufile,u);
        writeToFile(vfile,v);
    }

    // Get total populations from last time point
    vec uEnd(Nr), vEnd(Nr);
    uEnd = u.col(Nt-1); vEnd = v.col(Nt-1);
    double uTot = (R/Nr)*sum(uEnd % r);
    double vTot = (R/Nr)*sum(vEnd % r);
    double Etot = (R/Nr)*sum(E.col(Nt-1) % r);

    // Get useful statistics
    double uMedian, uMean, uVar, uSkew, uWaveEdge, uWaveMid,
    vMedian, vMean, vVar, vSkew, vWaveEdge, vWaveMid;

    // u-Stats and v-Stats
    getStatistics(uEnd, r, R, Nr, u0, uTot, uMedian, uMean, uVar,
    uSkew, uWaveEdge, uWaveMid);
    getStatistics(vEnd, r, R, Nr, v0, vTot, vMedian, vMean, vVar,
    vSkew, vWaveEdge, vWaveMid);
  
    // Check to see if file has been created
    bool fileExists = is_file_exist(statisticsFile);

    // Output file containing all parameters and interesting quantities to a
    // tab-delimited file (with header).
    ofstream paramFile;
    paramFile.open(statisticsFile, ofstream::out | ofstream::app);

    // If file does not exist we give the header
    if(!fileExists)
    {
        paramFile <<
        "#dt\tdr\tT\tR\tR0\tu0\tv0\teta\ta\txi_u\tk_u\txi_v\t" <<
        "k_v\tu_{tot}\tv_{tot}\tu_{edge}\tv_{edge}\tu_{mid}\tv_{mid}\tu_{med}"
        << "\tv_{med}\t" << "u_{mean}\tv_{mean}\t" <<
        "u_{var}\tv_{var}\tu_{skew}\tv_{skew}\tEnergyLeft?\tSteadyState? \n" <<
        "#------------------------------------------------------------" <<
        "-------------------------------------------------------------" <<
        "------------------------------------------------------------------" << 
        "------------------------------------------------------------------" << 
        "-----------" << endl;
    }
    paramFile << dt << "\t" << dr << "\t" << tFinal << "\t" << R << "\t"
    << R0 << "\t" << u0 << "\t" << v0 << "\t" << eta << "\t" << a << "\t" 
    << xi_u << "\t" << phi_u << "\t" << xi_v << "\t" << phi_v << "\t" << uTot
    << "\t"<< vTot << "\t" << uWaveEdge << "\t"<< vWaveEdge << "\t" <<
    uWaveMid << "\t" << vWaveMid << "\t" << uMedian << "\t" << vMedian << "\t"
    << uMean << "\t" << vMean << "\t" << uVar << "\t" << vVar << "\t" << uSkew
    << "\t" << vSkew << "\t" << energyLeft << "\t" << steadyState << endl;

    // Close the file
    paramFile.close();


    // Print the time and return.
    if(verbosity>1)
        printf("Time elapsed: %g seconds\n", (float)(clock()-start)/
        CLOCKS_PER_SEC);

    return 0;
}

/*

    make_laplacian(N,dx)

    This function creates the matrix that approximates d^2/dr^2 with No-flux
    boundary conditions
        
        INPUTS
            N:  The size of the matrix
            dx: The step size

        OUTPUTS
            laplacian: The laplacian matrix

*/
mat make_laplacian(int N, double dx)
{

    // Create matrix
    mat laplacian(N,N,fill::zeros);

    // We fill the matrix by indexing
    for(int i=0;i<N;i++)
    {

        // Don't fill past the boundary (sub/super diagonals)
        if(i<N-1)
            laplacian(i,i+1) = 1.0;
        if(i>0)
            laplacian(i,i-1) = 1.0;
        // Main diagonal
        laplacian(i,i) = -2.0;

    }

    // Boundary conditions (No-flux)
    laplacian(0,1) = 2.0;
    laplacian(N-1,N-2) = 2.0;

    // laplacian(0,0) = -1.0;
    // laplacian(0,1) = 1.0;
    // laplacian(N-1,N-2) = 1.0;
    // laplacian(N-1,N-1) = -1.0;

    // Divide by dx^2 for approximating the laplacian
    return laplacian/(dx*dx);

}

/*

    make_gradient(N,dx,r,dim)

    This function creates the matrix that approximates (1/r)*d/dr with No-flux
    boundary conditions
        
        INPUTS
            N:  The size of the matrix
            dx: The step size
            r:  The spatial vector

        OUTPUTS
            gradient: The gradient matrix

*/
mat make_gradient(int N, double dx, vec r, int dim)
{

    // Create matrix
    mat gradient(N,N,fill::zeros);

    int i;
    // Fill matrix (we skip the edges because no-flux implies the gradient is
    // zero there).
    for(i=1;i<N-1;i++)
    {
            gradient(i,i-1) = -(dim - 1)/r[i];
            gradient(i,i+1) = (dim - 1)/r[i];
    }

    return gradient/(2.0*dx);

}

/*

    make_gradient(N,dx)

    This function creates the matrix that approximates d/dr with No-flux
    boundary conditions
        
        INPUTS
            N:  The size of the matrix
            dx: The step size

        OUTPUTS
            gradient: The gradient matrix

*/
mat make_gradient(int N, double dx)
{

    // Create matrix
    mat gradient(N,N,fill::zeros);

    int i;
    // Fill matrix (we skip the edges because no-flux implies the gradient is
    // zero there).
    for(i=1;i<N-1;i++)
    {
            gradient(i,i-1) = -1.0;
            gradient(i,i+1) = 1.0;
    }

    return gradient/(2.0*dx);

}

/*

    f(E,xi)

    This function outputs the vector f(E,xi) - the chemotactic term
    Currently we use log(xi + E), but other functional forms can be given
        
        INPUTS
            E:  The amount of energy
            xi: The tendency to move (high xi = low movement)

        OUTPUTS
            f:  A vector outputting the values for a given energy.

*/
vec f(vec E, double xi)
{

    return log(xi + E);

}

/*

    g(E, u, v,k, c)

    This function outputs the vector g(E, u, v,k) - the growth term
    Currently we use E/(k + E), but other functional forms can be given
        
        INPUTS
            E:  The amount of energy
            u:  The amount of goers
            v:  The amount of growers
            k: The sensitivity to less energy (high k = high sensitivity)
            c: Vector of ones and 0 that restricts growth if ECM prohibits it

        OUTPUTS
            g:  A vector outputting the values for a given energy.

*/
vec g(vec E, vec u, vec v, double k, vec c)
{

    if(k==0.0)
        return (c - u - v);
    else
        return (E/(k+E)) % (c - u - v);

}

/*

    h(u,v,a)

    This function outputs the vector h(u,v,a) - the consumption term
    Currently we use -a(u+v), but other functional forms can be given
        
        INPUTS
            u:  The amount of cell u
            v:  The amount of cell v
            a:  Consumption rate

        OUTPUTS
            h:  A vector outputting the total consumption.

*/
vec h(vec u, vec v, double a)
{
    return -a*(u + v);
}



void writeToTerminal(string name, string str, int verbNeeded, int verbosity)
{
    if(verbosity >= verbNeeded)
    {
        cout << name << ": " << str << endl;
    }
}

/*

    writeToTerminal(name, variable, verbNeeded, verbosity)

    This overloaded function writes the given variable to terminal. Mostly used
    for error-checking or more interactive runs.

        INPUTS
            name:       Variable name
            variable:   the quantity to be outputted
            verbNeeded: The minimum verbosity number needed to print to terminal
            verbosity:  User-defined determines whether we output info to term
        
        OUTPUTS
            NONE

*/
void writeToTerminal(string name, mat M, int verbNeeded, int verbosity)
{

    if(verbosity >= verbNeeded)
    {
        cout << name << ":" << endl;
        cout << M << endl;
    }

}

void writeToTerminal(string name, vec v, int verbNeeded, int verbosity)
{

    if(verbosity >= verbNeeded)
    {
        cout << name << ":" << endl;
        cout << v << endl;
    }

}

void writeToTerminal(string name, int i, int verbNeeded, int verbosity)
{

    if(verbosity >= verbNeeded)
        cout << name << ": " << i << endl;

}

void writeToTerminal(string name, double d, int verbNeeded, int verbosity)
{
    if(verbosity >= verbNeeded)
    {
        cout << name << ": " << d << endl;
    } 
}

/*

    setInitialConditions(nrows, ncols, r, R, R0, u0, v0)

    This function sets initial conditions for the matrices E, u, v our solution
    that we need.

        INPUTS
            nrows:      Number of matrix rows (Nr)
            ncols:      Number of matrix columns (Nt)
            r:          The spatial vector
            R:          Radius of the dish
            R0:         Seeding radius
            u0:         Seeding density of cell type u
            v0:         Seeding density of cell type v
        
        OUTPUTS
            E:          Energy matrix with initial condition filled.
            u:          Cell u type with initial condition filled.
            v:          Cell v type with initial condition filled

*/
tuple<mat, mat, mat> setInitialConditions(int nrows, int ncols, vec r, double R,
double R0, double u0, double v0)
{

    // Make matrices
    mat E(nrows,ncols), u(nrows,ncols), v(nrows,ncols);

    // Iterate through space to set the initial conditions
    for(int n=0;n<nrows;n++)
    {
        // Assume E is uniformly set in space
        E(n,0) = 1.0;

        // Seed u, v if r < R0
        u(n,0) = (r[n]<=R0 ? u0 : 0);
        v(n,0) = (r[n]<=R0 ? v0 : 0);
    }

    return make_tuple(E, u, v);

}

/*

    writeToFile(filename, variable)

    This overloaded function writes the given variable to the file "filename".

*/
void writeToFile(string filename, vec v)
{
    ofstream myfile;
    myfile.open(filename, ofstream::out);
    myfile << v;
    myfile.close();
}

void writeToFile(string filename, mat M)
{
    ofstream myfile;
    myfile.open(filename, ofstream::out);
    myfile << M;
    myfile.close();
}

void getStatistics(vec x, vec r, double R, int Nr, double x0,
double total, double &median, double &mean, double &variance,
double &skewness, double &waveEdge, double &waveMid)
{

    waveMid = 0.0;
    waveEdge = 0.0;

    // Initialize vectors
    vec pdf(Nr), cdf(Nr);

    // Make PDF and CDF
    pdf = x/total;
    cdf = (R/Nr)*cumsum(pdf % r);

    // Calculate statistical values
    mean = (R/Nr)*sum(pdf % r % r);
    variance = (R/Nr)*sum(pdf % r % r % r) - mean*mean;
    skewness = ((R/Nr)*sum(pdf % r % r % r % r) - 3.0*mean*variance -
                pow(mean,3))/pow(variance,1.5);
    median = r(index_min(abs(cdf - 0.5)));

    // Get degree of infiltration
    if(x0>0)
    {
        waveMid = as_scalar(r(index_min((x-0.5)%(x-0.5))));
        waveEdge = as_scalar(r(index_min((x-0.01)%(x-0.01))));
    }

}