#ifndef SIMDATA_H
#define SIMDATA_H

#include <vector>
#include <string>
#include "platformEnv.h"
#include "dataArgs.h"


/*!
  Currently (and possibly permanently) a very hacky way of keeping singleCell cons2prims function
  general for the benefit of the IMEX integrator.
  SCC2P is an external void function pointer that will be set inside the model constructor to
  be the
  __device__ void getPrimitiveVarsSingleCellParallel()
  function. This way the IMEX integrator should be able to handle multiple models.
  Note: the reason this is required is because CUDA does not allow host classes to contain global
  or device member functions, and the integrator requires the member function of
  the model during the rootfind to perform the C2P conversion and to get the source
  vector.
*/
// __device__
// extern void (*SCC2P)(double *, double *, double *, double, double);
// extern void (*SCSource)();

//! <b> Data object </b>
/*!
  @par
    Class contains all the data of the simulation relevant to any of the other
  modules. Containing it in this way prevents issues of cyclic includes, also
  results in Simulation as more of an interface than a class that needs to be
  known to lower objects---good practice. <br>

  Usage <br>
  =========
  @par
    Call the constructor with at least: the number of cells and domain limits in
  the x, y, and z direction, and the simulation end time. As ghost cells are
  automatically added to the domain, there would be a minimum number of total cells
  in any direction of nine---as a result, running 2D simulations would be 9x slower
  than necessary. To combat this, you can request zero cells in the z-direction
  which is equivalent to a 2D domain where Nz=1. Keep this behaviour in mind when
  constructing new models.
  @par
    Other variables, such as the courant factor, can also be set in the constructor
  but by default have sensible values that should work for most set ups, these can
  largely be ignored unless the simulation is failing to converge. <br>
  @par
    Selecting a model will automatically set the number of cons prims and aux vars,
  both inside the model class and in the data class (although it is the data class
  values that are accessed by the various functions). <br>
  @par
    The elementID function data.id(var, i, j, k) is a useful shortcut for accessing
  data belonging to a specific cell, where var is the id of the variable in the array
  we wish to access, and \f$(i, j, k)\f$ corresond to the\f$ (x, y, z)\f$ coordinates of the
  cell we wish to access. Note: this includes ghost cells, so in practice these values
  can range from \f$ 0 \le (i,j,k) < (nx,ny,nz)+2 Ng \f$ or \f$0 \le (i,j,k) < (Nx,Ny,Nz) \f$
  assuming we are working in 3D. <br>
  For 2D simulations, k=0 at all times, and similarly j=0 for 1D simulations.
*/
class Data
{
  public:
    int
    //@{
    nx, ny, nz;            //!< Number of physical cells in specified direction
    //@}
    double
    //@{
    xmin, xmax,
    ymin, ymax,            //!< Positional limits of domain in specified direction
    zmin, zmax,
    //@}
    endTime,               //!< End time of simulation
    cfl;                   //!< Courant factor
    int Ng;                //!< Number of ghost cells
    double
    gamma,                 //!< Adiabatic index
    sigma;                 //!< Resistivity
    int
    memSet,                //!< Indicator that memory has been allocated for state vectors
    bcsSet,                //!< Indicator that boundary conditions have been created (before this information about the domain decomposition used in MPI version will not be correct).
    //@{
    Ncons, Nprims, Naux;   //!< Number of specified variables
    //@}
    double
    cp,                    //!< Constant divergence cleaning term
    //@{
    mu1, mu2;              //!< Charge mass ratio of specified fluid species, q/m (for two fluid model)
    //@}
    int
    frameSkip,             //!< Number of timesteps per file output
    tpb,                   //!< Threads per block for GPU code
    bpg;                   //!< Blocks per grid for GPU code
    double
    //@{
    *cons, *prims, *aux,
    *f, *fnet,             //!< Pointer to specified work array
    *source,
    //@}
    //@{
    *x, *y, *z;            //!< Specified coordinate location of the center of the compute cells (incl ghost cells)
    //@}
    double
    //@{
    alphaX, alphaY, alphaZ,//!< Max wave speed in specified direction. As we are evolving EM fields, this is always the speed of light.
    //@}
    t=-1,                     //!< Current time
    dt,                    //!< Width of current timestep
    //@{
    dx, dy, dz;            //!< Witdth of specified spatial step
    //@}
    int
    iters,                 //!< Number of interations that have been completed
    //@{
    Nx, Ny, Nz;            //!< Total number of compute cells in domain in the specified direction
    //@}
    std::vector<std::string>
    //@{
    consLabels,
    primsLabels,           //!< Vector of labels for the specified variables
    auxLabels;
    //@}
    int
    dims,                  //!< Number of dimensions of simulation
    //@{
    is, js, ks,
    ie, je, ke;            //!< Cell IDs for interior grid points
    //@}
    std::vector<double>
    optionalSimArgs;          //!< Array of optional arguments that depend on the simulation being run
    std::vector<std::string>
    optionalSimArgNames;     //!< Names of optionalSimArgs array elements
    int
    nOptionalSimArgs=0;      //!< Number of elements to include in optionalSimArgs array
    int
    GPUcount;              //!< Number of NVIDIA devices detected
    cudaDeviceProp
    prop;                  //!< Properties of NVIDIA device (assuming all are same)
    int
    Ncells,                //!< Total number of computational cells in domain
    Nstreams;              //!< Number of GPU streams for the fully paralellisable functions (time int, c2p, etc)

    //! Element ID function
    /*!
    @par
        To access the 2nd conserved variable at (x, y) = (12, 4) for example,
      we call elem=data.id(2, 12, 4) and use this in d.cons[elem].

      @param var the variable number, eg cons[0] is energy, cons[1] is x-momentum etc.
      @param i cell number in the x-direction
      @param j cell number in the y-direction
      @param k cell number in the z-direction
    */
    int id(int var, int i, int j, int k) {
      return var * this->Nx * this->Ny * this->Nz + i * this->Ny * this->Nz + j * this->Nz + k;
    }

    //! Initialiser
    /*!  
        @par
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object. Called automatically from constructors after setting object vars.
      This is separated from the constructor to avoid duplicated code between the two available
      constructors for Data.
     */
     void initData(PlatformEnv *env, int nOptionalSimArgs=0, std::vector<double> optionalSimArgs=std::vector<double>(), std::vector<std::string> optionalSimArgNames=std::vector<std::string>());


    //! Constructor   -- all vars specified by comma separated list
    /*!
      @par
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object.
      @param nx number of physical cells in x-direction
      @param ny number of physical cells in y-direction
      @param nz number of physical cells in z-direction
      @param xmin minimum value of x domain
      @param xmax maximum value of x domain
      @param ymin minimum value of y domain
      @param ymax maximum value of y domain
      @param zmin minimum value of z domain
      @param zmax maximum value of z domain
      @param endTime desired end time of the simulation
      @param cfl courant factor
      @param Ng number of ghost cells in each direction
      @param gamma adiabatic index
      @param sigma value of conductivity
      @param cp time scale for divergence cleaning. cp = 1 / kappa
      @param mu1 charge mass ratio of species 1
      @param mu2 charge mass ratio of species 2
    */
    Data(int nx, int ny, int nz,
         double xmin, double xmax,
         double ymin, double ymax,
         double zmin, double zmax,
         double endTime, PlatformEnv *env,
         double cfl=0.5, int Ng=4,
         double gamma=5.0/3.0, double sigma=1e3,
         double cp=0.1,
         double mu1=-1.0e4, double mu2=1.0e4,
         int frameskip=10);

    //! Constructor
    /*!
      @par
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object.
      @param args simulation arguments such as cfl, sigma etc, as read from checkpoint restart file
      @param mu1 charge mass ratio of species 1
      @param mu2 charge mass ratio of species 2
    */
    Data(DataArgsBase args, PlatformEnv *env);

    ~Data() {};

   
};

#endif
