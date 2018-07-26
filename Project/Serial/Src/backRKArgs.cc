#include "backRKArgs.h"

//! Additional arguments parameterized constructor
BackRKArguments::BackRKArguments(Data * data) : data(data)
{
  // Small arrays, no need to malloc
  constar    = new double[data->Ncons ];
  primstar   = new double[data->Nprims];
  auxstar    = new double[data->Naux  ];
  sourcestar = new double[data->Ncons ];
  allocd = 1;
}

BackRKArguments::~BackRKArguments()
{
  delete [] constar;
  delete [] primstar;
  delete [] auxstar;
  delete [] sourcestar;
}

//! Overload assignment operator
BackRKArguments& BackRKArguments::operator=(const BackRKArguments &args)
{
  // Set simulation data
  data = args.data;

  // If no memory has been allocated, allocate
  if (!allocd) {
    constar    = new double[data->Ncons ];
    primstar   = new double[data->Nprims];
    auxstar    = new double[data->Naux  ];
    sourcestar = new double[data->Ncons ];
  }

  // Copy accross data
  for (int i(0); i < data->Ncons ; i++) constar[i]    = args.constar[i];
  for (int i(0); i < data->Nprims; i++) primstar[i]   = args.primstar[i];
  for (int i(0); i < data->Naux  ; i++) auxstar[i]    = args.auxstar[i];
  for (int i(0); i < data->Ncons ; i++) sourcestar[i] = args.sourcestar[i];
  return *this;
}
