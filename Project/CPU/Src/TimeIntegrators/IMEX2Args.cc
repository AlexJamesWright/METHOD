#include "IMEX2Args.h"

//! Additional arguments parameterized constructor
IMEX2Arguments::IMEX2Arguments(Data * data) : data(data),
                                              gam(0.2928932188134525),
                                              om2gam(0.4142135623730949)
{
  // // Small arrays, no need to malloc
  cons     = new double[data->Ncons ];
  prims    = new double[data->Nprims];
  aux      = new double[data->Naux  ];
  source   = new double[data->Ncons ];
  source1  = new double[data->Ncons ];
  flux1    = new double[data->Ncons ];
  source2  = new double[data->Ncons ];
  flux2    = new double[data->Ncons ];
  allocd = 1;
}

IMEX2Arguments::~IMEX2Arguments()
{

  delete [] cons;
  delete [] prims;
  delete [] aux;
  delete [] source;
  delete [] source1;
  delete [] flux1;
  delete [] source2;
  delete [] flux2;
  allocd = 0;
}

//! Overload assignment operator
IMEX2Arguments& IMEX2Arguments::operator=(const IMEX2Arguments &args)
{
  // Set simulation data
  data = args.data;

  // If no memory has been allocated, allocate
  if (!allocd) {

    cons     = new double[data->Ncons ];
    prims    = new double[data->Nprims];
    aux      = new double[data->Naux  ];
    source   = new double[data->Ncons ];
    source1  = new double[data->Ncons ];
    flux1    = new double[data->Ncons ];
    source2  = new double[data->Ncons ];
    flux2    = new double[data->Ncons ];
    allocd = 1;
  }

  // Copy accross data
  for (int i(0); i < data->Ncons ; i++) cons [i]   = args.cons [i];
  for (int i(0); i < data->Nprims; i++) prims[i]   = args.prims[i];
  for (int i(0); i < data->Naux  ; i++) aux[i]     = args.aux[i];
  for (int i(0); i < data->Ncons ; i++) source[i]  = args.source[i];
  for (int i(0); i < data->Ncons ; i++) source1[i] = args.source1[i];
  for (int i(0); i < data->Ncons ; i++) flux1[i]   = args.flux1[i];
  for (int i(0); i < data->Ncons ; i++) source2[i] = args.source2[i];
  for (int i(0); i < data->Ncons ; i++) flux2[i]   = args.flux2[i];
  return *this;
}
