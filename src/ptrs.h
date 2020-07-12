

#pragma once

#include "mpi.h"

// LAMMPS include files
#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "library.h"
#include "memory.h"

#include "geo.h"

namespace GEO_NS
{
    class Ptrs
    {
    public:
        Ptrs(GEO *ptr) :
            geo(ptr),
            mem(ptr->mem),
            //poptimer(ptr->poptimer),
            //popinput(ptr->popinput),
            //config(ptr->config),
            lmp(ptr->lmp)
            {}

        virtual ~Ptrs() {}

    protected:
        GEO *geo;
        Mem *&mem;
        //PopTimer *&poptimer;
        //PopInput *&popinput;
        //Config *&config;
        LAMMPS_NS::LAMMPS *&lmp;
    };
}

