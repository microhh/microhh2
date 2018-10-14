/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "master.h"
#include "grid.h"
#include "boundary_outflow.h"

namespace
{
    template<typename TF>
    void compute_outflow(
            TF* const restrict a,
            const int iend,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii1 = 1;
        const int ii2 = 2;

        // Set the ghost cells using extrapolation.
        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
            {
                const int ijk = (iend-1) + j*icells + k*ijcells;
                a[ijk+ii1] = TF(3.)*a[ijk] - TF(3.)*a[ijk-ii1] + a[ijk-ii2];
                a[ijk+ii2] = TF(6.)*a[ijk] - TF(8.)*a[ijk-ii1] + TF(3.)*a[ijk-ii2];
            }
    }

    template<typename TF>
    void compute_inflow(
            TF* const restrict a, const TF value,
            const int istart,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        // Set the ghost cells using extrapolation.
        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
            {
                const int ijk = istart + j*icells + k*ijcells;
                a[ijk-ii1] = value;
                a[ijk-ii2] = value;
                a[ijk-ii3] = value;
            }
    }
}

template<typename TF>
Boundary_outflow<TF>::Boundary_outflow(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin)
{
}

template<typename TF>
Boundary_outflow<TF>::~Boundary_outflow()
{
}

template<typename TF>
void Boundary_outflow<TF>::exec(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Outflow
    if (md.mpicoordx == md.npx-1)
        compute_outflow(
                data,
                gd.iend,
                gd.icells, gd.jcells, gd.kcells,
                gd.ijcells);

    if (md.mpicoordx == 0)
        compute_inflow(
                data, TF(0.),
                gd.istart,
                gd.icells, gd.jcells, gd.kcells,
                gd.ijcells);
}

template class Boundary_outflow<double>;
template class Boundary_outflow<float>;
