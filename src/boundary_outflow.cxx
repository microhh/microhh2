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
    void compute_outflow_2nd(
            TF* const restrict a,
            const int iend,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii = 1;

        // Set the ghost cells using extrapolation.
        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
            {
                const int ijk = (iend-1) + j*icells + k*ijcells;
                a[ijk+ii] = a[ijk];
            }
    }

    template<typename TF>
    void compute_inflow_2nd(
            TF* const restrict a, const TF value,
            const int istart,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int ii = 1;

        // Set the ghost cells using extrapolation.
        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
            {
                const int ijk = istart + j*icells + k*ijcells;
                a[ijk-ii] = value - a[ijk];
            }
    }

    template<typename TF>
    void compute_outflow_4th(
            TF* const restrict a,
            const int iend,
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
                const int ijk = (iend-1) + j*icells + k*ijcells;
                a[ijk+ii1] = TF(15./8.)*a[ijk] - TF(10./8.)*a[ijk-ii1] + TF(3./8.)*a[ijk-ii2];
                a[ijk+ii2] = TF(15./8.)*a[ijk] - TF(10./8.)*a[ijk-ii1] + TF(3./8.)*a[ijk-ii2];
                a[ijk+ii3] = TF(15./8.)*a[ijk] - TF(10./8.)*a[ijk-ii1] + TF(3./8.)*a[ijk-ii2];
            }
    }

    template<typename TF>
    void compute_inflow_4th(
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
                a[ijk-ii1] = value + TF( 1./8.)*a[ijk] - TF( 1./4.)*a[ijk+ii1] + TF( 1./8.)*a[ijk+ii2];
                a[ijk-ii2] = value + TF( 9./8.)*a[ijk] - TF( 9./4.)*a[ijk+ii1] + TF( 9./8.)*a[ijk+ii2];
                a[ijk-ii3] = value + TF(25./8.)*a[ijk] - TF(25./4.)*a[ijk+ii1] + TF(25./8.)*a[ijk+ii2];
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
    {
        if (grid.get_spatial_order() == Grid_order::Second)
            compute_outflow_2nd(
                    data,
                    gd.iend,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);

        else if (grid.get_spatial_order() == Grid_order::Fourth)
            compute_outflow_4th(
                    data,
                    gd.iend,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);
    }

    // Inflow
    if (md.mpicoordx == 0)
    {
        if (grid.get_spatial_order() == Grid_order::Second)
            compute_inflow_2nd(
                    data, TF(0.),
                    gd.istart,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);

        else if (grid.get_spatial_order() == Grid_order::Fourth)
            compute_inflow_4th(
                    data, TF(0.),
                    gd.istart,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);
    }
}

template class Boundary_outflow<double>;
template class Boundary_outflow<float>;
