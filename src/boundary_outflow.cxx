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
            TF* const restrict data, TF* const restrict tmp)
    {
    }

    template<typename TF>
    void compute_inflow(
            TF* const restrict data, TF* const restrict tmp, const TF value)
    {
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
void Boundary_outflow<TF>::exec(
        TF* const restrict data, TF* const restrict tmp)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Outflow
    if (md.mpicoordx == md.npx-1)
        compute_outflow(data, tmp);

    if (md.mpicoordx == 0)
        compute_inflow(data, tmp, TF(0.));
}

template class Boundary_outflow<double>;
template class Boundary_outflow<float>;
