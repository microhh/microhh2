/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "defines.h"
#include "timeloop.h"
#include "timedep.h"
#include "finite_difference.h"
#include "netcdf_interface.h"
#include <iostream>
#include <vector>
#include <chrono>

#include "boundary_cyclic.h"

// Boundary schemes.
#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_surface_bulk.h"
// #include "boundary_surface_patch.h"
// #include "boundary_patch.h"


namespace
{
    template<typename TF>
    void set_bc(TF* const restrict a, TF* const restrict agrad, TF* const restrict aflux,
                const Boundary_type sw, const TF aval, const TF visc, const TF offset,
                const int icells, const int jcells)
    {
        const int jj = icells;

        if (sw == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    a[ij] = aval - offset;
                }
        }
        else if (sw == Boundary_type::Neumann_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    agrad[ij] = aval;
                    aflux[ij] = -aval*visc;
                }
        }
        else if (sw == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    aflux[ij] = aval;
                    agrad[ij] = -aval/visc;
                }
        }
    }
}

template<typename TF>
Boundary<TF>::Boundary(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin),
    grid(gridin),
    fields(fieldsin),
    boundary_cyclic(master, grid),
    field3d_io(master, grid)
{
    swboundary = "default";
}

template<typename TF>
Boundary<TF>::~Boundary()
{
    // empty the map
    // CvH: is this necessary?
    sbc.clear();

    // clean up time dependent data
    // for (auto& i : timedepdata)
    //     delete[] i.second;
}

template<typename TF>
std::string Boundary<TF>::get_switch()
{
    return swboundary;
}

template<typename TF>
void Boundary<TF>::process_bcs(Input& input)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    std::string swbot = input.get_item<std::string>("boundary", "mbcbot", "");
    std::string swtop = input.get_item<std::string>("boundary", "mbctop", "");

    ubot = input.get_item<TF>("boundary", "ubot", "", 0.);
    utop = input.get_item<TF>("boundary", "utop", "", 0.);
    vbot = input.get_item<TF>("boundary", "vbot", "", 0.);
    vtop = input.get_item<TF>("boundary", "vtop", "", 0.);

    // set the bottom bc
    if (swbot == "noslip")
        mbcbot = Boundary_type::Dirichlet_type;
    else if (swbot == "freeslip" )
        mbcbot = Boundary_type::Neumann_type;
    else if (swbot == "neumann" )
        mbcbot = Boundary_type::Neumann_type;
    else if (swbot == "ustar")
        mbcbot = Boundary_type::Ustar_type;
    else
    {
        std::string msg = swbot + " is an illegal value for mbcbot";
        throw std::runtime_error(msg);
    }

    // set the top bc
    if (swtop == "noslip")
        mbctop = Boundary_type::Dirichlet_type;
    else if (swtop == "freeslip")
        mbctop = Boundary_type::Neumann_type;
    else if (swtop == "neumann")
        mbctop = Boundary_type::Neumann_type;
    else if (swtop == "ustar")
        mbctop = Boundary_type::Ustar_type;
    else
    {
        std::string msg = swtop + " is an illegal value for mbctop";
        throw std::runtime_error(msg);
    }

    // read the boundaries per field
    for (auto& it : fields.sp)
    {
        sbc.emplace(it.first, Field3dBc<TF>());
        swbot = input.get_item<std::string>("boundary", "sbcbot", it.first);
        swtop = input.get_item<std::string>("boundary", "sbctop", it.first);
        sbc.at(it.first).bot = input.get_item<TF>("boundary", "sbot", it.first);
        sbc.at(it.first).top = input.get_item<TF>("boundary", "stop", it.first);

        // set the bottom bc
        if (swbot == "dirichlet")
            sbc.at(it.first).bcbot = Boundary_type::Dirichlet_type;
        else if (swbot == "neumann")
            sbc.at(it.first).bcbot = Boundary_type::Neumann_type;
        else if (swbot == "flux")
            sbc.at(it.first).bcbot = Boundary_type::Flux_type;
        else
        {
            std::string msg = swbot + " is an illegal value for sbcbot";
            throw std::runtime_error(msg);
        }

        // set the top bc
        if (swtop == "dirichlet")
            sbc.at(it.first).bctop = Boundary_type::Dirichlet_type;
        else if (swtop == "neumann")
            sbc.at(it.first).bctop = Boundary_type::Neumann_type;
        else if (swtop == "flux")
            sbc.at(it.first).bctop = Boundary_type::Flux_type;
        else
        {
            std::string msg = swbot + " is an illegal value for sbctop";
            throw std::runtime_error(msg);
        }
    }

    sbot_2d_list = input.get_list<std::string>("boundary", "sbot_2d_list", "", std::vector<std::string>());

    std::string swopenbc_in = input.get_item<std::string>("boundary", "swopenbc", "", "0");

    if (swopenbc_in == "0"){
        swopenbc = Openbc_type::disabled;
      }
    else if (swopenbc_in == "1")
    {
        swopenbc = Openbc_type::enabled;
        openbc_list = input.get_list<std::string>("boundary", "openbc_list", "", std::vector<std::string>());
        for (auto& it : openbc_list){
            openbc_profs[it] = std::vector<TF>(4*gd.kcells);
          }
    }
}

template<typename TF>
void Boundary<TF>::init(Input& input, Thermo<TF>& thermo)
{
    // Read the boundary information from the ini files, it throws at error.
    process_bcs(input);

    // there is no option (yet) for prescribing ustar without surface model
    if (mbcbot == Boundary_type::Ustar_type || mbctop == Boundary_type::Ustar_type)
    {
        throw std::runtime_error ("Cannot use ustar bc for default boundary");
    }

    // Initialize the boundary cyclic.
    boundary_cyclic.init();

    // Initialize the IO operators.
    field3d_io.init();
}

template<typename TF>
void Boundary<TF>::create(Input& input, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    process_openbc(input_nc);
    process_time_dependent(input, input_nc);
}

template<typename TF>
void Boundary<TF>::process_openbc(Netcdf_handle& input_nc)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    Netcdf_group& group_nc = input_nc.get_group("init");

    if (swopenbc == Openbc_type::enabled)
    {
        // Check whether the fields in the list exist in the prognostic fields.
        for (std::string& it : openbc_list)
        {
            if (!fields.ap.count(it))
            {
                std::string msg = "field " + it + " in [boundary][openbclist] is illegal";
                throw std::runtime_error(msg);
            }
            else
            {
                std::vector<TF> tmp(4*gd.ktot);
                group_nc.get_variable(tmp, it+"_bc", {0,0}, {4,gd.ktot});
                for (int n = 0; n < 4; n++){
                    std::copy(tmp.begin() + n * gd.ktot, tmp.begin() + (n+1) * gd.ktot, openbc_profs[it].begin() + gd.kstart + n * gd.kcells);
                  }
                // Modify subsidence

                // Modify initial field (done in field.cxx)
                // Subtract utrans and vtrans, if necessary
            }
        }
    }
}


template<typename TF>
void Boundary<TF>::process_time_dependent(Input& input, Netcdf_handle& input_nc)
{
    // get the list of time varying variables
    bool swtimedep = input.get_item<bool>("boundary", "swtimedep"  , "", false);
    std::vector<std::string> timedeplist = input.get_list<std::string>("boundary", "timedeplist", "", std::vector<std::string>());

    if (swtimedep)
    {
        if (!sbot_2d_list.empty())
            master.print_warning("Provided 2D sbot fields are potentially overwritten by timedep");

        // Create temporary list to check which entries are used.
        std::vector<std::string> tmplist = timedeplist;

        // See if there is data available for the surface boundary conditions.
        for (auto& it : fields.sp)
        {
            std::string name = it.first+"_sbot";
            if (std::find(timedeplist.begin(), timedeplist.end(), name) != timedeplist.end())
            {
                // Process the time dependent data.
                tdep_bc.emplace(it.first, new Timedep<TF>(master, grid, name, true));
                tdep_bc.at(it.first)->create_timedep(input_nc);

                // Remove the item from the tmplist.
                std::vector<std::string>::iterator ittmp = std::find(tmplist.begin(), tmplist.end(), name);
                if (ittmp != tmplist.end())
                    tmplist.erase(ittmp);
            }
        }

        // Display a warning for the non-supported.
        for (std::vector<std::string>::const_iterator ittmp=tmplist.begin(); ittmp!=tmplist.end(); ++ittmp)
            master.print_warning("%s is not supported (yet) as a time dependent parameter\n", ittmp->c_str());
    }

}
#ifndef USECUDA
template <typename TF>
void Boundary<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const TF no_offset = 0.;

    for (auto& it : tdep_bc)
    {
        it.second->update_time_dependent(sbc.at(it.first).bot, timeloop);
        set_bc<TF>(fields.sp.at(it.first)->fld_bot.data(), fields.sp.at(it.first)->grad_bot.data(), fields.sp.at(it.first)->flux_bot.data(),
                sbc.at(it.first).bcbot, sbc.at(it.first).bot, fields.sp.at(it.first)->visc, no_offset, gd.icells, gd.jcells);
    }
}
#endif

template<typename TF>
void Boundary<TF>::set_values()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    set_bc<TF>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
           mbcbot, ubot, fields.visc, grid.utrans,
           gd.icells, gd.jcells);
    set_bc<TF>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
           mbcbot, vbot, fields.visc, grid.vtrans,
           gd.icells, gd.jcells);

    set_bc<TF>(fields.mp.at("u")->fld_top.data(), fields.mp.at("u")->grad_top.data(), fields.mp.at("u")->flux_top.data(),
           mbctop, utop, fields.visc, grid.utrans,
           gd.icells, gd.jcells);
    set_bc<TF>(fields.mp.at("v")->fld_top.data(), fields.mp.at("v")->grad_top.data(), fields.mp.at("v")->flux_top.data(),
           mbctop, vtop, fields.visc, grid.vtrans,
           gd.icells, gd.jcells);

    const TF no_offset = 0.;

    for (auto& it : fields.sp)
    {
        // Load 2D fields for bottom boundary from disk.
        if (std::find(sbot_2d_list.begin(), sbot_2d_list.end(), it.first) != sbot_2d_list.end())
        {
            std::string filename = it.first + "_bot.0000000";
            master.print_message("Loading \"%s\" ... ", filename.c_str());

            auto tmp = fields.get_tmp();
            TF* fld_2d_ptr = nullptr;
            if (sbc.at(it.first).bcbot == Boundary_type::Dirichlet_type)
                fld_2d_ptr = it.second->fld_bot.data();
            else if (sbc.at(it.first).bcbot == Boundary_type::Neumann_type)
                fld_2d_ptr = it.second->grad_bot.data();
            else if (sbc.at(it.first).bcbot == Boundary_type::Flux_type)
                fld_2d_ptr = it.second->flux_bot.data();

            if (field3d_io.load_xy_slice(fld_2d_ptr, tmp->fld.data(), filename.c_str()))
            {
                master.print_message("FAILED\n");
                throw std::runtime_error("Error loading 2D field of bottom boundary");
            }
            fields.release_tmp(tmp);
        }
        else
        {
            set_bc<TF>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
                   sbc.at(it.first).bcbot, sbc.at(it.first).bot, it.second->visc, no_offset,
                   gd.icells, gd.jcells);
            set_bc<TF>(it.second->fld_top.data(), it.second->grad_top.data(), it.second->flux_top.data(),
                   sbc.at(it.first).bctop, sbc.at(it.first).top, it.second->visc, no_offset,
                   gd.icells, gd.jcells);
        }
    }
}

namespace
{
    template<typename TF>
    void calc_ghost_cells_bot_2nd(TF* const restrict a, const TF* const restrict dzh, Boundary_type boundary_type,
                                  TF* const restrict abot, TF* const restrict agradbot,
                                  const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    a[ijk-kk] = TF(2.)*abot[ij] - a[ijk];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
                }
        }
    }

    template<typename TF>
    void calc_ghost_cells_top_2nd(TF* const restrict a, const TF* const restrict dzh, Boundary_type boundary_type,
                                  TF* const restrict atop, TF* const restrict agradtop,
                                  const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    a[ijk+kk] = TF(2.)*atop[ij] - a[ijk];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
                }
        }
    }

    template<typename TF>
    void calc_ghost_cells_bot_4th(TF* const restrict a, const TF* const restrict z, Boundary_type boundary_type,
                                  TF* const restrict abot, TF* restrict agradbot,
                                  const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    a[ijk-kk1] = TF(8./3.)*abot[ij] - TF(2.)*a[ijk] + TF(1./3.)*a[ijk+kk1];
                    a[ijk-kk2] = TF(8.)*abot[ij] - TF(9.)*a[ijk] + TF(2.)*a[ijk+kk1];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            using Finite_difference::O4::grad4;

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    a[ijk-kk1] = TF(-1.)*grad4(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
                    a[ijk-kk2] = TF(-3.)*grad4(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
                }
        }
    }

    template<typename TF>
    void calc_ghost_cells_top_4th(TF* const restrict a, const TF* const restrict z, Boundary_type boundary_type,
                                  TF* const restrict atop, TF* const restrict agradtop,
                                  const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk1;
                    a[ijk+kk1] = TF(8./3.)*atop[ij] - TF(2.)*a[ijk] + TF(1./3.)*a[ijk-kk1];
                    a[ijk+kk2] = TF(8.)*atop[ij] - TF(9.)*a[ijk] + TF(2.)*a[ijk-kk1];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            using Finite_difference::O4::grad4;

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk1;
                    a[ijk+kk1] = TF(1.)*grad4(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
                    a[ijk+kk2] = TF(3.)*grad4(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
                }
        }
    }

    // BOUNDARY CONDITIONS FOR THE VERTICAL VELOCITY (NO PENETRATION)
    template<typename TF>
    void calc_ghost_cells_botw_cons_4th(TF* const restrict w,
            const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kstart*kk1;
                w[ijk-kk1] = -w[ijk+kk1];
                w[ijk-kk2] = -w[ijk+kk2];
            }
    }

    template<typename TF>
    void calc_ghost_cells_topw_cons_4th(TF* const restrict w,
            const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kend*kk1;
                w[ijk+kk1] = -w[ijk-kk1];
                w[ijk+kk2] = -w[ijk-kk2];
            }
    }

    template<typename TF>
    void calc_ghost_cells_botw_4th(TF* restrict w,
            const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kstart*kk1;
                w[ijk-kk1] = TF(-6.)*w[ijk+kk1] + TF(4.)*w[ijk+kk2] - w[ijk+kk3];
            }
    }

    template<typename TF>
    void calc_ghost_cells_topw_4th(TF* restrict w,
            const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kend*kk1;
                w[ijk+kk1] = TF(-6.)*w[ijk-kk1] + TF(4.)*w[ijk-kk2] - w[ijk-kk3];
            }
    }
    template<typename TF>
    void calc_openbc(TF* restrict data, TF* corners,const int xsize, const int ysize, const int dx, const int dy,
            const int igc, const int jgc, const int itot, const int jtot, const int iend, const int jend, const int kend, const int istart,const int jstart,const int kstart, const int icells, const int jcells, const int kcells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // West boundary
        for (int k=kstart; k<kend; ++k)
        {
          int n=jend-jstart;
          #pragma ivdep
          /*We impose the gradient on the qhost cells and the first grid cells, which is why
          we loop until i<igc+1 (to include first grid cell)*/
          for(int i=0;i<igc+1;i++){
            /*The first part of the loop consists of using linear regression to
            determine the slope and intercept of the data on this side */
            TF avg_x=0;
            TF sum_squared_xdiff=0;
            TF inprod_xy=0;
            TF avg_y=0;
            /*This loop computes the average of the data, needed to use linear regression*/
            for(int j=jstart;j<jend;j++){
              avg_x+=j;
              const int ijk0 = i + (j)*jj + k*kk;
              avg_y+=data[ijk0];
            }
            avg_x/=n; /*Computing the average of the x values, i.e. the cell numbers*/
            avg_y/=n; /*Computing the average of the y values, i.e. the data points, e.g. Thl or Qt*/
            for(int j=jstart;j<jend;j++){
              const int ijk0 = i + (j)*jj + k*kk;
              inprod_xy+=(j-avg_x)*(data[ijk0]-avg_y); /*Computing the covariance of the x and y values */
              sum_squared_xdiff+=(j-avg_x)*(j-avg_x); /*Computing the variance of the x values, needed for LR*/
            }
            TF slopeLRW=inprod_xy / sum_squared_xdiff; /*Formula for slope by using linear regression*/
            TF interceptLRW=avg_y - slopeLRW * avg_x; /*Formula for intercept by using linear regression*/

            TF DeltaSlopeW = (corners[k + 2 *kcells] - corners[k + 0 *kcells]) / (jtot-1)
                        - slopeLRW; /*Slope difference between imposed gradient and linear regression slope of actual data*/
            TF DeltaInterceptW = (corners[k + 0 *kcells]-(corners[k + 2 *kcells] - corners[k + 0 *kcells])*jstart / (jtot-1))
                        - interceptLRW;/*Intercept difference between imposed gradient and linear regression intercept of actual data*/
            for(int j=jstart;j<jend;j++){
              const int ijk0=i + j*jj + k*kk;
              data[ijk0] += DeltaSlopeW*j + DeltaInterceptW; /*Adding the slope and intercept difference to match the imposed gradient*/
            }
          }
        }
        /*Updating the other boundaries happens in the exact same way by using linear regression. Since it is completely
        analogous to the west boundary, no comments are added to the following section of the code.*/
        //East boundary

        for (int k=kstart; k<kend; ++k)
        {
          int n=jend-jstart;
          #pragma ivdep
          for(int i=0;i<igc+1;i++){

            TF avg_x=0;
            TF sum_squared_xdiff=0;
            TF inprod_xy=0;
            TF avg_y=0;

            for(int j=jstart;j<jend;j++){
              avg_x+=j;
              const int ijk0 = i + iend - 1 + (j)*jj + k*kk;
              avg_y+=data[ijk0];
            }
            avg_x/=n;
            avg_y/=n;
            for(int j=jstart;j<jend;j++){
              const int ijk0 = i + iend - 1 + (j)*jj + k*kk;
              inprod_xy+=(j-avg_x)*(data[ijk0]-avg_y);
              sum_squared_xdiff+=(j-avg_x)*(j-avg_x);
            }
            TF slopeLRE=inprod_xy / sum_squared_xdiff;
            TF interceptLRE=avg_y - slopeLRE * avg_x;
            TF DeltaSlopeE = (corners[k + 3 *kcells] - corners[k + 1 *kcells]) / (jtot-1)
                        - slopeLRE;
            TF DeltaInterceptE = (corners[k + 1 *kcells]-(corners[k + 3 *kcells] - corners[k + 1 *kcells])*jstart / (jtot-1))
                        - interceptLRE;
            for(int j=jstart;j<jend;j++){
              const int ijk0= i + iend - 1 + j*jj + k*kk;
              data[ijk0] += DeltaSlopeE*j + DeltaInterceptE;
            }
          }
        }

        // North boundary
        for (int k=kstart; k<kend; ++k)
        {
          int n=iend-istart;
          #pragma ivdep
          for(int j=0;j<jgc+1;j++){
            TF avg_x=0;
            TF sum_squared_xdiff=0;
            TF inprod_xy=0;
            TF avg_y=0;

            for(int i=istart;i<iend;i++){
              avg_x+=i;
              const int ijk0 = i + (j+jend-1)*jj + k*kk;
              avg_y+=data[ijk0];
            }
            avg_x/=n;
            avg_y/=n;
            for(int i=istart;i<iend;i++){
              const int ijk0 = i + (j+jend-1)*jj + k*kk;
              inprod_xy+=(i-avg_x)*(data[ijk0]-avg_y);
              sum_squared_xdiff+=(i-avg_x)*(i-avg_x);
            }
            TF slopeLRN=inprod_xy / sum_squared_xdiff;
            TF interceptLRN=avg_y - slopeLRN * avg_x;

            TF DeltaSlopeN = (corners[k + 3 *kcells] - corners[k + 2 *kcells]) / (itot-1)
                          - slopeLRN;
            TF DeltaInterceptN = (corners[k + 2 *kcells]-(corners[k + 3 *kcells] - corners[k + 2 *kcells])*istart / (itot-1))
                          - interceptLRN;

            if(j==jend-1){
              for(int i=istart+1;i<iend-1;i++){/*To not add to the first gridpoints that were
                already updated in the east-west boundary part, this if statement is included. Note:
                the j=jend-1 is exactly the northern first gridpoints.*/
                const int ijk0=i + (j+jend-1)*jj + k*kk;
                data[ijk0] += DeltaSlopeN*i + DeltaInterceptN;
              }
            }
            else{
              for(int i=istart;i<iend;i++){/*For the ghost cells, the endpoints are included as usual.*/
                const int ijk0=i + (j+jend-1)*jj + k*kk;
                data[ijk0] += DeltaSlopeN*i + DeltaInterceptN;
              }
            }
          }
        }

        //South boundary
        for (int k=kstart; k<kend; ++k)
        {
          int n=iend-istart;
          #pragma ivdep
          for(int j=0;j<jgc+1;j++){
            TF avg_x=0;
            TF sum_squared_xdiff=0;
            TF inprod_xy=0;
            TF avg_y=0;

            for(int i=istart;i<iend;i++){
              avg_x+=i;
              const int ijk0 = i + (j)*jj + k*kk;
              avg_y+=data[ijk0];
            }
            avg_x/=n;
            avg_y/=n;
            for(int i=istart;i<iend;i++){
              const int ijk0 = i + (j)*jj + k*kk;
              inprod_xy+=(i-avg_x)*(data[ijk0]-avg_y);
              sum_squared_xdiff+=(i-avg_x)*(i-avg_x);
            }
            TF slopeLRS=inprod_xy / sum_squared_xdiff;
            TF interceptLRS=avg_y - slopeLRS * avg_x;

            TF DeltaSlopeS = (corners[k + 1 *kcells] - corners[k + 0 *kcells]) / (itot-1)
                          - slopeLRS;
            TF DeltaInterceptS = (corners[k + 0 *kcells]-(corners[k + 1 *kcells] - corners[k + 0 *kcells])*istart / (itot-1))
                          - interceptLRS;
            if(j==jstart){
              for(int i=istart+1;i<iend-1;i++){
                const int ijk0=i + (j)*jj + k*kk;
                data[ijk0] += DeltaSlopeS*i + DeltaInterceptS;
              }
            }
            else{
              for(int i=istart;i<iend;i++){
                const int ijk0=i + (j)*jj + k*kk;
                data[ijk0] += DeltaSlopeS*i + DeltaInterceptS;
              }
            }
          }
        }
    }
}

#ifndef USECUDA
template<typename TF>
void Boundary<TF>::exec(Thermo<TF>& thermo)
{
      const Grid_data<TF>& gd = grid.get_grid_data();
      /*We now want to apply the cyclic boundary conditions to the rest of the fields. If openbc is disabled,
      then no fields were updated previously, so update all the fields. On the other hand, if openbc is enabled,
      then fields in openbc_list were already updated, so only update the ones that are not in openbc_list.*/
      for (auto& it : fields.sp){
          if (swopenbc == Openbc_type::disabled||!(std::find(openbc_list.begin(), openbc_list.end(), it.first) != openbc_list.end())){
            boundary_cyclic.exec(it.second->fld.data());
            }
        }
      for (auto& it : fields.sp){
        /*Include if statement to ensure that the cyclic bc happens before open_bc, only when openbc is enabled
        (first part of if statement) and only for the fields in open_bc_list (second part of if statement checks
        if it is in openbc_list)*/
          if (swopenbc == Openbc_type::enabled && (std::find(openbc_list.begin(), openbc_list.end(), it.first) != openbc_list.end())){
            boundary_cyclic.exec(it.second->fld.data());
          }
        }
      /*If openbc is enabled, then perform the code to compute the openbc only for the items in openbc_list*/
      if (swopenbc == Openbc_type::enabled){
          for (auto& it : openbc_list)
          {
              calc_openbc(fields.sp.at("thl")->fld.data(), openbc_profs.at(it).data(),gd.xsize,gd.ysize,gd.dx,gd.dy,
                  gd.igc, gd.jgc, gd.itot, gd.jtot, gd.iend, gd.jend, gd.kend, gd.istart, gd.jstart,gd.kstart, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
          }
        }
      /*Execute the cyclic boundary part for the velocities*/
      boundary_cyclic.exec(fields.mp.at("u")->fld.data());
      boundary_cyclic.exec(fields.mp.at("v")->fld.data());
      boundary_cyclic.exec(fields.mp.at("w")->fld.data());


    // Update the boundary values.
    update_bcs(thermo);

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        calc_ghost_cells_bot_2nd<TF>(fields.mp.at("u")->fld.data(), gd.dzh.data(), mbcbot,
                fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_2nd<TF>(fields.mp.at("u")->fld.data(), gd.dzh.data(), mbctop,
                fields.mp.at("u")->fld_top.data(), fields.mp.at("u")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        calc_ghost_cells_bot_2nd<TF>(fields.mp.at("v")->fld.data(), gd.dzh.data(), mbcbot,
                fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_2nd<TF>(fields.mp.at("v")->fld.data(), gd.dzh.data(), mbctop,
                fields.mp.at("v")->fld_top.data(), fields.mp.at("v")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_2nd<TF>(it.second->fld.data(), gd.dzh.data(),
                    sbc.at(it.first).bcbot, it.second->fld_bot.data(), it.second->grad_bot.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_top_2nd<TF>(it.second->fld.data(), gd.dzh.data(),
                    sbc.at(it.first).bctop, it.second->fld_top.data(), it.second->grad_top.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        calc_ghost_cells_bot_4th<TF>(fields.mp.at("u")->fld.data(), gd.z.data(), mbcbot,
                fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_4th<TF>(fields.mp.at("u")->fld.data(), gd.z.data(), mbctop,
                fields.mp.at("u")->fld_top.data(), fields.mp.at("u")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        calc_ghost_cells_bot_4th<TF>(fields.mp.at("v")->fld.data(), gd.z.data(), mbcbot,
                fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_4th<TF>(fields.mp.at("v")->fld.data(), gd.z.data(), mbctop,
                fields.mp.at("v")->fld_top.data(), fields.mp.at("v")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        calc_ghost_cells_botw_4th<TF>(fields.mp.at("w")->fld.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_topw_4th<TF>(fields.mp.at("w")->fld.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_4th<TF>(it.second->fld.data(), gd.z.data(), sbc.at(it.first).bcbot,
                    it.second->fld_bot.data(), it.second->grad_bot.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_top_4th<TF>(it.second->fld.data(), gd.z.data(), sbc.at(it.first).bctop,
                    it.second->fld_top.data(), it.second->grad_top.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
    }

    // Update the boundary fields that are a slave of the boundary condition.
    update_slave_bcs();
}

template<typename TF>
void Boundary<TF>::set_ghost_cells_w(const Boundary_w_type boundary_w_type)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        if (boundary_w_type == Boundary_w_type::Normal_type)
        {
            calc_ghost_cells_botw_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_topw_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
        else if (boundary_w_type == Boundary_w_type::Conservation_type)
        {
            calc_ghost_cells_botw_cons_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_topw_cons_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
    }
}
#endif

/*
template<typename TF>
void Boundary<TF>::exec_cross()
{
}
*/

template<typename TF>
void Boundary<TF>::exec_stats(Stats<TF>&)
{
}

// Computational kernel for boundary calculation.
namespace
{
    template<typename TF, int spatial_order>
    void calc_slave_bc_bot(TF* const restrict abot, TF* const restrict agradbot, TF* const restrict afluxbot,
                           const TF* const restrict a, const TF* const restrict dzhi,
                           const Boundary_type boundary_type, const TF visc,
                           const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        using namespace Finite_difference;

        // Variable dzhi in this case is dzhi for 2nd order and dzhi4 for 4th order.
        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    if (spatial_order == 2)
                        agradbot[ij] = O2::grad2(a[ijk-kk1], a[ijk]) * dzhi[kstart];
                    else if (spatial_order == 4)
                        agradbot[ij] = O4::grad4(a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1]) * dzhi[kstart];
                    afluxbot[ij] = -visc*agradbot[ij];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    if (spatial_order == 2)
                        abot[ij] = O2::interp2(a[ijk-kk1], a[ijk]);
                    else if (spatial_order == 4)
                        abot[ij] = O4::interp4c(a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1]);
                }
        }
    }
}

template<typename TF>
void Boundary<TF>::update_slave_bcs()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        calc_slave_bc_bot<TF,2>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
                                fields.mp.at("u")->fld.data(), gd.dzhi.data(),
                                mbcbot, fields.mp.at("u")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        calc_slave_bc_bot<TF,2>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
                                fields.mp.at("v")->fld.data(), gd.dzhi.data(),
                                mbcbot, fields.mp.at("v")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
            calc_slave_bc_bot<TF,2>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
                                    it.second->fld.data(), gd.dzhi.data(),
                                    sbc.at(it.first).bcbot, it.second->visc,
                                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        calc_slave_bc_bot<TF,4>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
                                fields.mp.at("u")->fld.data(), gd.dzhi4.data(),
                                mbcbot, fields.mp.at("u")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        calc_slave_bc_bot<TF,4>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
                                fields.mp.at("v")->fld.data(), gd.dzhi4.data(),
                                mbcbot, fields.mp.at("v")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
            calc_slave_bc_bot<TF,4>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
                                    it.second->fld.data(), gd.dzhi4.data(),
                                    sbc.at(it.first).bcbot, it.second->visc,
                                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
    }
}

template<typename TF>
void Boundary<TF>::update_bcs(Thermo<TF>& thermo)
{
}

template<typename TF>
std::shared_ptr<Boundary<TF>> Boundary<TF>::factory(Master& master, Grid<TF>& grid, Fields<TF>& fields, Input& input)
{
    std::string swboundary;
    swboundary = input.get_item<std::string>("boundary", "swboundary", "", "default");

    // else if (swboundary == "surface_patch")
    //     return new Boundary_surface_patch(modelin, inputin);
    // else if (swboundary == "patch")
    //     return new Boundary_patch(modelin, inputin);
    // else if (swboundary == "default")
    if (swboundary == "default")
        return std::make_shared<Boundary<TF>>(master, grid, fields, input);
    else if (swboundary == "surface")
        return std::make_shared<Boundary_surface<TF>>(master, grid, fields, input);
    else if (swboundary == "surface_bulk")
        return std::make_shared<Boundary_surface_bulk<TF>>(master, grid, fields, input);
    else
    {
        std::string msg = swboundary + " is an illegal value for swboundary";
        throw std::runtime_error(msg);
    }
}
/*
template<typename TF>
void Boundary<TF>::get_mask(Field3d* field, Field3d* fieldh, Mask* m)
{
    // Set surface mask
    for (int i=0; i<grid.ijcells; ++i)
        fieldh->fld_bot[i] = 1;

    // Set atmospheric mask
    for (int i=0; i<grid.ncells; ++i)
    {
        field ->data[i] = 1;
        fieldh->data[i] = 1;
    }
}

template<typename TF>
void Boundary<TF>::get_surface_mask(Field3d* field)
{
    // Set surface mask
    for (int i=0; i<grid.ijcells; ++i)
        field->fld_bot[i] = 1;
}

*/
/* template<typename TF>
void Boundary<TF>::prepare_device()
{
}

template<typename TF>
void Boundary<TF>::forward_device()
{
}

template<typename TF>
void Boundary<TF>::backward_device()
{
}
 */
template class Boundary<double>;
template class Boundary<float>;
