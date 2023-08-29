//===========================================================================
//
// File: LookUpCellCentroid.hh
//
// Created: Wed July 26 10:48:00 2023
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2023 Equinor ASA.

  This file is part of The Open Porous Media project  (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <dune/grid/common/mcmgmapper.hh>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>



#include <type_traits>

namespace Dune
{
class CpGrid;
}

namespace Opm
{

class EclipseGrid;

/// LookUpCellCentroid struct - To search cell centroids via element index
///                             Using a specialitation for Dune::CpGrid.
template<typename Grid, typename GridView>
struct LookUpCellCentroid {
    /// \brief:     Constructor taking a GridView, CartesianMapper
    ///
    /// \param [in] GridView
    /// \param [in] CartesianIndexMapper
    /// \param [in] EclipseGrid
    explicit LookUpCellCentroid(const GridView& gridView,
                                const Dune::CartesianIndexMapper<Grid>& cartMapper,
                                const Opm::EclipseGrid* eclgrid) :
        gridView_(gridView),
        cartMapper_(&cartMapper),
        eclGrid_(eclgrid)
    {
    }

    const GridView& gridView_;
    const Dune::CartesianIndexMapper<Grid>* cartMapper_;
    const Opm::EclipseGrid* eclGrid_;

    /// \brief: Call operator
    ///
    ///         For grids different from Dune::CpGrid, it takes an element index, and
    ///         returns its cell centroid, from an EclipseGrid.
    ///
    /// \param [in] elemIdx      Element Index.
    /// \return     centroid     Centroid of the element, computed as in Eclipse.
    std::array<double,3> operator()(std::size_t elemIdx) const
    {
        if (std::is_same_v<Grid,Dune::CpGrid>)
        {
            OPM_THROW(std::logic_error, "Specialization for CpGrid must be used!");
        }
        return this -> eclGrid_ -> getCellCenter(this -> cartMapper_->cartesianIndex(elemIdx));
    }

}; // end struct LookUpCellCentroid

/// Spetialization for Dune::CpGrid
template<typename GridView>
struct LookUpCellCentroid<Dune::CpGrid,GridView>
{
    /// \brief:     Constructor taking a GridView, CartesianMapper.
    ///             Same argument as the general grid, to unify instantiations.
    ///
    /// \param [in] GridView
    /// \param [in] CartesianIndexMapper
    /// \param [in] EclipseGrid
    explicit LookUpCellCentroid(const GridView& gridView,
                                const Dune::CartesianIndexMapper<Dune::CpGrid>& cartMapper __attribute__((unused)),
                                const Opm::EclipseGrid* eclgrid __attribute__((unused))) :
        gridView_(gridView)
    {
    }

    const GridView& gridView_;
    static constexpr Dune::CartesianIndexMapper<Dune::CpGrid>* cartMapper_{nullptr};
    static constexpr Opm::EclipseGrid* eclGrid_{nullptr};

    /// \brief: Call operator
    ///
    ///         For Dune::CpGrid, it returns a function, taking an integer,
    ///         returning cell centroid, computed as in Eclipse.
    ///
    /// \param [in] elemIdx       Element Index.
    /// \return     centroid      Element centroid, computed as in Eclipse.
    std::array<double,3> operator()(std::size_t elemIdx) const
    {
        return this -> gridView_.grid().getEclCentroid(elemIdx);
    }
}; // end struct LookUpCellCentroid-Specialization for Dune::CpGrid.
}
// end namespace Opm
