//===========================================================================
//
// File: LookUpData.hpp
//
// Created: Tue May 23 14:44:00 2023
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
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp> // "CartesianIndexMapper not specialized for given grid"

namespace Dune
{
template <typename GridType>
class LookUpData
{
public:
    // Constructor taking a CpGrid object
    LookUpData(const GridType& grid) :
        gridView_(grid.leafGridView()),
        elemMapper_(gridView_, Dune::mcmgElementLayout()),
        cartMapper_(grid)
    {
    }

    template<typename feature_type>
    int operator()(const Dune::cpgrid::Entity<0>& elem, const std::vector<feature_type>& feature_vec)
    {
        // Assuming there is no LGR, so level 0 = leafview = "GLOBAL"
        return feature_vec[elemMapper_.index(elem)];
    }
protected:
    typename GridType::LeafGridView gridView_;
    Dune::MultipleCodimMultipleGeomTypeMapper<typename GridType::LeafGridView> elemMapper_;
    Dune::CartesianIndexMapper<GridType> cartMapper_; // currently not supported for general grid?


}; // end LookUpData class
}
// end namespace Dune
