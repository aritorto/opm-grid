//===========================================================================
//
// File: lookupdataCpGrid_test.cpp
//
// Created: Thurs 25.05.2023 16:05:00
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2023 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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
#include "config.h"

#define BOOST_TEST_MODULE LookUpDataCpGridTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/LookUpData.hh>

#include <dune/grid/common/mcmgmapper.hh>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }

    static int rank()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        return Dune::MPIHelper::instance(m_argc, m_argv).rank();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

void lookup_check(const Dune::CpGrid& grid)
{
    const auto& data = grid.data_;
    
    std::vector<int> fake_feature(data[0]->size(0), 0);
    std::iota(fake_feature.begin(), fake_feature.end(), 3);

    std::vector<double> fake_feature_double(data[0]->size(0), 0.);
    std::iota(fake_feature_double.begin(), fake_feature_double.end(), .5);

    // LookUpData
    const auto& leaf_view = grid.leafGridView();
    const Opm::LookUpData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>> lookUpData(leaf_view);
    // LookUpCartesianData
    const Dune::CartesianIndexMapper<Dune::CpGrid> cartMapper(grid);
    const Opm::LookUpCartesianData<Dune::CpGrid, Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>
        lookUpCartesianData(leaf_view, cartMapper);

    const auto& level0_view = grid.levelGridView(0);
    const Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView> leafMapper(leaf_view, Dune::mcmgElementLayout());
    const Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView> level0Mapper(level0_view, Dune::mcmgElementLayout());

    const auto& leaf_idSet = (*data.back()).local_id_set_;
    const auto& level0_idSet = (*data[0]).local_id_set_;

    for (const auto& elem : elements(leaf_view)) {
        // Search via Entity/Element
        const auto featureInElem = lookUpData(elem, fake_feature);
        const auto featureInElemDouble = lookUpData(elem, fake_feature_double);
        const auto featureInElemCartesian = lookUpCartesianData(elem, fake_feature);
        const auto featureInElemDoubleCartesian = lookUpCartesianData(elem,fake_feature_double);
        BOOST_CHECK(featureInElem == lookUpData.getOriginIndexFromEntity(elem) +3);
        BOOST_CHECK(featureInElemDouble == lookUpData.getOriginIndexFromEntity(elem) +.5);
        BOOST_CHECK(featureInElemCartesian == lookUpData.getOriginIndexFromEntity(elem) +3);
        BOOST_CHECK(featureInElemDoubleCartesian == lookUpData.getOriginIndexFromEntity(elem) +.5);
        BOOST_CHECK(elem.getOrigin().index() == lookUpData.getOriginIndexFromEntity(elem));
        // Search via INDEX
        const auto featureInElemIDX = lookUpData(elem.index(), fake_feature);
        const auto featureInElemDoubleIDX = lookUpData(elem.index(), fake_feature_double);
        const auto featureInElemCartesianIDX = lookUpCartesianData(elem.index(), fake_feature);
        const auto featureInElemDoubleCartesianIDX = lookUpCartesianData(elem.index(),fake_feature_double);
        BOOST_CHECK(featureInElemIDX == lookUpData.getOriginIndex(elem.index()) +3);
        BOOST_CHECK(featureInElemDoubleIDX == lookUpData.getOriginIndex(elem.index()) +.5);
        BOOST_CHECK(featureInElemCartesianIDX == lookUpData.getOriginIndex(elem.index()) +3);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == lookUpData.getOriginIndex(elem.index()) +.5);
        BOOST_CHECK(elem.getOrigin().index() == lookUpData.getOriginIndex(elem.index()));
        BOOST_CHECK(featureInElemIDX == featureInElem);
        BOOST_CHECK(featureInElemDoubleIDX == featureInElemDouble);
        BOOST_CHECK(featureInElemCartesianIDX == featureInElemCartesian);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == featureInElemDoubleCartesian);
        // Extra checks related to ElemMapper
        BOOST_CHECK(featureInElem == level0Mapper.index(elem.getOrigin()) +3);
        BOOST_CHECK(featureInElem == fake_feature[lookUpData.getOriginIndexFromEntity(elem)]);
        // Extra checks for  element index
        BOOST_CHECK(level0Mapper.index(elem.getOrigin()) == lookUpData.getOriginIndexFromEntity(elem));
        BOOST_CHECK(level0Mapper.index(elem.getOrigin()) == lookUpData.getOriginIndex(elem.index()));
        // Extra checks for cartesian element index
        BOOST_CHECK(cartMapper.cartesianIndex(elem.getOrigin().index()) == lookUpCartesianData.getCartesianOriginIndexFromEntity(elem));
        BOOST_CHECK(cartMapper.cartesianIndex(elem.getOrigin().index()) == lookUpCartesianData.getCartesianOriginIndex(elem.index()));
        if (elem.hasFather()) { // leaf_cell has a father!
            const auto& id = (*leaf_idSet).id(elem);
            const auto& parent_id = (*level0_idSet).id(elem.father());
            BOOST_CHECK(elem.index() == id);
            BOOST_CHECK(elem.index() == leafMapper.index(elem));
            BOOST_CHECK(elem.father().index() == featureInElem -3);
            BOOST_CHECK(elem.father().index() == parent_id);
            BOOST_CHECK(elem.father().index() == level0Mapper.index(elem.father()));
            BOOST_CHECK(elem.father().index() == lookUpData.getOriginIndexFromEntity(elem));
            // Extra checks for cartesian element index
            BOOST_CHECK(cartMapper.cartesianIndex(elem.father().index()) == lookUpCartesianData.getCartesianOriginIndexFromEntity(elem));
            BOOST_CHECK(cartMapper.cartesianIndex(elem.father().index()) == lookUpCartesianData.getCartesianOriginIndex(elem.index()));
        }
    }
}


BOOST_AUTO_TEST_CASE(one_lgr_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {3,2,3};  // patch_dim = {3-1, 2-0, 3-1} ={2,2,2}
    const std::string lgr_name = {"LGR1"};
    grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(single_cell_lgr_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::array<int, 3> cells_per_dim = {2,2,2};
    const std::array<int, 3> startIJK = {1,0,1};
    const std::array<int, 3> endIJK = {2,1,2};  // patch_dim = {2-1, 1-0, 2-1} ={1,1,1} -> Single Cell!
    const std::string lgr_name = {"LGR1"};
    grid.addLgrsUpdateLeafView({cells_per_dim}, {startIJK}, {endIJK}, {lgr_name});

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(lgrs_grid_A)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}, {4,4,4}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {0,0,2}, {3,2,2}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,1,1}, {1,1,3}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(lgrs_grid_B)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>> cells_per_dim_vec = {{2,2,2}, {3,3,3}};
    const std::vector<std::array<int,3>> startIJK_vec = {{0,0,0}, {3,2,0}};
    const std::vector<std::array<int,3>> endIJK_vec = {{2,2,1}, {4,3,3}};
    const std::vector<std::string> lgr_name_vec = {"LGR1", "LGR2"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(grid);
}


BOOST_AUTO_TEST_CASE(lgrs_grid_C)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {5,4,4};
    grid.createCartesian(grid_dim, cell_sizes);
    // Add LGRs and update LeafGridView
    const std::vector<std::array<int,3>>& cells_per_dim_vec = {{2,3,4}, {3,2,4}, {4,3,2}};
    const std::vector<std::array<int,3>>& startIJK_vec = {{0,0,0}, {4,0,0}, {4,3,3}};
    const std::vector<std::array<int,3>>& endIJK_vec = {{3,2,2}, {5,2,1}, {5,4,4}};
    const std::vector<std::string>& lgr_name_vec = {"LGR1", "LGR2", "LGR3"};
    grid.addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgr_name_vec);

    lookup_check(grid);
}

BOOST_AUTO_TEST_CASE(no_lgrs_grid)
{
    // Create a grid
    Dune::CpGrid grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    grid.createCartesian(grid_dim, cell_sizes);
    lookup_check(grid);
}
