//===========================================================================
//
// File: globalRefine_test.cpp
//
// Created: Fri 2 Feb 2023
//
// Author(s): Antonella Ritorto   <antonella.ritorto@opm-op.com>
//
// $Date$
//
// $Revision$
//
//===========================================================================
/*
  Copyright 2024 Equinor ASA.

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

#define BOOST_TEST_MODULE GlobalRefineTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/DefaultGeometryPolicy.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/LookUpData.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <sstream>
#include <iostream>

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

#define CHECK_COORDINATES(c1, c2)                                       \
    for (int c = 0; c < 3; c++) {                                       \
        BOOST_TEST(c1[c] == c2[c], boost::test_tools::tolerance(1e-12)); \
    }




void check_global_refine(const Dune::CpGrid& refined_grid, const Dune::CpGrid& equiv_fine_grid)
{

    const auto& refined_data = refined_grid.data_;
    const auto& equiv_data = equiv_fine_grid.data_;

    const auto& refined_leaf = *refined_data.back();
    const auto& equiv_leaf = *equiv_data[0];

    // Check the container sizes
    BOOST_CHECK_EQUAL(refined_leaf.face_to_cell_.size(), equiv_leaf.face_to_cell_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_to_point_.size(), equiv_leaf.face_to_point_.size());
    BOOST_CHECK_EQUAL(refined_leaf.cell_to_point_.size(), equiv_leaf.cell_to_point_.size());
    BOOST_CHECK_EQUAL(refined_leaf.face_normals_.size(), equiv_leaf.face_normals_.size());

    // Check that the points (ordering/coordinates) matches
    auto equiv_point_iter = equiv_leaf.geomVector<3>().begin();
    for(const auto& point: refined_leaf.geomVector<3>())
    {
        CHECK_COORDINATES(point.center(), equiv_point_iter->center());
        for(const auto& coord: point.center())
            BOOST_TEST(std::isfinite(coord));
        ++equiv_point_iter;
    }
    auto equiv_cell_iter = equiv_leaf.geomVector<3>().begin();
    for(const auto& cell: refined_leaf.geomVector<3>())
    {
        CHECK_COORDINATES(cell.center(), equiv_cell_iter->center());
        for(const auto& coord: cell.center())
            BOOST_TEST(std::isfinite(coord));
        BOOST_CHECK_CLOSE(cell.volume(), equiv_cell_iter->volume(), 1e-6);
        ++equiv_cell_iter;
    }
    
    /////
    const auto& grid_view = refined_grid.leafGridView();
    const auto& equiv_grid_view = equiv_fine_grid.leafGridView();

    auto equiv_element_iter = equiv_grid_view.begin<0>();
    for(const auto& element: elements(grid_view))
    {
        BOOST_CHECK( element.getOrigin().level() == 0);
        //BOOST_CHECK( element.getOrigin().index() == element.index());
        for(const auto& intersection: intersections(grid_view, element))
        {
            // find matching intersection (needed as ordering is allowed to be different
            bool matching_intersection_found = false;
            for(auto& intersection_match: intersections(equiv_grid_view, *equiv_element_iter))
            {
                if(intersection_match.indexInInside() == intersection.indexInInside())
                {
                    BOOST_CHECK(intersection_match.neighbor() == intersection.neighbor());

                    if(intersection.neighbor())
                    {
                        BOOST_CHECK(intersection_match.indexInOutside() == intersection.indexInOutside());
                    }

                    CHECK_COORDINATES(intersection_match.centerUnitOuterNormal(), intersection.centerUnitOuterNormal());
                    const auto& geom_match = intersection_match.geometry();
                    BOOST_TEST(0.0 == 1e-11, boost::test_tools::tolerance(1e-8));
                    const auto& geom =  intersection.geometry();
                    BOOST_CHECK_CLOSE(geom_match.volume(), geom.volume(), 1e-6);
                    CHECK_COORDINATES(geom_match.center(), geom.center());
                    BOOST_CHECK(geom_match.corners() == geom.corners());

                    decltype(geom.corner(0)) sum_match{}, sum{};

                    for(int cor = 0; cor < geom.corners(); ++cor)
                    {
                        sum += geom.corner(cor);
                        sum_match += geom_match.corner(1);
                    }
                    CHECK_COORDINATES(sum, sum_match);
                    matching_intersection_found = true;
                    break;
                }
            }
            BOOST_CHECK(matching_intersection_found);
        }
        ++equiv_element_iter;
        }
    /////
}

BOOST_AUTO_TEST_CASE(globalRefineOneLgr)
{
    // Create a 4x4x2 grid with length 4x4x2
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,4,2};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    coarse_grid.globalRefine(1);
    // coarse_grid.addLgrUpdateLeafView({4,4,2}, {0,0,0}, {4,4,2},"LGR");
    
    // Create a 16x16x4 grid with length 4x4x2
    Dune::CpGrid fine_grid;
    const std::array<double, 3> fine_cell_sizes = {0.25, 0.25, 0.5};
    const std::array<int, 3> fine_grid_dim = {16,16,4};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}

BOOST_AUTO_TEST_CASE(globalRefineNo)
{
    // Create a 4x3x3 grid with length 4x3x3
    Dune::CpGrid coarse_grid;
    const std::array<double, 3> cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> grid_dim = {4,3,3};
    coarse_grid.createCartesian(grid_dim, cell_sizes);
    coarse_grid.globalRefine(0);

    // Create a 4x3x3 grid with length 4x3x3
    Dune::CpGrid fine_grid;
    const std::array<double, 3> fine_cell_sizes = {1.0, 1.0, 1.0};
    const std::array<int, 3> fine_grid_dim = {4,3,3};
    fine_grid.createCartesian(fine_grid_dim, fine_cell_sizes);

    check_global_refine(coarse_grid, fine_grid);
}

