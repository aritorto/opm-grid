//===========================================================================
//
// File: test_lookupdata_polyhedral.cpp
//
// Created: Monday 24.07.2023 14:30:00
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

#define BOOST_TEST_MODULE LookUpDataPolyhedralGridTest
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/LookUpData.hh>

#include <dune/common/unused.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/grid/cpgrid/dgfparser.hh>
#include <opm/grid/polyhedralgrid/dgfparser.hh>

#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

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

void lookup_check(const Dune::PolyhedralGrid<3,3>& grid)
{
    std::vector<int> fake_feature(grid.size(0), 0);
    std::iota(fake_feature.begin(), fake_feature.end(), 3);

    std::vector<double> fake_feature_double(grid.size(0), 0.);
    std::iota(fake_feature_double.begin(), fake_feature_double.end(), .5);

    const auto& leaf_view = grid.leafGridView();
    using GridView = std::remove_cv_t< typename std::remove_reference<decltype(grid.leafGridView())>::type>;
    // LookUpData
    const Opm::LookUpData<Dune::PolyhedralGrid<3,3>, GridView> lookUpData(leaf_view);
    // LookUpCartesianData
    const Dune::CartesianIndexMapper<Dune::PolyhedralGrid<3,3>> cartMapper(grid);
    const Opm::LookUpCartesianData<Dune::PolyhedralGrid<3,3>, GridView> lookUpCartesianData(leaf_view, cartMapper);
    // Mapper
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(grid.leafGridView(), Dune::mcmgElementLayout());

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
        // Search via INDEX
        const auto idx = mapper.index(elem);
        const auto featureInElemIDX = lookUpData(idx, fake_feature);
        const auto featureInElemDoubleIDX = lookUpData(idx, fake_feature_double);
        const auto featureInElemCartesianIDX = lookUpCartesianData(idx, fake_feature);
        const auto featureInElemDoubleCartesianIDX = lookUpCartesianData(idx, fake_feature_double);
        BOOST_CHECK(featureInElemIDX == (lookUpData.getOriginIndex<Dune::PolyhedralGrid<3,3>>(idx))+3);
        BOOST_CHECK(featureInElemDoubleIDX == (lookUpData.getOriginIndex<Dune::PolyhedralGrid<3,3>>(idx)) +.5);
        BOOST_CHECK(featureInElemCartesianIDX == (lookUpData.getOriginIndex<Dune::PolyhedralGrid<3,3>>(idx)) +3);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == (lookUpData.getOriginIndex<Dune::PolyhedralGrid<3,3>>(idx)) +.5);
        BOOST_CHECK(idx == (lookUpData.getOriginIndex<Dune::PolyhedralGrid<3,3>>(idx)));
        BOOST_CHECK(featureInElemIDX == featureInElem);
        BOOST_CHECK(featureInElemDoubleIDX == featureInElemDouble);
        BOOST_CHECK(featureInElemCartesianIDX == featureInElemCartesian);
        BOOST_CHECK(featureInElemDoubleCartesianIDX == featureInElemDoubleCartesian);
        // Extra check for  element index
        BOOST_CHECK(idx == lookUpData.getOriginIndexFromEntity(elem));
        // Extra checks for cartesian element index
        BOOST_CHECK(cartMapper.cartesianIndex(idx) == lookUpCartesianData.getCartesianOriginIndexFromEntity(elem));
        BOOST_CHECK(cartMapper.cartesianIndex(idx) == lookUpCartesianData.getCartesianOriginIndex(idx));
    }

}

BOOST_AUTO_TEST_CASE(PolyGridFromEcl)
{
#if HAVE_ECL_INPUT
    const char *deckString =
        "RUNSPEC\n"
        "METRIC\n"
        "DIMENS\n"
        "4 4 4 /\n"
        "GRID\n"
        "DXV\n"
        "4*1 /\n"
        "DYV\n"
        "4*1 /\n"
        "DZ\n"
        "16*1 /\n"
        "TOPS\n"
        "16*100.0 /\n";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseGrid eclgrid( deck);
    std::vector<double> porv;

    Dune::PolyhedralGrid<3,3> grid(eclgrid, porv);
    lookup_check(grid);
#endif
}

