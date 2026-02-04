/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE RefineMoreThanSixIntersectionsCellTest
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

struct Fixture
{
    Fixture()
    {
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        Dune::MPIHelper::instance(m_argc, m_argv);
        Opm::OpmLog::setupSimpleDefaultLogging();
    }
};

BOOST_GLOBAL_FIXTURE(Fixture);

template <class LeafView, class Entity>
int countIntersections(const LeafView& leafView, const Entity& element)
{
    int intersection_count = 0;
    for ([[maybe_unused]]const auto& intersection : Dune::intersections(leafView, element)){
        ++intersection_count;
    }
    return intersection_count;
}

void basicSizes(const Dune::CpGrid& grid)
{
    for (const auto& element : Dune::elements(grid.leafGridView())) {
       
        const auto& cell_to_face = grid.currentLeafData().cellToFace(element.index());
        BOOST_CHECK_EQUAL(countIntersections(grid.leafGridView(), element), cell_to_face.size());

        for (const auto& face : cell_to_face) {
            const auto& face_to_point = grid.currentLeafData().faceToPoint(face.index());
            BOOST_CHECK_EQUAL(face_to_point.size(), 4);
        }
    }
}

// Create a test grid from a simple deckstring.
// Refine it with one of:
// 0-> addLgrsUpdateLeafView
// 1-> globalRefine
// 2-> autoRefine
// 3-> adapt
void createAndCheckRefinedTestGrid(const std::string& deckString,
                                   int selectMethod)
{
    // Create the grid from string
    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseGrid eclGrid(deck);
    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    basicSizes(grid);


    for (const auto& element : Dune::elements(grid.leafGridView())) {
        // Check that each element has 7 intersections (see deck_string definition)
        BOOST_CHECK_EQUAL( countIntersections(grid.leafGridView(), element), 7 );

        // If selectMethod == 3-> adapt -> mark at least one element
        if ((selectMethod == 3) && (element.index()==0)){
            grid.mark(1,element, /* throwOnFailure = */ true);
        }
    }

    if (selectMethod == 0) {

        grid.addLgrsUpdateLeafView({{3,2,4}}, // cells_per_dim_vec
                                   {{0,0,0}}, // startIJK_vec
                                   {{1,1,1}}, // endIJK_vec
                                   {"LGR1"}); // lgr_name_vec

        Opm::checkGridWithLgrs(grid,
                               {{3,2,4}},  // cells_per_dim_vec
                               {"LGR1"},   // lgr_name_vec
                               false);    // gridHasBeenGlobalRefined

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            if (element.level() == 1) {
                BOOST_CHECK_EQUAL( countIntersections(grid.leafGridView(), element), 6 );
            }
            else {
                std::cout<< countIntersections(grid.leafGridView(), element) << " intersection for elemIdx " << element.index() <<std::endl;
                // prints 14!!! for element index 24 (equivalent leaf cell to level zero cell with index 1)
            }
        }
    }
    else if (selectMethod == 1){

        grid.globalRefine(1);

        Opm::checkGridWithLgrs(grid,
                               {{2,2,2}},  // cells_per_dim_vec
                               {"GR1"},   // lgr_name_vec
                               true);    // gridHasBeenGlobalRefined

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            BOOST_CHECK_EQUAL( countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
    }
    else if (selectMethod == 2) {

        grid.autoRefine(std::array{3,5,7}); // refinement factors in (x,y,z) directions

        // Extract the refined level grids name, excluding level zero grid name ("GLOBAL").
        // Note: in this case there is only one refined level grid, storing the global
        // refinement.
        std::vector<std::string> lgrNames(grid.maxLevel());
        for (const auto& [name, level] : grid.getLgrNameToLevel()) {
            if (level==0) { // skip level zero grid name for the checks
                continue;
            }
            lgrNames[level-1] = name; // Shift the index since level zero has been removed.
        }
        Opm::checkGridWithLgrs(grid,
                               /* cells_per_dim_vec = */ {{3,5,7}},
                               /* lgr_name_vec = */ lgrNames,
                               /* gridHasBeenGlobalRefined = */ true);

         for (const auto& element : Dune::elements(grid.leafGridView())) {
            BOOST_CHECK_EQUAL( countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
    }
    else if (selectMethod == 3) {

        grid.preAdapt();
        grid.adapt();
        grid.postAdapt();

        const auto& data = grid.currentData();
        BOOST_CHECK(static_cast<int>(data.size()) == grid.maxLevel() +2); // + level zero and leaf grids

        Opm::checkVertexAndFaceIndexAreNonNegative(grid); 
        Opm::checkGridBasicHiearchyInfo(grid, {{2,2,2}}, /*preAdaptMaxLevel*/ 0);
        Opm::checkGridLocalAndGlobalIdConsistency(grid, data);
        Opm::checkGlobalCellBounds(grid, data, /*lgrsHaveBlockShape*/ true, /*gridHasBeenGlobalRefined*/  false);

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            if (element.level() == 1) {
                BOOST_CHECK_EQUAL( countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
            }
            else {
                std::cout<< countIntersections(grid.leafGridView(), element) << " intersection for elemIdx " << element.index() <<std::endl;
                // 10 faces for element index 8 (equivalent leaf cell to level zero cell with index 1)
            }
        }
    }
    else if (selectMethod == 4) {

        // Cell 0 is the parent cell of LGR1 and Cell 1 is the parent cell of LGR2.
        // They "share I_FACE(s)" then the number of subdivisions in the y- and z- directions
        // should coincide. 
        grid.addLgrsUpdateLeafView({{3,2,4}, {2,2,4}}, // cells_per_dim_vec
                                   {{0,0,0}, {1,0,0}}, // startIJK_vec
                                   {{1,1,1}, {2,1,1}}, // endIJK_vec
                                   {"LGR1", "LGR2"}); // lgr_name_vec

        Opm::checkGridWithLgrs(grid,
                               {{3,2,4}, {2,2,4}},  // cells_per_dim_vec
                               {"LGR1", "LGR2"},   // lgr_name_vec
                               false);    // gridHasBeenGlobalRefined

         for (const auto& element : Dune::elements(grid.leafGridView())) {
            BOOST_CHECK_EQUAL( countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
    }
}


BOOST_AUTO_TEST_CASE(cellsWithMoreThanSixIntersections_fromDeck)
{
    // Level zero grid dims = 2x1x1
    //
    // cell 0
    // bottom face corners (0,0,0), (2,0,0), (0,2,0), (2,2,0)
    //    top face corners (0,0,2), (2,0,2), (0,2,2), (2,2,2)
    //
    // cell 1
    // bottom face corners (2,0,1), (4,0,1), (2,2,1), (4,2,1)
    //    top face corners (2,0,3), (4,0,3), (2,2,3), (4,2,3)

    const std::string deckString =
        R"(RUNSPEC
DIMENS
 2 1 1 /

GRID

COORD
0 0 0   0 0 3
2 0 0   2 0 3
4 0 0   4 0 3

0 2 0   0 2 3
2 2 0   2 2 3
4 2 0   4 2 3
/

ZCORN
0 0 1 1  0 0 1 1
2 2 3 3  2 2 3 3
/

ACTNUM
2*1
/
)";

    createAndCheckRefinedTestGrid(deckString, 0); // 0-> refinement via addLgrsUpdateLeafView (1 LGR)
    createAndCheckRefinedTestGrid(deckString, 1); // 1-> refinement via globalRefine
    createAndCheckRefinedTestGrid(deckString, 2); // 2-> refinement via autoRefine
    createAndCheckRefinedTestGrid(deckString, 3); // 3-> refinement via adpat
    createAndCheckRefinedTestGrid(deckString, 4); // 4-> refinement via adddLgrsUpdateLeafView (2 LGRs)
}


BOOST_AUTO_TEST_CASE(cellsWithMoreThanSixIntersections, *boost::unit_test::disabled())
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {2,1,1}, /* cell_sizes = */ {2.,2.,2.});
    //                          LGR1 parent cells
    // --------------------|    ------------------
    // k = 0       | 0  1  |    | - 1|
    //---------------------|
    // Refine cell 1 into 2x2x2 children (2 subdivisions per x-,y-,z- direction).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{1,1,2}},
                               /* startIJK_vec = */ {{1,0,0}},
                               /* endIJK_vec = */ {{2,1,1}},
                               /* lgr_name_vec = */ {"LGR1"});
    // Leaf grid view
    //  ----------
    //  |    | 2 |
    //  | 0  |---|
    //  |    | 1 |
    //  ----------

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        if (element.index()!=0)
            continue;
        grid.mark(1, element, true);
    }
    const auto& cell_to_face = grid.currentLeafData().cellToFace(0);
    std::cout<< "count faces: " << cell_to_face.size() << std::endl;

    const auto& cell_to_point = grid.currentLeafData().cellToPoint(0);
    std::cout<< "count point: " << cell_to_point.size() << std::endl;

    grid.preAdapt();
    grid.adapt();
    grid.postAdapt();

    const auto& data = grid.currentData();
    BOOST_CHECK(static_cast<int>(data.size()) == grid.maxLevel() +2);

    Opm::checkVertexAndFaceIndexAreNonNegative(grid); // FAILS - for now
    Opm::checkGridBasicHiearchyInfo(grid, {{2,2,2}}, /*preAdaptMaxLevel*/ 1);
    Opm::checkGridLocalAndGlobalIdConsistency(grid, data);
    Opm::checkGlobalCellBounds(grid, data, /*lgrsHaveBlockShape*/ true, /*gridHasBeenGlobalRefined*/  false);
}

