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

#define BOOST_TEST_MODULE FaultFullyInLgrTest
#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/CpGridUtilities.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <array>
#include <algorithm>

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

/*bool facePointsBelongToCellPoints(const Dune::CpGrid& grid, const Dune::cpgrid::Entity<0>& element)
{
    const auto& cellFaces = grid.currentLeafData().cellToFace(element.index());
    const auto& cellPoints = grid.currentLeafData().cellToPoint(element.index());

    bool allBelong = false;
    for (const auto& face : cellFaces) {
        int facePointsCount = grid.numFaceVertices(face.index());
        for (int i = 0; i < facePointsCount; ++i) {
            int pointIdx  = grid.faceVertex(face.index(), i);
            allBelong = std::find(cellPoints.begin(), cellPoints.end(),pointIdx) != cellPoints.end();
            if (!allBelong)
                return false;
        }
    }
    return true;
    }*/

void basicSizes(const Dune::CpGrid& grid)
{
    for (const auto& element : Dune::elements(grid.leafGridView())) {
       
        const auto& cell_to_face = grid.currentLeafData().cellToFace(element.index());
        BOOST_CHECK_EQUAL(Opm::Lgr::countIntersections(grid.leafGridView(), element), cell_to_face.size());

        for (const auto& face : cell_to_face) {
            const auto& face_to_point = grid.currentLeafData().faceToPoint(face.index());
            BOOST_CHECK_EQUAL(face_to_point.size(), 4); // not necessarily 
        }
    }
}

void checkParentCellMoreThanSixIntersections(const auto& classified_faces, int elemIdx)
{
    
        // Check that element with index 0 has 2 {I_FACE, true} faces
        if (elemIdx == 0) {
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE false */ 0].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE true  */ 1].size(), 2);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE false */ 2].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE true  */ 3].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE false */ 4].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE true  */ 5].size(), 1);

            for (const auto& itrue : classified_faces[1])
            {
                std::cout<< "soy una iTrue in level 0, elem 0: " << itrue << std::endl;
            }
            
        }
        else if (elemIdx == 1) { // Check that element with index 1 has 2 {I_FACE, false} faces
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE false */ 0].size(), 2);
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE true  */ 1].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE false */ 2].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE true  */ 3].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE false */ 4].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE true  */ 5].size(), 1);

             for (const auto& ifalse : classified_faces[0])
            {
                std::cout<< "soy una ifalse in level 0, elem 1: " << ifalse << std::endl;
            }
        }
}

void createAndCheckRefinedTestGrid(const std::string& deckString)
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
        // (Grid dim 2x1x1 with cell0 2 iTrue faces and cell1 with 2 iFalse faces)
        BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 7 );

        const auto classified_faces = Opm::Lgr::classifyAndCollectFaceIndices(grid.currentLeafData(), element);
        checkParentCellMoreThanSixIntersections(classified_faces, element.index());
      
    }
        
   
    
        

        // Cell 0 is the parent cell of LGR1 and Cell 1 is the parent cell of LGR2.
        // They "share I_FACE(s)" then the number of subdivisions in the y- and z- directions
        // should coincide. 
        grid.addLgrsUpdateLeafView({{3,2,2}, {2,2,2}}, // cells_per_dim_vec
                                   {{0,0,0}, {1,0,0}}, // startIJK_vec
                                   {{1,1,1}, {2,1,1}}, // endIJK_vec
                                   {"LGR1", "LGR2"}); // lgr_name_vec

        const auto& [nx1, ny1, nz1] = grid.currentData()[1]->logicalCartesianSize();
         const auto& [nx2, ny2, nz2] = grid.currentData()[2]->logicalCartesianSize();
         const auto& [nx, ny, nz] = grid.logicalCartesianSize();
         std::cout<< nx << " " << ny << " " << nz << " grid cart size " <<std::endl;
         std::cout<< nx1 << " " << ny1 << " " << nz1 << " LGR1 cart size " <<std::endl;
         std::cout<< nx2 << " " << ny2 << " " << nz2 << " LGR2 cart size " <<std::endl;
             
        
        Opm::checkGridWithLgrs(grid,
                               {{3,2,2}, {2,2,2}},  // cells_per_dim_vec
                               {"LGR1", "LGR2"},   // lgr_name_vec
                               false);    // gridHasBeenGlobalRefined

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
        basicSizes(grid);

        for (const auto& element : Dune::elements(grid.leafGridView())) {
          
                std::cout<< Opm::Lgr::countIntersections(grid.leafGridView(), element) << " intersection for elemIdx " << element.index() <<std::endl;
                // prints 14!!! for element index 24 (equivalent leaf cell to level zero cell with index 1)

                const auto classified_faces = Opm::Lgr::classifyAndCollectFaceIndices(grid.currentLeafData(), element);
      
                std::cout<< classified_faces[/* I_FACE false */ 0].size() << " I false"<<std::endl; // 7
                std::cout<< classified_faces[/* I_FACE true  */ 1].size() << " I true" <<std::endl; // 1
                std::cout<< classified_faces[/* J_FACE false */ 2].size() << " J false"<<std::endl; // 1
                std::cout<< classified_faces[/* J_FACE true  */ 3].size() << " J true"<<std::endl;  // 1
                std::cout<< classified_faces[/* K_FACE false */ 4].size() << " K false"<<std::endl; // 1
                std::cout<< classified_faces[/* K_FACE true  */ 5].size() << " K true"<<std::endl;  // 1
            }
        basicSizes(grid);
}


BOOST_AUTO_TEST_CASE(cellsWithMoreThanSixIntersections_fromDeck)
{
    // Level zero grid dims = 2x1x1
    //
    // cell 0
    // bottom face corners (0,0,0), (6,0,0), (0,6,0), (6,6,0)
    //    top face corners (0,0,8), (6,0,8), (0,6,8), (6,6,8)
    //
    // cell 1
    // bottom face corners (6,0,1), (12,0,1), (6,6,1), (12,6,1)
    //    top face corners (6,0,9), (12,0,9), (12,6,9), (12,6,9)
    

    const std::string deckString =
        R"(RUNSPEC
DIMENS
 2 1 1 /

GRID

COORD
 0 0 0     0 0 9
 6 0 0     6 0 9
12 0 0    12 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9
12 6 0   12 6 9
/

ZCORN
0 0 1 1  0 0 1 1
8 8 9 9  8 8 9 9
/

ACTNUM
2*1
/
)";
    createAndCheckRefinedTestGrid(deckString); // refinement via adddLgrsUpdateLeafView (2 LGRs)
}
