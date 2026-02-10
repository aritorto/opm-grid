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
#include <opm/grid/cpgrid/CpGridData.hpp>
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





/*template<typename FaceToPoint>
  bool coarserFaceContainsFinerFace(const Dune::cpgrid::CpGridData& coarserCpData,
                                  const FaceToPoint& coarserFace_to_point,
                                  const FaceToPoint& finerace_to_point,
                                  const Dune::cpgrid::CpGridData& finerCpData)
{

    const auto c0 =  coarserCpData.vertexPosition(coarserCpData.faceToPoint(face.index())[0]);
                 const auto c1 =  coarserCpData.vertexPosition(.faceToPoint(face.index())[1]);
                  const auto c2 =  grid.vertexPosition(grid.currentLeafData().faceToPoint(face.index())[2]);
                   const auto c3 =  grid.vertexPosition(grid.currentLeafData().faceToPoint(face.index())[3]);
                std::cout<< c0[0] << " " << c0[1] << " " << c0[2] << std::endl;
                 std::cout<< c1[0] << " " << c1[1] << " " << c1[2] << std::endl;
                  std::cout<< c2[0] << " " << c2[1] << " " << c2[2] << std::endl;
                  std::cout<< c3[0] << " " << c3[1] << " " << c3[2] << std::endl;
    
    return true;
    }*/

/*std::array<std::vector<int>, 6> classifyAndCollectFaceIndices(const Dune::CpGrid& grid,
                                                              const Dune::cpgrid::Entity<0>& element)
{
    std::array<std::vector<int>, 6> classified_face_idxs{};
    // clasified_face_idxs[0] stores I false face indices
    // clasified_face_idxs[1] stores I true  face indices
    // clasified_face_idxs[2] stores J false face indices
    // clasified_face_idxs[3] stores J true  face indices
    // clasified_face_idxs[4] stores K false face indices
    // clasified_face_idxs[5] stores K true  face indices

    for (int i = 0; i < 6; ++i){
        classified_face_idxs[i].reserve(Opm::Lgr::countIntersections(grid.leafGridView(), element)); // more than needed
    }
    
    
    const auto& cell_to_face = grid.currentLeafData().cellToFace(element.index());
    for (const auto& face : cell_to_face) {
        const auto tag = grid.currentLeafData().faceTag(face.index());
        const bool orientation = face.orientation();
        if ((tag == I_FACE) && !orientation) { 
            classified_face_idxs[0].push_back(face.index());
            if (element.level() == 0) 
            {
                std::cout<< face.index() << " faceIdx, " << tag << " faceTag, " << orientation << ", elemIdx: " << element.index() <<std::endl;
                std::cout<< grid.faceArea(face.index()) << " faceArea " << std::endl;
                const auto c0 =  grid.vertexPosition(grid.currentLeafData().faceToPoint(face.index())[0]);
                 const auto c1 =  grid.vertexPosition(grid.currentLeafData().faceToPoint(face.index())[1]);
                  const auto c2 =  grid.vertexPosition(grid.currentLeafData().faceToPoint(face.index())[2]);
                   const auto c3 =  grid.vertexPosition(grid.currentLeafData().faceToPoint(face.index())[3]);
                std::cout<< c0[0] << " " << c0[1] << " " << c0[2] << std::endl;
                 std::cout<< c1[0] << " " << c1[1] << " " << c1[2] << std::endl;
                  std::cout<< c2[0] << " " << c2[1] << " " << c2[2] << std::endl;
                   std::cout<< c3[0] << " " << c3[1] << " " << c3[2] << std::endl;
            }
        }
        else if ((tag == I_FACE) && orientation) {
            classified_face_idxs[1].push_back(face.index());
        }
        else if ((tag == J_FACE) && !orientation) {
            classified_face_idxs[2].push_back(face.index());
        }
        else if ((tag == J_FACE) && orientation) {
            classified_face_idxs[3].push_back(face.index());
        }
        else if ((tag == K_FACE) && !orientation) {
            classified_face_idxs[4].push_back(face.index());
        }
        else if ((tag == K_FACE) && orientation) {
            classified_face_idxs[5].push_back(face.index());
        }
        else {
            std::cout<< face.index() << " faceIdx, " << tag << " faceTag, " << orientation << std::endl;
            std::cout<< "why are we here?" << std::endl;
        }
    }
    return classified_face_idxs;
}
*/
bool facePointsBelongToCellPoints(const Dune::CpGrid& grid, const Dune::cpgrid::Entity<0>& element)
{
    const auto& cellFaces = grid.currentLeafData().cellToFace(element.index());
    const auto& cellPoints = grid.currentLeafData().cellToPoint(element.index());
     
    for (const auto& face : cellFaces) {
        int facePointsCount = grid.numFaceVertices(face.index());
        for (int i = 0; i < facePointsCount; ++i) {
            int pointIdx  = grid.faceVertex(face.index(), i);
            //  std::cout<<   (std::find(cellPoints.begin(), cellPoints.end(),
            //                 pointIdx) != cellPoints.end() ) << " pointIdx " << pointIdx << ", faceIdx: "<< face.index() <<std::endl;
        }
    }
    return true;
}

void basicSizes(const Dune::CpGrid& grid)
{
    for (const auto& element : Dune::elements(grid.leafGridView())) {
       
        const auto& cell_to_face = grid.currentLeafData().cellToFace(element.index());
        BOOST_CHECK_EQUAL(Opm::Lgr::countIntersections(grid.leafGridView(), element), cell_to_face.size());

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
        BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 7 );

        const auto classified_faces = Opm::Lgr::classifyAndCollectFaceIndices(grid.currentLeafData(), element);
        // Check that element with index 0 has 2 {I_FACE, true} faces
        if (element.index() == 0) {
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE false */ 0].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE true  */ 1].size(), 2);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE false */ 2].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE true  */ 3].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE false */ 4].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE true  */ 5].size(), 1);

            facePointsBelongToCellPoints(grid, element);
        }
        else if (element.index() == 1) { // Check that element with index 1 has 2 {I_FACE, false} faces
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE false */ 0].size(), 2);
            BOOST_CHECK_EQUAL( classified_faces[/* I_FACE true  */ 1].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE false */ 2].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* J_FACE true  */ 3].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE false */ 4].size(), 1);
            BOOST_CHECK_EQUAL( classified_faces[/* K_FACE true  */ 5].size(), 1);

            facePointsBelongToCellPoints(grid, element);
        }
      
      

        // If selectMethod == 3-> adapt -> mark at least one element
        if ((selectMethod == 3) && (element.index()==0)){
            grid.mark(1,element, /* throwOnFailure = */ true);
        }
    }

    if (selectMethod == 0) {

        grid.addLgrsUpdateLeafView({{2,2,3}}, // cells_per_dim_vec
                                   {{0,0,0}}, // startIJK_vec
                                   {{1,1,1}}, // endIJK_vec
                                   {"LGR1"}); // lgr_name_vec

        Opm::checkGridWithLgrs(grid,
                               {{2,2,3}},  // cells_per_dim_vec
                               {"LGR1"},   // lgr_name_vec
                               false);    // gridHasBeenGlobalRefined

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            if (element.level() == 1) {
                BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 6 );
            }
            else {
                std::cout<< Opm::Lgr::countIntersections(grid.leafGridView(), element) << " intersection for elemIdx " << element.index() <<std::endl;
                // prints 14!!! for element index 24 (equivalent leaf cell to level zero cell with index 1)

                const auto classified_faces = Opm::Lgr::classifyAndCollectFaceIndices(grid.currentLeafData(), element);
      
                std::cout<< classified_faces[/* I_FACE false */ 0].size() << " I false"<<std::endl; // 7
                std::cout<< classified_faces[/* I_FACE true  */ 1].size() << " I true" <<std::endl; // 1
                std::cout<< classified_faces[/* J_FACE false */ 2].size() << " J false"<<std::endl; // 1
                std::cout<< classified_faces[/* J_FACE true  */ 3].size() << " J true"<<std::endl;  // 1
                std::cout<< classified_faces[/* K_FACE false */ 4].size() << " K false"<<std::endl; // 1
                std::cout<< classified_faces[/* K_FACE true  */ 5].size() << " K true"<<std::endl;  // 1

                facePointsBelongToCellPoints(grid, element);
            }
        }
        basicSizes(grid);
    }
    else if (selectMethod == 1){

        grid.globalRefine(1);

        Opm::checkGridWithLgrs(grid,
                               {{2,2,2}},  // cells_per_dim_vec
                               {"GR1"},   // lgr_name_vec
                               true);    // gridHasBeenGlobalRefined

        for (const auto& element : Dune::elements(grid.leafGridView())) {
            BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
        basicSizes(grid);
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
            BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
        basicSizes(grid);
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
                BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
            }
            else {
                std::cout<< Opm::Lgr::countIntersections(grid.leafGridView(), element) << " intersection for elemIdx " << element.index() <<std::endl;
                // 10 faces for element index 8 (equivalent leaf cell to level zero cell with index 1)
            }
        }
        basicSizes(grid);
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
            BOOST_CHECK_EQUAL( Opm::Lgr::countIntersections(grid.leafGridView(), element), 6 ); // suspicious?
        }
        basicSizes(grid);
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
    //  createAndCheckRefinedTestGrid(deckString, 1); // 1-> refinement via globalRefine
    //  createAndCheckRefinedTestGrid(deckString, 2); // 2-> refinement via autoRefine
    // createAndCheckRefinedTestGrid(deckString, 3); // 3-> refinement via adpat
    // createAndCheckRefinedTestGrid(deckString, 4); // 4-> refinement via adddLgrsUpdateLeafView (2 LGRs)
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

