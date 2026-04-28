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

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/CpGridUtilities.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
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

template<typename Coordinate> // Local/Global
Coordinate crossProduct(const Coordinate& v, const Coordinate& w) {
    return { v[1]*w[2] - v[2] * w[1],
             v[2] * w[0] - v[0] * w[2],
             v[0] * w[1] - v[1] * w[0]};
}

template<typename Coordinate>
double dotProduct(const Coordinate& v, const Coordinate& w) {
    return v[0]*w[0]+ v[1]*w[1] + v[2]*w[2];
}

template<typename Coordinate>
std::pair<bool,Coordinate> computeSegmentIntersection(const Coordinate& startSegmentA, const Coordinate& endSegmentA,
                                                      const Coordinate& startSegmentB, const Coordinate& endSegmentB,
                                                      bool& isInteriorInSegmentA,
                                                      bool& isInteriorInSegmentB)
{
    const auto directionA = endSegmentA - startSegmentA;
    const auto directionB = endSegmentB - startSegmentB;

    const auto c = crossProduct(directionA, directionB); // dotProduct against directionA/B is zero

    // Segment A: startSegmentA + t*directionA, t in [0,1]
    // Segment B: startSegmentB + s*directionB, s in [0,1]

    const auto coplanar = dotProduct(startSegmentB - startSegmentA, c);
    if (coplanar > 1e-8){ // segments are skew (no intersection)
        return std::make_pair<bool,Coordinate>(false, Coordinate{});
    }
    else { // If segments are parallel, there is no intersection
        if (c.two_norm() < 1e-8) {
            return std::make_pair<bool, Coordinate>(false, Coordinate{});
        }
        else {
            double t = dotProduct(crossProduct(startSegmentB - startSegmentA, directionB), c) / dotProduct(c,c);
            double s = dotProduct(crossProduct(startSegmentB - startSegmentA, directionA), c) / dotProduct(c,c);

            if ((t >= 0) && (t <= 1) && (s >= 0) && (s <= 1)) { // segments intersect
                isInteriorInSegmentA = (t > 0) && (t < 1);
                isInteriorInSegmentB = (s > 0) && (s < 1);
                return std::make_pair<bool,Coordinate>(true, startSegmentA + (t*directionA));
            } else { // lines intersect, but not the segments
                return std::make_pair<bool,Coordinate>(false, Coordinate{});
            }
        }
    }
}

template<typename Coordinate>
std::vector<std::array<Coordinate,2>> createEdges(const Coordinate& f0, const Coordinate& f1,
                                                  const Coordinate& f2, const Coordinate& f3)
{
    return std::vector<std::array<Coordinate,2>>{ {f0,f1}, {f1,f2}, {f2,f3}, {f3,f0} };
}

void printIfNewVertex(const Dune::cpgrid::CpGridData& singleCellRefinementData, // (refinement of the parent cell)
                      const Dune::cpgrid::Entity<0>& refinedElem,
                      const Dune::cpgrid::CpGridData& parentGridData,
                      const Dune::cpgrid::Entity<0>& parentElem)
{
    const auto& refinedCellToFace = singleCellRefinementData.cellToFace(refinedElem.index());
    const auto& parentCellToFace = parentGridData.cellToFace(parentElem.index());

    for (const auto& face : parentCellToFace) {

        const auto& faceToPoint = parentGridData.faceToPoint(face.index());

        const auto& f0 =  Dune::cpgrid::Entity<3>(parentGridData, faceToPoint[0], true).geometry().center();
        const auto& f1 =  Dune::cpgrid::Entity<3>(parentGridData, faceToPoint[1], true).geometry().center();
        const auto& f2 =  Dune::cpgrid::Entity<3>(parentGridData, faceToPoint[2], true).geometry().center();
        const auto& f3 =  Dune::cpgrid::Entity<3>(parentGridData, faceToPoint[3], true).geometry().center();

        const auto edges = createEdges(f0,f1,f2,f3);

        const auto& faceTag = parentGridData.faceTag(face.index());
        const auto& faceOrientation = face.orientation();

        for (const auto& refinedFace : refinedCellToFace) {
            // Skip face if it is not on the boundary of the single cell refinement grid
            if (singleCellRefinementData.faceToCellSize(refinedFace.index()) != 1)
                continue;


            const auto& refinedFaceToPoint = singleCellRefinementData.faceToPoint(refinedFace.index());

            const auto& rf0 =  Dune::cpgrid::Entity<3>(singleCellRefinementData, refinedFaceToPoint[0], true).geometry().center();
            const auto& rf1 =  Dune::cpgrid::Entity<3>(singleCellRefinementData, refinedFaceToPoint[1], true).geometry().center();
            const auto& rf2 =  Dune::cpgrid::Entity<3>(singleCellRefinementData, refinedFaceToPoint[2], true).geometry().center();
            const auto& rf3 =  Dune::cpgrid::Entity<3>(singleCellRefinementData, refinedFaceToPoint[3], true).geometry().center();

            const auto edges_rf = createEdges(rf0,rf1,rf2,rf3);

            const auto& refinedFaceTag = singleCellRefinementData.faceTag(refinedFace.index());
            const auto& refinedFaceOrientation = refinedFace.orientation();

            if ((refinedFaceTag == faceTag) && (refinedFaceOrientation == faceOrientation)) {



                for (const auto& edge : edges) {
                    for (const auto& edge_rf : edges_rf) {
                        bool isInteriorInEdge{};
                        bool isInteriorInEdgeRf{};

                        const auto [found, segmentInter] = computeSegmentIntersection(edge[0], edge[1],
                                                                                      edge_rf[0], edge_rf[1],
                                                                                      isInteriorInEdge,
                                                                                      isInteriorInEdgeRf);

                        if (found && isInteriorInEdge && isInteriorInEdgeRf){
                            std::cout<< segmentInter[0] << " " << segmentInter[1] << " " << segmentInter[2] << "  segment intersection found " <<std::endl;
                            std::cout<< face.index() << " coarse face index, refined face index: " << refinedFace.index() << std::endl;
                        }

                    }
                }
            }
        }
    }
}


BOOST_AUTO_TEST_CASE(cellsWithMoreThanSixIntersections_I_FACE)
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

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{2,3,2}},
                              /* startIJK_vec */ {{0,0,0}},
                              /* endIJK_vec */ {{1,1,1}},
                              /* lgr_name_vec */ {"LGR1"});

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];

    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {
        printIfNewVertex(refinedGridData, refinedElem, parentGridData, parentElem);
    }
}


BOOST_AUTO_TEST_CASE(cellsWithMoreThanSixIntersections_J_FACE)
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
 1 2 1 /

GRID

COORD
 0 0 0    0 0 9
 6 0 0    6 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9

 0 12 0   0 12 9
 6 12 0   6 12 9
/

ZCORN
0 0 1 1  0 0 1 1
8 8 9 9  8 8 9 9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridAndAddLgrs(grid,
                              deckString,
                              /* cells_per_dim_vec */ {{3,2,2}},
                              /* startIJK_vec */ {{0,0,0}},
                              /* endIJK_vec */ {{1,1,1}},
                              /* lgr_name_vec */ {"LGR1"});

    const auto& refinedGridData = *grid.currentData()[1];
    const auto& parentGridData = *grid.currentData()[0];

    const auto parentElem = Dune::cpgrid::Entity<0>(parentGridData, 0, true);

    for (const auto& refinedElem : Dune::elements(grid.levelGridView(1))) {
        printIfNewVertex(refinedGridData, refinedElem, parentGridData, parentElem);
    }
}
