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

#define BOOST_TEST_MODULE LgrWithFaultsTests
#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/LgrFaultHelpers.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <tests/cpgrid/lgr/LgrChecks.hpp>

#include <string>



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


BOOST_AUTO_TEST_CASE(skewPillarsSimpleGrid)
{

    const std::string deckString =
        R"(RUNSPEC
DIMENS
 2 1 1 /

GRID

COORD
 0 0 0     0 0 9
 6 0 0     6 3 9
12 0 0    12 0 9

 0 6 0    0 6 9
 6 6 0    6 6 9
12 6 0   12 6 9
/

ZCORN
0 0 1 1  0 0 1 1
8 8.2 9 9  7.1 5 9 9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid,
                                  deckString);

    // Element faces
    //
    // I- Face index 0:
    //
    // {0, 0, 8} 
    //     |      \
    //     |         {0, 6, 7.1}
    //     |             |
    //     |             |
    // {0, 0, 0} --- {0, 6, 0}
    //
    //
    //          {6, 2.73333, 8.2}
    //              /           \
    //             /             \
    //            /               \
    //           /                  {6, 6, 5}
    //          /                    |
    //         /  I+ Face index 2    |
    //        /                      |
    //   {6, 0.33333, 1}----------  {6, 6, 1}
    //       /                       |
    //      /  I+ Face index 1       |
    //     /                         |
    // {6, 0, 0} -----------------  {6, 6, 0}
    

    const auto parentCell = Dune::cpgrid::Entity<0>(grid.currentLeafData(), 0, true);
    
    for (const auto& intersection : Dune::intersections(grid.leafGridView(), parentCell)) {
        
        const auto& faceToPoint = grid.currentLeafData().faceToPoint(intersection.id());
        const auto faceTag =  grid.currentLeafData().faceTag(intersection.id());
       
        if (faceTag == 0) {// I face
        std::cout<< "Face index: " << intersection.id() << std::endl;
        for (const auto& point : faceToPoint) {
            const auto v = Dune::cpgrid::Entity<3>( grid.currentLeafData(), point, true).geometry().center();
            std::cout<< v[0] << " " << v[1] << " " << v[2] << std::endl;
        }
        std::cout<<std::endl;
        }
    }
    
    grid.addLgrsUpdateLeafView({{1,2,1}}, // cells_per_dim
                               {{0,0,0}}, // startIJK
                               {{1,1,1}}, // endIJK
                               {"LGR1"}); // lgr name

    std::cout<< "Faces after refinement element zero " <<std::endl;
    std::cout<<std::endl;

    for (const auto& element : Dune::elements(grid.levelGridView(1))) {
        if (Opm::Lgr::isAtGridBoundary(*grid.currentData()[1], element)) {

            for (const auto& intersection : Dune::intersections(grid.levelGridView(1), element)) {

                if (!intersection.neighbor()) {
        
                const auto& faceToPoint = grid.currentData()[1]->faceToPoint(intersection.id());
                const auto faceTag =  grid.currentData()[1]->faceTag(intersection.id());
       
                if (faceTag == 0) {// I face
                    std::cout<< "Face index: " << intersection.id() << std::endl;
                    for (const auto& point : faceToPoint) {
                        const auto v = Dune::cpgrid::Entity<3>( *grid.currentData()[1], point, true).geometry().center();
                        std::cout<< v[0] << " " << v[1] << " " << v[2] << std::endl;
                    }
                    std::cout<<std::endl;
                }
                }
            }
        }
    }



     bool isInteriorInA, isInteriorInB;
    
     const auto seg =  Opm::Lgr::computeSegmentIntersection(/* startA */ {6., 3., 0.}, /* endA */ {6., 4.36667, 6.6},
                                                           /* startB */ {6., 0.33333, 1.}, /* endB */ {6., 6., 1.},
                                                           isInteriorInA,
                                                           isInteriorInB);
    if (seg.has_value()) {
        const auto [p,q] = seg.value();
        std::cout<< p[0] << " " << p[1] << " " << p[2] << std::endl;
        std::cout<< q[0] << " " << q[1] << " " << q[2] << std::endl;
    }

}
