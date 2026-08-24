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

#define BOOST_TEST_MODULE LocalCoordinatesTests
#include <boost/test/unit_test.hpp>

#include <dune/common/fvector.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Geometry.hpp>

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

std::array<Dune::FieldVector<double,3>,8> unitCube = {{
        {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {1.,1.,0.},
        {0.,0.,1.}, {1.,0.,1.}, {0.,1.,1.}, {1.,1.,1.}
    }};


bool areClose(const Dune::FieldVector<double,3>& v,
              const Dune::FieldVector<double,3>& w)
{
    return (std::abs(v[0] - w[0]) < 1e-12) && (std::abs(v[1] - w[1]) < 1e-12) && (std::abs(v[2] - w[2])< 1e-12);
}

void checkUnitCubeCorners(const Dune::cpgrid::CpGridData& grid,
                          const std::array<Dune::FieldVector<double,3>,8>& expectedLocalCorners)
{
    for (int elemIdx = 0; elemIdx < grid.size(0); ++elemIdx) {

        const auto cellGeom = Dune::cpgrid::Entity<0>( grid, elemIdx, true).geometry();
        const auto& cellToPoint = grid.cellToPoint(elemIdx);

        BOOST_CHECK_EQUAL(cellToPoint.size(), 8);

        for (int corner = 0; corner < 8; ++corner) {
            const auto vertex = Dune::cpgrid::Entity<3>( grid, cellToPoint[corner], true).geometry().center();
            const auto localVtx = cellGeom.local(vertex);
            
            BOOST_CHECK( areClose(cellGeom.local(vertex), expectedLocalCorners[corner]) );

            BOOST_CHECK( areClose(cellGeom.global(localVtx), vertex) );
        }
    }
}

BOOST_AUTO_TEST_CASE(cellWithFewerThanEightCornersDoesNotMapToUnitCube)
{
    const std::string deckString = R"(
RUNSPEC

DIMENS
 1  1  1/

GRID

-- COORD: 4 pillars × 6 values (x1 y1 z1  x2 y2 z2)
COORD
0 0 0  0 0 1   -- Pillar 1: bottom at (0,0,0), top at (0,0,1)
1 0 0  1 0 0   -- Pillar 2: bottom at (1,0,0) = same at top
0 1 0  0 1 1   -- Pillar 3: bottom at (0,1,0), top at (0,1,1)
1 1 0  1 1 0   -- Pillar 4: bottom at (1,1,0) = same at top
/

-- ZCORN: eight Z values: top 4, then bottom 4
ZCORN
0 0 0 0 -- bottom end of the pillars
1 0 1 0 -- top end of the pillars
/

ACTNUM
1*1
/

PORO
1*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    /** Global cell [collapsed] corners = { {0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
                                            {0,0,1}, {1,0,0}, {0,1,1}, {1,1,0} }
        Local cell [collapsed] corners  = { {0,0,0}, {1,0,.5}, {0,1,0}, {1,1,.5},
                                            {0,0,1}, {1,0,.5}, {0,1,1}, {1,1,.5} }
    */
    std::array<Dune::FieldVector<double,3>,8> expected = {{ {0,0,0}, {1,0,.5}, {0,1,0}, {1,1,.5},
                                                            {0,0,1}, {1,0,.5}, {0,1,1}, {1,1,.5} }};
    
    checkUnitCubeCorners(grid.currentLeafData(), expected); // Collaped corners also in reference element
}



BOOST_AUTO_TEST_CASE(skewPillarsSimpleGrid)
{
    /*  For a grid with dimensions 2x1x1, in ZCORN: z-coord of
        corner0-cell0 corner1-cell0   corner0-cell1 corner1-cell1   corner2-cell0 corner3-cell0   corner2-cell1 corner3-cell1 (bottom layer 1)
        corner4-cell0 corner5-cell0   corner4-cell1 corner5-cell1   corner6-cell0 corner7-cell0   corner6-cell1 corner7-cell1 (   top layer 1)
    */

    const std::string deckString =
        R"(RUNSPEC
DIMENS
 2 1 1 /  -- nx ny nz

GRID

COORD             -- (nx+1)*(ny+1) = (2+1)*(1+1) = 6 pillars
 0 0 0     0 0 9  -- bottom - top pillar 0
 6 0 0     6 3 9  -- bottom - top pillar 1
12 0 0    12 0 9  -- bottom - top pillar 2

 0 6 0    0 6 9  --  bottom - top pillar 3
 6 6 0    6 6 9  --  bottom - top pillar 4
12 6 0   12 6 9  --  bottom - top pillar 5
/

ZCORN                     -- 8*nx*ny values = 8*2*1 = 16
0   0  1  1    0  0  1  1 -- bottom layer 1
8 8.2  9  9  7.1  5  9  9 -- top layer 1
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckString);

    checkUnitCubeCorners(grid.currentLeafData(), unitCube);
}
