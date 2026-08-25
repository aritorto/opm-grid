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


// Level zero grid dims = 2x2x1
//
// Cell 0 (x=0..6,   y=0..6)
// bottom face corners (0,0,0), (6,0,0), (0,6,1), (6,6,1)
//    top face corners (0,0,8), (6,0,8), (0,6,9), (6,6,9)
//
// Cell 1 (x=6..12,  y=0..6)
// bottom face corners (6,0,1), (12,0,1), (6,6,2), (12,6,2)
//    top face corners (6,0,9), (12,0,9), (6,6,10), (12,6,10)
//
// Cell 2 (x=0..6,   y=6..12)
// bottom face corners (0,6,0), (6,6,0), (0,12,1), (6,12,1)
//    top face corners (0,6,8), (6,6,8), (0,12,9), (6,12,9)
//
// Cell 3 (x=6..12,  y=6..12)
// bottom face corners (6,6,1), (12,6,1), (6,12,2), (12,12,2)
//    top face corners (6,6,9), (12,6,9), (6,12,10), (12,12,10)

const std::string deckTwoColumnsGrid =
    R"(RUNSPEC
DIMENS
 2 1 2 /

GRID

COORD
 0  0 0     0  0 17
 6  0 0     6  0 17
12  0 0    12  0 17

 0  6 0     0  6 17
 6  6 0     6  6 17
12  6 0    12  6 17
/

ZCORN
-- cell0   cell1    cell0   cell1
-- p0 p1   p0 p1    p2 p3   p2 p3     
    0  0    1  1     0  0    1  1    -- bottom layer 1

-- cell2   cell3    cell2   cell3
-- p4 p5   p4 p5    p6 p7   p6 p7 
    8  8    9  9     8  8    9  9    --    top layer 1

-- cell2   cell3    cell2   cell3
-- p0 p1   p0 p1    p2 p3   p2 p3
    8  8    9  9     8  8    9  9    -- bottom layer 2
-- cell2   cell3    cell2   cell3
-- p4 p5   p4 p5    p6 p7   p6 p7 
   16 16   17 17    16 16   17 17    --    top layer 2
/

ACTNUM
4*1
/

PORO
4*0.15
/
)";


BOOST_AUTO_TEST_CASE(faultBetweenLgrAndGlobalCells)
{
    /* Add in deck:

       CARFIN
       -- name  I1 I2 J1 J2 K1 K2 NX NY NZ
       'LGR1'    1  2  1  1  1  1  2  2  1
       ENDFIN

       This together with deck deckTwoColumnsGrid
       translates into

       startIJK      = {0,0,0}
       endIJK        = {2,1,1}
       cells_per_dim = {1,2,1}
       parent-cell-block-dimensions = 2x1x1
       LGR1 dimensions = {NX, NY, NZ} = {2,2,1}

       To be able to invoke
       processEclipseFormat(....);
       to create a CpGridData object containing the geoemtries representing LGR1,
       I need:
       1. eclipseLgr1Grid,
       2. lgr1EclState

       QUESTION 1: Can I avoid creating the EclipseState object?
                   What's the simplest EclipseState object I could use here?
      
       Assuming that we can get from the modified input deck-created grid
       the COORD, ZCORN, and ACTNUM values for the CARFIN block, we might
       be able to create an EclipseGrid representing only the LGR1, that
       we can process into a CpGridData object.

       In this simple case, for LGR1 defined as above in the CARFIN block:
       
      const std::vector<double> lgr1_coord = {
              0, 0, 0,     0, 0, 9,
              6, 0, 0,     6, 0, 9,
             12, 0, 0,    12, 0, 9,
              0, 3, 0,     0, 3, 9,
              6, 3,0,     6, 3, 9,
             12, 3, 0,    12, 3, 9,
              0, 6, 0,     0, 6, 9,
              6, 6, 0,     6, 6, 9,
             12, 6, 0,    12, 6, 9
        };

    const std::vector<double> lgr1_zcorn = {
        8, 8, 9, 9,        8,8, 9, 9,        8,8,9,9,      8,8, 9, 9,
        0, 0, 1, 1,        0, 0, 1, 1,       0,0, 1, 1,      0,0,1,1
    };

    I couldn't make this construction work:

    Opm::EclipseGrid lgr1Grid({2,2,1}, lgr1_coord, lgr1_zcorn);

    Alternatively, I created an lgr1DeckString (see below).
    There should be an easier way to do this, instead of going through a deck. 
     
     */

    // -------- Potential steps for refinement from deck ---------
    // 
    // 1. From input grid *.DATA file create a CpGrid
    //
    // 2. From input grid *.DATA collect, for each CARFIN block,
    //    the necessary information to create
    //    - lgrEclipseGrid
    //    - lgrEclipseState
    //
    // 3. Create, for each {lgrEclipGrid, lgrEclipseState} "pair"
    //    a CpGridData object where processEclipseFormat(...)
    //    can be invoked.
    //
    // 4. Even though step 3. provides a collection of CpGridData
    //    objects that each of them represents a different LGR
    //    defined via CARFIN in the input deck, all the relationships
    //    parent-to-children need to be defined.
    //
    //    For example:
    //    level_to_leaf_cells_       container mapping {level, level cell index} -> equivalent leaf cell index
    //    parent_to_children_cells_  ...
    //    cell_to_idxInParentCell_   ...
    //
    //    QUESTION: Is some of this information already avalible in opm-common?
    //              Maybe in some kind of equivalent "maps"
    //
    // 5. Push all the level grids (CpGridData objects, "completly populated")
    //    into the CpGrid grid. They become new entries in grid.current_data_
    //
    // 6. Build the leaf grid discarting repeated entities and creating new intersections,
    //    if needed, at boundary of two cells of different levels.
    //    It can be:
    //    - coarse cell (level zero grid)  and refined cell (LGR cell)
    //    - refined cell (LGR1 with nx1, ny1, nz1 subdivisions) and
    //      refined cell from other level  (LGR2 with nx2, ny2, nz2 subdivisions)
    //
    //
    // How to deal with different levels sharing boundary (an idea):
    //
    // For each coarse intersection involved in one of these situations:
    // - situation coarse-refined:   coarse cell (level zero grid)  and refined cell (LGR cell)
    // - situtation refined-refined: refined cell (LGR1 with nx1, ny1, nz1 subdivisions) and
    //                               refined cell from other level  (LGR2 with nx2, ny2, nz2 subdivisions)
    // create an auxiliary CpGridData aux_grid with the 2 cells attached to such face, repeat the cheapest
    // version of the "deck-refinement" process, to obtain again through
    // aux_grid.processEclipseFormat(...) the correct leaf grid new intersections.
    // For that, we need (as before) COORD, ZCORN, and ACTNUM, ... stuff to be able to create
    // aux_eclipseGrid and aux_eclipseState, to be able to call processEclipseFormat(...).


    
    // Incomplete example to illustrate steps 1.-6.
    // -----------------------------------------------------------------------------------
    // 1. From input grid *.DATA file create a CpGrid
    Dune::CpGrid grid;
    Opm::createGridFromDeckString(grid, deckTwoColumnsGrid);

    auto& data = grid.currentData();

    // 2. From input grid *.DATA collect, for each CARFIN block,
    //    the necessary information to create
    //    - lgrEclipseGrid
    //    - lgrEclipseState
    // (Pendent, not done here; perhaps already existing in opm-common or Halvor's work)


    // 3. Create, for each {lgrEclipGrid, lgrEclipseState} "pair"
    //    a CpGridData object where processEclipseFormat(...)
    //    can be invoked.
    //
    // Unfortunately, for now, I have to create a deck.
    // (Check better options from opm-common or Halvor's work)
    // Auxiliary deck string - hopefully there is an easier way to do this
    const std::string lgr1AuxDeckString =
        R"(RUNSPEC
DIMENS
 2 2 1 /

GRID

COORD
 0 0 0     0 0 9
 6 0 0     6 0 9
12 0 0    12 0 9

 0 3 0     0 3 9
 6 3 0     6 3 9
12 3 0    12 3 9

 0 6 0     0 6 9
 6 6 0     6 6 9
12 6 0    12 6 9
/

ZCORN
0 0 1 1  0 0 1 1  0 0 1 1  0 0 1 1
8 8 9 9  8 8 9 9  8 8 9 9  8 8 9 9 
/

ACTNUM
4*1
/

PORO
4*0.15
/
)";

       Opm::Parser parser;  
       const auto lgr1Deck = parser.parseString(lgr1AuxDeckString);
       Opm::EclipseState lgr1EclState(lgr1Deck);
       Opm::EclipseGrid eclipseLgr1Grid = lgr1EclState.getInputGrid();
    
       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& lgr1Data = data;
       std::shared_ptr<Dune::cpgrid::CpGridData> lgr1Grid_ptr = std::make_shared<Dune::cpgrid::CpGridData>(lgr1Data); // ccobj_
       auto& lgr1Grid = *lgr1Grid_ptr;

       lgr1Grid.processEclipseFormat(&eclipseLgr1Grid, &lgr1EclState, false, false, false);


       // 4. Even though step 3. provides a collection of CpGridData
       //    objects that each of them represents a different LGR
       //    defined via CARFIN in the input deck, all the relationships
       //    parent-to-children need to be defined.
       //
       // (Pendent, not done here; perhaps already existing in opm-common or Halvor's work)


       
       // 5. Push all the level grids (CpGridData objects, "completly populated")
       //    into the CpGrid grid. They become new entries in grid.current_data_
       data.push_back(lgr1Grid_ptr);


       // 6. Build the leaf grid discarting repeated entities and creating new intersections,
       //    if needed, at boundary of two cells of different levels.
       //    It can be:
       //    - coarse cell (level zero grid)  and refined cell (LGR cell)
       //    - refined cell (LGR1 with nx1, ny1, nz1 subdivisions) and
       //      refined cell from other level  (LGR2 with nx2, ny2, nz2 subdivisions)
       std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& leafData = data;
       auto leafGrid_ptr = std::make_shared<Dune::cpgrid::CpGridData>(leafData);
       auto& leafGrid = *leafGrid_ptr;
       // (Pendent, not done here; perhaps already existing in opm-common or Halvor's work)
       //
       // -----------------------------------------------------------------------------------





       
       // These are the expected faces in LGR1
       const std::vector<std::vector<Dune::FieldVector<double,3>>> expectedLgr1Faces = {
           {{0.,0.,0.}, {0.,3.,0.}, {0.,3.,8.}, {0.,0.,8.}},      // I_FACE   x = 0,  face 0
           {{6.,0.,0.}, {6.,3.,0.}, {6.,3.,1.}, {6.,0.,1.}},      // I_FACE   x = 6,  face 1
           {{6.,0.,1.}, {6.,3.,1.}, {6.,3.,8.}, {6.,0.,8.}},      // I_FACE   x = 6,  face 2
           {{6.,0.,8.}, {6.,3.,8.}, {6.,3.,9.}, {6.,0.,9.}},      // I_FACE   x = 6,  face 3
           {{12.,0.,1.}, {12.,3.,1.}, {12.,3.,9.}, {12.,0.,9.}},  // I_FACE   x = 12, face 4
           {{0.,3.,0.}, {0.,6.,0.}, {0.,6.,8.}, {0.,3.,8.}},      // I_FACE   x = 0,  face 5
           {{6.,3.,0.}, {6.,6.,0.}, {6.,6.,1.}, {6.,3.,1.}},      // I_FACE   x = 6,  face 6
           {{6.,3.,1.}, {6.,6.,1.}, {6.,6.,8.}, {6.,3.,8.}},      // I_FACE   x = 6,  face 7
           {{6.,3.,8.}, {6.,6.,8.}, {6.,6.,9.}, {6.,3.,9.}},      // I_FACE   x = 6,  face 8
           {{12.,3.,1.}, {12.,6.,1.}, {12.,6.,9.}, {12.,3.,9.}},  // I_FACE   x = 12, face 9
           {{ 6.,0.,0.}, {0.,0.,0.}, {0.,0.,8.}, { 6.,0.,8.}},    // J_FACE   y = 0,  face 10
           {{12.,0.,1.}, {6.,0.,1.}, {6.,0.,9.}, {12.,0.,9.}},    // J_FACE   y = 0,  face 11
           {{ 6.,3.,0.}, {0.,3.,0.}, {0.,3.,8.}, { 6.,3.,8.}},    // J_FACE   y = 3,  face 12        
           {{12.,3.,1.}, {6.,3.,1.}, {6.,3.,9.}, {12.,3.,9.}},    // J_FACE   y = 3,  face 13
           {{ 6.,6.,0.}, {0.,6.,0.}, {0.,6.,8.}, { 6.,6.,8.}},    // J_FACE   y = 6,  face 14
           {{12.,6.,1.}, {6.,6.,1.}, {6.,6.,9.}, {12.,6.,9.}},    // J_FACE   y = 6,  face 15
           {{0.,0.,0.}, { 6.,0.,0.}, { 6.,3.,0.}, {0.,3.,0.}},    // K_FACE   z = 0,  face 16
           {{0.,0.,8.}, { 6.,0.,8.}, { 6.,3.,8.}, {0.,3.,8.}},    // K_FACE   z = 8,  face 17
           {{6.,0.,1.}, {12.,0.,1.}, {12.,3.,1.}, {6.,3.,1.}},    // K_FACE   z = 1,  face 18
           {{6.,0.,9.}, {12.,0.,9.}, {12.,3.,9.}, {6.,3.,9.}},    // K_FACE   z = 9,  face 19
           {{0.,3.,0.}, { 6.,3.,0.}, { 6.,6.,0.}, {0.,6.,0.}},    // K_FACE   z = 0,  face 20
           {{0.,3.,8.}, { 6.,3.,8.}, { 6.,6.,8.}, {0.,6.,8.}},    // K_FACE   z = 8,  face 21
           {{6.,3.,1.}, {12.,3.,1.}, {12.,6.,1.}, {6.,6.,1.}},    // K_FACE   z = 1,  face 22
           {{6.,3.,9.}, {12.,3.,9.}, {12.,6.,9.}, {6.,6.,9.}},    // K_FACE   z = 9,  face 23
       };

       Opm::checkFaces(lgr1Grid, expectedLgr1Faces);

       std::cout<< lgr1Grid.size(0) << " lgr1 cells "<< std::endl;
       std::cout<< lgr1Grid.numFaces() << " lgr1 faces "<< std::endl;
       std::cout<< lgr1Grid.size(3) << " lgr1 corners "<< std::endl;  
    
}


BOOST_AUTO_TEST_CASE(faultBetweenTwoColumnsSameLgr)
{

    /* Add in deck:

       CARFIN
       -- name  I1 I2 J1 J2 K1 K2 NX NY NZ
       'LGR1'    1  2  1  1  1  2  2  2  2
       ENDFIN
       
     */
   
}

BOOST_AUTO_TEST_CASE(faultBetweenTwoColumnsDifferentVerticalLgrs)
{
    /* Add in deck:

       CARFIN
       -- name  I1 I2 J1 J2 K1 K2 NX NY NZ
       'LGR1'    1  1  1  1  1  2  1  2  2
       ENDFIN

        CARFIN
       -- name  I1 I2 J1 J2 K1 K2 NX NY NZ
       'LGR2'    2  2  1  1  1  2  1  2  2
       ENDFIN
       
     */

}

BOOST_AUTO_TEST_CASE(faultBetweenTwoColumnsDifferentHorizontalLgrs)
{
     /* Add in deck:

       CARFIN
       -- name  I1 I2 J1 J2 K1 K2 NX NY NZ
       'LGR1'    1  2  1  1  1  1  1  2  1
       ENDFIN

        CARFIN
       -- name  I1 I2 J1 J2 K1 K2 NX NY NZ
       'LGR2'    1  2  1  1  2  2  1  2  1
       ENDFIN
       
     */

}

