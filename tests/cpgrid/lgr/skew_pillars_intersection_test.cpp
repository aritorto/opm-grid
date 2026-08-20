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
#include <opm/grid/cpgrid/CpGridUtilities.hpp>
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

std::array<Dune::FieldVector<double,3>,8> unitCube = {{
        {0.,0.,0.}, {1.,0.,0.}, {0.,1.,0.}, {1.,1.,0.},
        {0.,0.,1.}, {1.,0.,1.}, {0.,1.,1.}, {1.,1.,1.}
    }};

void checkUnitCubeCorners(const Dune::cpgrid::CpGridData& grid)
{
    for (int elemIdx = 0; elemIdx < grid.size(0); ++elemIdx) {
        
        const auto cellGeom = Dune::cpgrid::Entity<0>( grid, elemIdx, true).geometry();
        const auto& cellToPoint = grid.cellToPoint(elemIdx);
        
        BOOST_CHECK_EQUAL(cellToPoint.size(), 8);
        
        for (int corner = 0; corner < 8; ++corner) {
            const auto vertex = Dune::cpgrid::Entity<3>( grid, cellToPoint[corner], true).geometry().center();
            const auto localVtx = cellGeom.local(vertex);
            BOOST_CHECK( Opm::Lgr::areClose(cellGeom.local(vertex), unitCube[corner]) );

            BOOST_CHECK( Opm::Lgr::areClose(cellGeom.global(localVtx), vertex) );
        }
    }
}


BOOST_AUTO_TEST_CASE(skewPillarsSimpleGrid)
{
    /* DIMENS
       2 1 1 /    nx ny nz

       COORD      (nx+1)*(ny+1) pillars
 0 0 0     0 0 9 bottom - top pillar 0
 6 0 0     6 3 9 bottom - top pillar 1
12 0 0    12 0 9 bottom - top pillar 2

 0 6 0    0 6 9  bottom - top pillar 3
 6 6 0    6 6 9  bottom - top pillar 4
12 6 0   12 6 9  bottom - top pillar 5

     */
    
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
0   0  1  1    0  0  1  1
8 8.2  9  9  7.1  5  9  9
/

ACTNUM
2*1
/

PORO
2*0.15
/
)";

    Opm::Parser parser;
    const auto deck = parser.parseString(deckString);
    Opm::EclipseState ecl_state(deck);
    Opm::EclipseGrid eclipse_grid = ecl_state.getInputGrid();

    // I. Collect input info about the pillars
    const auto& input_coord = eclipse_grid.getCOORD();
    std::vector<std::pair<Dune::FieldVector<double,3>, Dune::FieldVector<double,3>>> pillarsBottomTop{};
    pillarsBottomTop.resize(6); // 6 = (nx+1)*(ny+1)

    int vertices_count = input_coord.size()/6;
    
    for (int i = 0; i < vertices_count; ++i)
    {
        pillarsBottomTop[i] = std::make_pair<Dune::FieldVector<double,3>, Dune::FieldVector<double,3>>({input_coord[6*i], input_coord[6*i+1], input_coord[6*i+2]},
                                          {input_coord[6*i+3], input_coord[6*i+4], input_coord[6*i+5]});
    }
    for (const auto& [b, t] : pillarsBottomTop) {
        std::cout<< b[0] << " " << b[1] << " " << b[2] << " bottom " <<std::endl;
        std::cout<< t[0] << " " << t[1] << " " << t[2] << " top " <<std::endl;

        std::cout<<std::endl;
    }
    // II. Create a Geometry "cell" for each "column" (4 pillars: (i,j), (i+1,j), (i, j+1), (i+1,j+1))
    //     Total amount of column-cells: nx*ny
    // pillarsBottomTop = pillar_0,       ..., pillar_nx       (j=0)
    //                    pillar_(nx+1), ...., pillar_2*nx (j=1)
    // ...
    //                    pillar_((nx+1)*ny)
    int nx = eclipse_grid.getNX();
    int ny = eclipse_grid.getNY();

    static constexpr std::array<int,8> corner_indices = {0,1,2,3,4,5,6,7};

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // Calculate center and volume for  "column" geometry (4 pillars: (i,j), (i+1,j), (i, j+1), (i+1,j+1))
            const int p0 = (j*(nx + 1)) + i;
            const int p1 = (j*(nx + 1)) + (i + 1);
            const int p2 = ((j + 1)*(nx + 1)) + i;
            const int p3 = ((j + 1)*(nx + 1)) + (i +1);

            const std::array<Dune::FieldVector<double,3>,8> corners = {{pillarsBottomTop[p0].first, pillarsBottomTop[p0].second,
                                                                            pillarsBottomTop[p1].first, pillarsBottomTop[p1].second,
                                                                            pillarsBottomTop[p2].first, pillarsBottomTop[p2].second,
                                                                            pillarsBottomTop[p3].first, pillarsBottomTop[p3].second}};
            
            const auto [center, volume] = Opm::computeCenterAndVolume(corners);

            
            std::cout<< center[0] << " " << center[1] << " " << center[2] << " center, vol : " << volume <<std::endl;

            /*     auto pillarCell_corners = std::make_shared<EntityVariable<cpgrid::Geometry<0, 3>, 3>>();
        EntityVariableBase<cpgrid::Geometry<0, 3>>& mutable_in_father_reference_elem_corners = *in_father_reference_elem_corners;
        // Assign the corners. Make use of the fact that pointers behave like iterators.
        mutable_in_father_reference_elem_corners.assign(corners_in_father_reference_elem_temp,
        corners_in_father_reference_elem_temp + 8);*/

            const auto pillarCell = Dune::cpgrid::Geometry<3,3>(center, volume, corners, corner_indices.data());
        }
    }
    

   

    Dune::CpGrid grid;
    grid.processEclipseFormat(&eclipse_grid, &ecl_state, false, false, false);
    /*  Opm::createGridFromDeckString(grid,
        deckString);*/

    /* Element faces
    //
    // I- Face index 0:
    //
    // {0, 0, 8} 
    ///    |      \
    //     |         {0, 6, 7.1}
    //     |             |
    //     |             |
    // {0, 0, 0} --- {0, 6, 0}
    //
    //
    //          {6, 2.73333, 8.2}
    ///             /           \
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
    */


    // 1. Extract the corner computation from Geomtry header of a single-cell-refinement
    // 2. Create refined-cell-pillars with extra vertices of parent faces
    // 3. Rearrange the info from step 2 to have a coord vector 
    // 4. zcorn can be taken from Opm::Lgr::lgrCOORDandZCORN(...) [here coord is wrong for parent cell with >6 intersections, but zcorn should be okay]
    // 5. Create an actnum for this refinement
    // 6. Create an eclipGrid(dims, coord, zcorn, actnum)
    // 7. Create a CpGridData with processEclipseGrid();

    checkUnitCubeCorners(grid.currentLeafData());
    
    const auto parentCell = Dune::cpgrid::Entity<0>(grid.currentLeafData(), 0, true);
    
    auto input_vertices = Opm::computeBasicRefinedCorners(parentCell,
                                                     {1,2,1}, // nxnynz 
                                                     {1.},{.5, .5},{1.}); //  widthsX,lengthsY,heightsZ

    Opm::addAllParentCellFaceVertices(grid.currentLeafData(),
                                 parentCell,
                                 input_vertices);
    


    
    const auto& cellPillars = Opm::Lgr::extendCellPillars(grid.currentLeafData(), parentCell.index());
    for (const auto& pillar :cellPillars)
    {
        for (const auto& p : pillar) {
            std::cout<< p << " pillar " <<std::endl;
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    
    /* const auto& extendedCellToPoint = Opm::Lgr::buildExtendedCellPointVertexMap(grid.currentLeafData(), parentCell.index());
    for (const auto& [v, idx] : extendedCellToPoint) {
        std::cout<< v[0] << " " << v[1] << " " << v[2] << " has index: " << idx<<  std::endl;
        }*/


    std::cout<< " Vertices input " << std::endl;
    for (const auto& v : input_vertices) {
         std::cout<< v[0] << " " << v[1] << " " << v[2] << std::endl;
    }


    

    /* const auto p = parentCellGeom.local({6, 2.73333, 8.2}); //{6., 0.3333, 1.} ); // {6., 4.36667, 6.6})
    std::cout<< p[0] << " " << p[1] << " " << p[2] << " local!" << std::endl;

     const auto g = parentCellGeom.global({1, 0, 1}); //{6., 0.3333, 1.} ); // {6., 4.36667, 6.6})
    std::cout<< g[0] << " " << g[1] << " " << g[2] << " global!" << std::endl;

     const auto l = parentCellGeom.local(g); //{6., 0.3333, 1.} ); // {6., 4.36667, 6.6})
    std::cout<< l[0] << " " << l[1] << " " << l[2] << " local!" << std::endl;
 
    std::set<Dune::FieldVector<double,3>, Opm::Lgr::FieldVectorLess> output_vertices{};
    
    for (int i = 0; i < cellRefGrid.size(3); ++i) {
        const auto v = Dune::cpgrid::Entity<3>( cellRefGrid, i, true).geometry().center();
        output_vertices.insert(v);
    }
    std::cout<<std::endl;

     std::cout<< " Vertices output " << std::endl;
    for (const auto& v : output_vertices) {
         std::cout<< v[0] << " " << v[1] << " " << v[2] << std::endl;
    }

    const auto& [coord, zcorn] = Opm::lgrCOORDandZCORN(cellRefGrid,  {1,2,1});
    for (const auto& c : coord)
    {
        std::cout<< c << std::endl;
    }
    */
}
