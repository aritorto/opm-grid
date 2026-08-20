/*
  Copyright 2025 Equinor ASA.

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

#include <config.h>

#include <opm/grid/cpgrid/CpGridUtilities.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>
#include <opm/grid/cpgrid/LgrFaultHelpers.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm
{

std::pair<std::unordered_map<int, int>, std::vector<std::array<int, 3>>>
lgrIJK(const Dune::CpGrid& grid, const std::string& lgr_name)
{
    // Check if lgr_name exists in lgr_names_
    const auto& lgr_names = grid.getLgrNameToLevel();
    auto it = lgr_names.find(lgr_name);
    if (it == lgr_names.end()) {
        OPM_THROW(std::runtime_error, "LGR name not found: " + lgr_name);
    }

    const auto level = it->second;
    const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapper(grid);
    const auto& levelView = grid.levelGridView(level);
    const auto numCells = levelView.size(0);

    std::vector<std::array<int, 3>> lgrIJK(numCells);
    std::unordered_map<int, int> lgrCartesianIdxToCellIdx;
    lgrCartesianIdxToCellIdx.reserve(numCells);

    // Iterate over (active) elements in the grid and populate the structures
    for (const auto& element : Dune::elements(grid.levelGridView(level))) {
        std::array<int, 3> ijk;
        levelCartMapper.cartesianCoordinate(element.index(), ijk, level);

        const int cellIndex = element.index();
        const int cartesianIdx = element.getLevelCartesianIdx();

        lgrIJK[cellIndex] = ijk;
        lgrCartesianIdxToCellIdx[cartesianIdx] = cellIndex;
    }

    return std::make_pair(lgrCartesianIdxToCellIdx, lgrIJK);
}

std::pair<std::vector<double>, std::vector<double>>
lgrCOORDandZCORN(const Dune::CpGrid& grid,
                 int level,
                 const std::unordered_map<int, int>&  lgrCartesianIdxToCellIdx,
                 const std::vector<std::array<int, 3>>& lgrIJK)
{
    const auto& levelGrid = *(grid.currentData()[level]);

    // Check not all cells are inactive
    const auto numCells = levelGrid.size(0);
    if (numCells == 0) {
        OPM_THROW(std::logic_error, "LGR in level " + std::to_string(level) + " has no active cells.\n");
    }

    // LGR dimensions
    const auto& lgr_dim = grid.currentData()[level]->logicalCartesianSize();
    const int nx = lgr_dim[0];
    const int ny = lgr_dim[1];
    const int nz = lgr_dim[2];

    // Initialize all pillars as inactive (setting COORD values to std::numeric_limits<double>::max()).
    std::vector<double> lgrCOORD(6*(nx+1)*(ny+1), std::numeric_limits<double>::max());

    // Initialize all ZCORN as inactive (setting values to std::numeric_limits<double>::max()).
    std::vector<double> lgrZCORN(8*nx*ny*nz, std::numeric_limits<double>::max());

    // Map to determine min and max k per cell column (i, j) (min/max_k = 0, ..., nz-1).
    // Initialized as {nz, -1} to detect inactive cell columns.
    std::vector<std::array<int,2>> minMaxPerCellPillar(nx*ny, {nz, -1});

    for (const auto& ijk : lgrIJK) {

        // Compute the bottom and top k per cell pillar (i, j).
        int cell_pillar_idx = ijk[1] * nx + ijk[0];
        auto& minMax = minMaxPerCellPillar[cell_pillar_idx];

        minMax[0] = std::min(ijk[2], minMax[0]);
        minMax[1] = std::max(ijk[2], minMax[1]);
    }

    for (const auto& elem : elements(grid.levelGridView(level))) {
        const auto& elemIJK = lgrIJK[elem.index()];

        // For a grid with nz layers, ZCORN values are ordered:
        //
        //      top layer nz-1
        //   bottom layer nz-1
        //      top layer nz-2
        //   bottom layer nz-2
        // ...
        //      top layer 1
        //   bottom layer 1
        //      top layer 0
        //   bottom layer 0

        int zcorn_top_00_idx = ((nz-1-elemIJK[2])*8*nx*ny) + (elemIJK[1]*4*nx) + (2*elemIJK[0]); // assoc. w. elem corner 4

        // Bottom indices
        int zcorn_top_10_idx = zcorn_top_00_idx + 1;  // assoc. w. elem corner 5
        int zcorn_top_01_idx = zcorn_top_00_idx + (2*nx);  // assoc. w. elem corner 6
        int zcorn_top_11_idx = zcorn_top_01_idx + 1; // assoc. w. elem corner 7

        // Top indices
        int zcorn_bottom_00_idx = zcorn_top_00_idx + (4*nx*ny); // assoc. w. elem corner 0
        int zcorn_bottom_10_idx = zcorn_bottom_00_idx + 1;  // assoc. w. elem corner 1
        int zcorn_bottom_01_idx = zcorn_bottom_00_idx + (2*nx); // assoc. w. elem corner 2
        int zcorn_bottom_11_idx = zcorn_bottom_01_idx + 1;  // assoc. w. elem corner

        // Note: zcorn_idx + 1 moves to the next position along the x-axis (i+1, j, k)
        //       zcorn_idx + (2*nx) moves to the next position along the y-axis (i, j+1, k)
        //       zcorn_idx + (4*nx*ny) moves to the next position along the z-axis (i,j, k+1)

        // Assign ZCORN values
        lgrZCORN[zcorn_top_00_idx] = elem.subEntity<3>(4).geometry().center()[2];
        lgrZCORN[zcorn_top_10_idx] = elem.subEntity<3>(5).geometry().center()[2];
        lgrZCORN[zcorn_top_01_idx] = elem.subEntity<3>(6).geometry().center()[2];
        lgrZCORN[zcorn_top_11_idx] = elem.subEntity<3>(7).geometry().center()[2];

        lgrZCORN[zcorn_bottom_00_idx] = elem.subEntity<3>(0).geometry().center()[2];
        lgrZCORN[zcorn_bottom_10_idx] = elem.subEntity<3>(1).geometry().center()[2];
        lgrZCORN[zcorn_bottom_01_idx] = elem.subEntity<3>(2).geometry().center()[2];
        lgrZCORN[zcorn_bottom_11_idx] = elem.subEntity<3>(3).geometry().center()[2];
    }

    // Rewrite values for active pillars
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int cell_pillar_idx = (j*nx) + i;

            // Get min/max k for pillar at (i,j)
            const auto& [bottom_k, top_k] = minMaxPerCellPillar[cell_pillar_idx];

            if ( bottom_k == nz ) {
                continue; // no active pillar at (i,j)
            }

            const auto bottom_lgr_cartesian_idx = (bottom_k*nx*ny) + cell_pillar_idx;
            const auto top_lgr_cartesian_idx = (top_k*nx*ny) + cell_pillar_idx;

            const auto& bottomElemIdx = lgrCartesianIdxToCellIdx.at(bottom_lgr_cartesian_idx);
            const auto& topElemIdx = lgrCartesianIdxToCellIdx.at(top_lgr_cartesian_idx);

            const auto& bottomElem = Dune::cpgrid::Entity<0>(levelGrid, bottomElemIdx, true);
            const auto& topElem = Dune::cpgrid::Entity<0>(levelGrid, topElemIdx, true);

            Opm::processPillars(i,j, nx, topElem, bottomElem, lgrCOORD);
        }
    }
    return std::make_pair(lgrCOORD, lgrZCORN);
}

void setPillarCoordinates(int i, int j, int nx,
                          int topCorner, int bottomCorner, int positionIdx,
                          const Dune::cpgrid::Entity<0>& topElem,
                          const Dune::cpgrid::Entity<0>& bottomElem,
                          std::vector<double>& lgrCOORD)
{
    // positionIdx (0, 1, 2, or 3) is used to distinguish the 4 corner pillars in a cell column.
    //
    // Corner pillar position mapping:
    //
    //   positionIdx   corresponding pillar   positionIdx / 2   positionIdx % 2
    //   ---------------------------------------------------------------------
    //       0          (i, j)                     0                  0
    //       1          (i+1, j)                   0                  1
    //       2          (i, j+1)                   1                  0
    //       3          (i+1, j+1)                 1                  1
    //
    // - positionIdx / 2 determines the position at the y-axis (0 for j, 1 for j+1).
    // - positionIdx % 2 determines the position at the x-axis (0 for i, 1 for i+1).

    const int pillar = ((j + positionIdx / 2) *6* (nx + 1)) + 6*(i + positionIdx % 2);

    // Top pillar's COORD values
    const auto& top_point = topElem.subEntity<3>(topCorner).geometry().center();
    std::ranges::copy(top_point, lgrCOORD.begin()+pillar);

    // Bottom pillar's COORD values
    const auto& bottom_point = bottomElem.subEntity<3>(bottomCorner).geometry().center();
    std::ranges::copy(bottom_point, lgrCOORD.begin() + pillar + 3);
}

void processPillars(int i, int j, int nx,
                    const Dune::cpgrid::Entity<0>& topElem,
                    const Dune::cpgrid::Entity<0>& bottomElem,
                    std::vector<double>& lgrCOORD)
{
    // Recall that a cell has 8 corners:
    //        6 --- 7
    //       /     /   TOP FACE
    //      4 --- 5
    //        2 --- 3
    //       /     /   BOTTOM FACE
    //      0 --- 1

    // To take into account inactive cells, consider for each (i,j) column of cells, 4 pillars:
    // (i,j)     pillar associated with bottom element corner 0 and top element corner 4
    // (i+1,j)   pillar associated with bottom element corner 1 and top element corner 5
    // (i,j+1)   pillar associated with bottom element corner 2 and top element corner 6
    // (i+1,j+1) pillar associated with bottom element corner 3 and top element corner 7
    setPillarCoordinates(i, j, nx, 4 /*topCorner*/, 0 /*bottomCorner*/,  0 /*to select pillar (i,j)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 5 /*topCorner*/, 1 /*bottomCorner*/,  1 /*to select pillar (i+1,j)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 6 /*topCorner*/, 2 /*bottomCorner*/,  2 /*to select pillar (i,j+1)*/, topElem, bottomElem, lgrCOORD);
    setPillarCoordinates(i, j, nx, 7 /*topCorner*/, 3 /*bottomCorner*/,  3 /*to select pillar (i+1,j+1)*/, topElem, bottomElem, lgrCOORD);
}


std::pair<std::vector<double>, std::vector<double>>
lgrCOORDandZCORN(const Dune::cpgrid::CpGridData& cellRefGrid,
                 const std::array<int, 3>& cellRefGrid_dim)
{
    // Check not all cells are inactive
    const auto numCells = cellRefGrid.size(0);
    if (numCells == 0) {
        OPM_THROW(std::logic_error, "Grid has no active cells.\n");
    }
    
    const int nx = cellRefGrid_dim[0];
    const int ny = cellRefGrid_dim[1];
    const int nz = cellRefGrid_dim[2];

    // Initialize all pillars as inactive (setting COORD values to std::numeric_limits<double>::max()).
    std::vector<double> coord(6*(nx+1)*(ny+1), std::numeric_limits<double>::max());

    // Initialize all ZCORN as inactive (setting values to std::numeric_limits<double>::max()).
    std::vector<double> zcorn(8*nx*ny*nz, std::numeric_limits<double>::max());

    // Map to determine min and max k per cell column (i, j) (min/max_k = 0, ..., nz-1).
    // Initialized as {nz, -1} to detect inactive cell columns.
    std::vector<std::array<int,2>> minMaxPerCellPillar(nx*ny, {nz, -1});


    for (int elemIdx = 0; elemIdx < cellRefGrid.size(0); ++elemIdx) {

       const auto ijk = Opm::Lgr::getIJK(elemIdx, cellRefGrid_dim);

         // Compute the bottom and top k per cell pillar (i, j).
        int cell_pillar_idx = ijk[1] * nx + ijk[0];
        auto& minMax = minMaxPerCellPillar[cell_pillar_idx];

        minMax[0] = std::min(ijk[2], minMax[0]);
        minMax[1] = std::max(ijk[2], minMax[1]);
    }

    for (int elemIdx = 0; elemIdx < cellRefGrid.size(0); ++elemIdx) {
         
        const auto elem = Dune::cpgrid::Entity<0>(cellRefGrid, elemIdx, true);
        const auto elemIJK = Opm::Lgr::getIJK(elemIdx, cellRefGrid_dim);

        // For a grid with nz layers, ZCORN values are ordered:
        //
        //      top layer nz-1
        //   bottom layer nz-1
        //      top layer nz-2
        //   bottom layer nz-2
        // ...
        //      top layer 1
        //   bottom layer 1
        //      top layer 0
        //   bottom layer 0

        int zcorn_top_00_idx = ((nz-1-elemIJK[2])*8*nx*ny) + (elemIJK[1]*4*nx) + (2*elemIJK[0]); // assoc. w. elem corner 4

        // Bottom indices
        int zcorn_top_10_idx = zcorn_top_00_idx + 1;  // assoc. w. elem corner 5
        int zcorn_top_01_idx = zcorn_top_00_idx + (2*nx);  // assoc. w. elem corner 6
        int zcorn_top_11_idx = zcorn_top_01_idx + 1; // assoc. w. elem corner 7

        // Top indices
        int zcorn_bottom_00_idx = zcorn_top_00_idx + (4*nx*ny); // assoc. w. elem corner 0
        int zcorn_bottom_10_idx = zcorn_bottom_00_idx + 1;  // assoc. w. elem corner 1
        int zcorn_bottom_01_idx = zcorn_bottom_00_idx + (2*nx); // assoc. w. elem corner 2
        int zcorn_bottom_11_idx = zcorn_bottom_01_idx + 1;  // assoc. w. elem corner

        // Note: zcorn_idx + 1 moves to the next position along the x-axis (i+1, j, k)
        //       zcorn_idx + (2*nx) moves to the next position along the y-axis (i, j+1, k)
        //       zcorn_idx + (4*nx*ny) moves to the next position along the z-axis (i,j, k+1)

        // Assign ZCORN values
        zcorn[zcorn_top_00_idx] = elem.subEntity<3>(4).geometry().center()[2];
        zcorn[zcorn_top_10_idx] = elem.subEntity<3>(5).geometry().center()[2];
        zcorn[zcorn_top_01_idx] = elem.subEntity<3>(6).geometry().center()[2];
        zcorn[zcorn_top_11_idx] = elem.subEntity<3>(7).geometry().center()[2];

        zcorn[zcorn_bottom_00_idx] = elem.subEntity<3>(0).geometry().center()[2];
        zcorn[zcorn_bottom_10_idx] = elem.subEntity<3>(1).geometry().center()[2];
        zcorn[zcorn_bottom_01_idx] = elem.subEntity<3>(2).geometry().center()[2];
        zcorn[zcorn_bottom_11_idx] = elem.subEntity<3>(3).geometry().center()[2];
    }

    // Rewrite values for active pillars
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const int cell_pillar_idx = (j*nx) + i;

            // Get min/max k for pillar at (i,j)
            const auto& [bottom_k, top_k] = minMaxPerCellPillar[cell_pillar_idx];

            if ( bottom_k == nz ) {
                continue; // no active pillar at (i,j)
            }

            const auto bottom_lgr_cartesian_idx = (bottom_k*nx*ny) + cell_pillar_idx;
            const auto top_lgr_cartesian_idx = (top_k*nx*ny) + cell_pillar_idx;
            // Active parent cell-> all active children (at least before processing MINPV(V) in CARFIN block???)
            const auto bottomElem = Dune::cpgrid::Entity<0>(cellRefGrid, bottom_lgr_cartesian_idx, true);
            const auto topElem = Dune::cpgrid::Entity<0>(cellRefGrid, top_lgr_cartesian_idx, true);

            Opm::processPillars(i,j, nx, topElem, bottomElem, coord);
        }
    }
    return std::make_pair(coord, zcorn);
}

std::set<Dune::FieldVector<double,3>,/*Opm::Lgr::FieldVectorLess*/ PillarLess>
computeBasicRefinedCorners(const Dune::cpgrid::Entity<0>& parentCell,
                           const std::array<int,3>& nxnynz, // or nxfin, nyfin, nzfin
                           const std::vector<double>& widthsX, // hxfin
                           const std::vector<double>& lengthsY, // hyfin
                           const std::vector<double>& heightsZ) // hzfin
{
    std::set<Dune::FieldVector<double,3>, /*Opm::Lgr::FieldVectorLess*/ PillarLess> coords{};

    const auto parentCellGeom = parentCell.geometry();

    int nx = nxnynz[0];
    int ny = nxnynz[1];
    int nz = nxnynz[2];

    // Initialize all ZCORN as inactive (setting values to std::numeric_limits<double>::max()).
    std::vector<double> zcorn(8*nx*ny*nz, std::numeric_limits<double>::max());
    
    assert(static_cast<int>(widthsX.size()) == nx);
    assert(static_cast<int>(lengthsY.size()) == ny);
    assert(static_cast<int>(heightsZ.size()) == nz);
    
    const auto localCoordNumerator = []( const std::vector<double>& vec,
                                         int sumLimit,
                                         double multiplier) {
        double lcn = 0;
        assert(!vec.empty());
        assert(sumLimit < static_cast<int>(vec.size()));
        lcn += multiplier*vec[sumLimit];
        for (int m = 0; m < sumLimit; ++m) {
            lcn += vec[m];
        }
        return lcn;
    };
    // E.g. localCoordNumerator( dx, 3, 0.25) =  x0 + x1 + x2 + 0.25.x3
    //
    const double sumWidths = std::accumulate(widthsX.begin(), widthsX.end(), double(0));
    // x0 + x1 + ... + xL, if dx = {x0, x1, ..., xL}
    const double sumLengths = std::accumulate(lengthsY.begin(), lengthsY.end(), double(0));
    // y0 + y1 + ... + yM, if dy = {y0, y1, ..., yM}
    const double sumHeights = std::accumulate(heightsZ.begin(), heightsZ.end(), double(0));
    // z0 + z1 + ... + zN, if dz = {z0, z1, ..., zN}

    for (int j = 0; j < ny +1; ++j) {
        double local_y = 0;
        for (int i = 0; i < nx +1; ++i) {
            double local_x = 0.;
            for (int k = 0; k < nz +1; ++k) {
                double local_z = 0.;

                // int refined_corner_idx = (j*(nx+1)*(nz+1)) + (i*(nz+1)) + k;
          
                if ( i == nx) { // last corner in the x-direction
                    local_x = sumWidths;
                } else {
                    local_x = localCoordNumerator(widthsX, i/nx, double((i % nx)) / nx);
                }
                if ( j == ny) { // last corner in the y-direction
                    local_y = sumLengths;
                } else {
                    local_y = localCoordNumerator(lengthsY, j/ny, double((j % ny)) / ny);
                }
                if ( k == nz) { // last corner in the z-direction
                    local_z = sumHeights;
                } else {
                    local_z = localCoordNumerator(heightsZ, k/nz, double((k % nz)) /nz);
                }

                const Dune::FieldVector<double,3> local_refined_corner = { local_x/sumWidths, local_y/sumLengths, local_z/sumHeights };
                assert(local_x/sumWidths <= 1.);
                assert(local_y/sumLengths <= 1.);
                assert(local_z/sumHeights <= 1.);
                
                coords.insert(parentCellGeom.global(local_refined_corner));
            } // end k-for-loop
        } // end i-for-loop
    } // end j-for-loop


    // Populate zcorn

    


    
    return coords;
}

void addAllParentCellFaceVertices(const Dune::cpgrid::CpGridData& grid,
                                  const Dune::cpgrid::Entity<0>& parentCell,
                                  std::set<Dune::FieldVector<double,3>, PillarLess>& input_vertices)
{
    for (const auto& face : grid.cellToFace(parentCell.index())) {
        for (const auto& point : grid.faceToPoint(face.index())) {
            input_vertices.insert(Dune::cpgrid::Entity<3>( grid, point, true).geometry().center());
        }
    }
}

std::pair<Dune::FieldVector<double,3>, double> computeCenterAndVolume(const std::array<Dune::FieldVector<double,3>,8>& corners)
{
    Dune::FieldVector<double,3> center = {0., 0., 0.};
    for (int i = 0; i < 8; ++i) {
        center += corners[i];
    }
    center /= 8.;

    double volume = 0.;
    
    static const std::vector<std::array<int,4>> cellFacesVtxIndices = {
        {0,2,6,4}, // I minus
        {1,3,7,5}, // I plus
        {2,3,7,6}, // J minus 
        {0,1,5,4}, // J plus
        {0,1,3,2}, // K minus
        {4,5,7,6}  // K plus
    };
    
    std::vector<std::vector<std::array<int,2>>> tetra_edge_indices;
    tetra_edge_indices.reserve(6);

    for (const auto& faceVtxIndices : cellFacesVtxIndices) {
        tetra_edge_indices.push_back(Opm::Lgr::createEdges(faceVtxIndices));
    }
    
    // Sum of the 24 volumes to get the volume of the hexahedron
 
    // Calculate the volume of each hexahedron, by adding
    // the 4 tetrahedra at each face (4x6 = 24 tetrahedra).
    for (int face = 0; face < 6; ++face) {
        
        const auto faceVtxIndices = cellFacesVtxIndices[face];
        const auto faceCenter = Opm::Lgr::computeFaceCenter({corners[faceVtxIndices[0]],
                corners[faceVtxIndices[1]],
                corners[faceVtxIndices[2]],
                corners[faceVtxIndices[3]]});
        
        for (int edge = 0; edge < 4; ++edge) {
            // Construction of each tetrahedron based on "face" with one
            // of its edges equal to "edge".
            const Dune::FieldVector<double,3> tetra_corners[4] = {
                corners[tetra_edge_indices[face][edge][0]],  
                corners[tetra_edge_indices[face][edge][1]],  
                faceCenter,
                center };  
            volume += std::fabs(simplex_volume(tetra_corners));
        } // end edge-for-loop
    } // end face-for-loop

    return {center, volume};
}

} // namespace Opm
