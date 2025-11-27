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

#ifndef OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED
#define OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <cstddef>      // for std::size_t
#include <utility>      // for std::move
#include <type_traits>  // for std::is_same_v
#include <vector>

namespace Opm
{
namespace Lgr
{

/// @brief Builds a mapping from level element indices to the Cartesian ordering required by output files.
///
/// Output file solution data containers for LGRs expect elements to be ordered in strict Cartesian order,
/// with the i-direction varying fastest, followed by j, then k.
///
/// @param [in] grid          The CpGrid instance representing the simulation grid.
/// @param [in] levelCartMapp The LevelCartesianIndexMapper providing Cartesian index information across all levels.
/// @param [in] level         The grid level (integer) for which to build the index mapping.
/// @return A vector where each position corresponds to the output ordering, and the value is the element index
///         in that position. This ensures compatibility with the Cartesian ordering expected by output files.
std::vector<int> mapLevelIndicesToCartesianOutputOrder(const Dune::CpGrid& grid,
                                                       const Opm::LevelCartesianIndexMapper<Dune::CpGrid>& levelCartMapp,
                                                       int level);

/// @brief Reorder data from a simulation container into the order assumed by output for refined level grids.
///
/// @param [in] simulatorContainer  Container with simulation data ordered by compressed indices.
/// @param [in] toOutput            A vector where each position corresponds to the output ordering, and the
///                                 value is the element index in that position.
/// @return container with reordered data as expected by output, i.e. in strict Cartesian order, with the
///         i-direction varying fastest, followed by j, then k.
template <typename Container>
Container reorderForOutput(const Container& simulatorContainer,
                           const std::vector<int>& toOutput);

/// @brief Creates a mapping from level-element indices to corresponding leaf-element indices.
///
/// For each level of the grid, this function maps every element index to the
/// index of its corresponding element index on the leaf grid. If a level-element
/// does not exist on the leaf (i.e., it is a parent cell with no leaf representation),
/// the value std::numeric_limits<int>::max() is stored to indicate an invalid
/// non-leaf entry.
///
/// @param [in] grid
/// @return A vector each entry represents the map between level indices
///         and leaf indices, per level grid.
///         map[ level ][ level element idx ] = leaf index, if leaf; rubbish otherwise.
std::vector<std::vector<int>> levelIdxToLeafIdxMaps(const Dune::CpGrid& grid);

/// @brief Collect level and level indices of leaf descendants.
///
/// This function gathers the {level, level-index} pairs of all *leaf* descendants
/// of a given element. It may only be called on elements that are *not* leaf
/// themselves. If the input element already appears on the leaf grid, this
/// method throws.
///
/// Example: Supposed an element is refined into nx.ny.nz children in level l1, with indices
/// in level l1 { i1, i2, ..., iN}. One of its chidldren in l1, let's say i2, got refined
/// into n1x.n1y.n1z in level l2 with level l2 indices {j1, j2, ..., jM}. No further refinement
/// is performed. Then, the resulting list of level and level indices of leaf descendants is:
/// { {l1, i1}, {l2, j1}, {l2, j2}, ...,{l2,jM}, {l1, i3}, ..., {l1, iN} }
///
/// @param [in] element   Element whose leaf descendants should be collected.
/// @parem [in] maxLevel  Maximum refinement level to consider.
/// @return A vector of {level, level index} pairs for all leaf descendants.
std::vector<std::array<int,2>> getLevelAndLevelIdxOfLeafDescendants(const Dune::cpgrid::Entity<0>& element,
                                                                    int maxLevel);

/// @breif Compute the average of the children data values
///
/// @param [in] levelVectors A collection of per-level data vectors.
///                          levelVectors[l][i] contains the data value for the element at level l
///                          with level index i. To extract the children data.
/// @param [in] element      A non-leaf element (parent cell). Its children data
///                          will be used to compute the average.
/// @param [in] grid         The grid from which the element is taken.
/// @return Average of children data values. If ScalarTypeis double, it's the pore volume average.
template <typename ScalarType>
ScalarType averageChildrenData(const std::vector<std::vector<ScalarType>>& levelVectors,
                               const Dune::cpgrid::Entity<0>& element,
                               const Dune::CpGrid& grid,
                               const std::vector<double>& porv_levelZero);

/// @brief Populate level data vectors based on leaf vector, for a specific named data field.
///
/// @param [in]       grid
/// @param [in]       maxLevel
/// @param [in]       leafVector Containing one named data field (e.g., pressure, saturation).
/// @param [in]       toOutput_refinedLevels For level grids 1,2,..,maxLevel, a map to store
///                                          data in the order expected by outout files
///                                          (increasing level Cartesian indices).
/// @param [in]       porv_levelZero Vector containing pore volume cell data for level-zero grid cells.
/// @param [out]      levelVectors A collection of per-level data vectors.
///                                levelVectors[l][i] contains the data value for the element at level l
///                                with level index i.
template <typename ScalarType>
void populateDataVectorLevelGrids(const Dune::CpGrid& grid,
                                  int maxLevel,
                                  const std::vector<ScalarType>& leafVector,
                                  const std::vector<std::vector<int>>& toOutput_refinedLevels,
                                  const std::vector<double>& porv_levelZero,
                                  std::vector<std::vector<ScalarType>>& levelVectors);

/// @brief Extracts and organizes solution data for all grid refinement levels.
///
/// It derives these level-specific solutions from a given leaf-solution. For cells
/// that no longer exist in the leaf grid (i.e., parent cells that were refined away),
/// the average of its children values is assigned.
/// This behavior is temporary and will be replaced with proper restriction in a future
/// implementation.
///
/// Each resulting ScalarBuffer (std::vector<Scalar>) follows the ordering expected by
/// output files, i.e., increasing level Cartesian indices.
///
/// @param [in]       grid
/// @param [in]       toOutput_refinedLevels For level grids 1,2,..,maxLevel, a map to store
///                                          data in the order expected by outout files
///                                          (increasing level Cartesian indices)
/// @param [in]       leafSolution The complete solution defined on the leaf grid, containing
///                                one or more named data fields (e.g., pressure, saturation).
/// @param [in] porv_levelZero Vector containing pore volume cell data for level-zero grid active cells.
/// @param [out]   A vector of Opm::data::Solution objects, one for each refinement level
///                (from level 0 to grid.maxLevel()), where each entry contains data reordered
///                according to increasing level Cartesian indices for output.
void extractSolutionLevelGrids(const Dune::CpGrid& grid,
                               const std::vector<std::vector<int>>& toOutput_refinedLevels,
                               const Opm::data::Solution& leafSolution,
                               const std::vector<double>& porv_levelZero,
                               std::vector<Opm::data::Solution>&);

/// @brief Constructs restart-value containers for all grid refinement levels.
///
/// The level-specific solution data are first derived from the leaf solution
/// using extractSolutionLevelGrids(...). Other data components (such as wells,
/// group/network values, and aquifers) are passed unchanged to each level.
///
/// @param [template] Grid The function has no effect for grids other than CpGrid.
/// @param [in]       grid
/// @param [in]       leafRestartValue
/// @param [in]       porv_levelZero Vector containing cell pore volume data for level-zero grid cells.
/// @param [out]      A vector of RestartValue objects, one for each refinement level
///                   (from level 0 to grid.maxLevel()).
template <typename Grid>
void extractRestartValueLevelGrids(const Grid& grid,
                                   const Opm::RestartValue& leafRestartValue,
                                   const std::vector<double>& porv_levelZero,
                                   std::vector<Opm::RestartValue>& restartValue_levels);

} // namespace Lgr
} // namespace Opm

template <typename Container>
Container Opm::Lgr::reorderForOutput(const Container& simulatorContainer,
                                     const std::vector<int>& toOutput)
{
    // Use toOutput to reorder simulatorContainer
    Container outputContainer;
    outputContainer.resize(toOutput.size());
    for (std::size_t i = 0; i < toOutput.size(); ++i) {
        outputContainer[i] = simulatorContainer[toOutput[i]];
    }
    return outputContainer;
}

template <typename ScalarType>
ScalarType Opm::Lgr::averageChildrenData(const std::vector<std::vector<ScalarType>>& levelVectors,
                                         const Dune::cpgrid::Entity<0>& element,
                                         const Dune::CpGrid& grid,
                                         const std::vector<double>& porv_levelZero)
{
    ScalarType partialSum{};

    auto termToAdd = [&levelVectors,
                      &grid,
                      &porv_levelZero](int level, int levelIdx)
                      {
                          ScalarType term{};
                          if constexpr( std::is_same_v<ScalarType, int>) {
                              term = levelVectors[ level ][ levelIdx ];
                          }
                          else {
                              const auto& child = Dune::cpgrid::Entity<0>(*grid.currentData()[level], levelIdx, true);
                              // (!) for nested refinement: instead of child.father(), use child.getOrigin()
                              const auto& origin = child.getOrigin();
                              // child_porv =  parent-cell-porv * child-cell-volume/parent-cell-volume
                              const auto child_porv = porv_levelZero[ origin.index() ]*child.geometry().volume()/origin.geometry().volume();
                              term = levelVectors[ level ][ levelIdx ]*child_porv;
                          }
                          return term;
                      };

    int count = 0;
    for (const auto& level_levelIdx : getLevelAndLevelIdxOfLeafDescendants(element, grid.maxLevel())) {
        partialSum += termToAdd(level_levelIdx[0], level_levelIdx[1]);
        ++count;
    }
    // If ScalarType == int, when dividing by int 'count', it'd be int division.
    if constexpr( std::is_same_v<ScalarType, int> ) {
        return partialSum/count; // (sum_{i in children} child_i_data) / total_children
    }
    else { // element pore volume = parent-cell-porv * element-cell-volume/parent-cell-volume
        double elemPorv = porv_levelZero[ element.getOrigin().index() ]*element.geometry().volume()/element.getOrigin().geometry().volume();
        return partialSum/elemPorv;  // (sum_{i in children} child_i_data * child_i_pore_volume)/ element_pore_volume
    }
}

template <typename ScalarType>
void Opm::Lgr::populateDataVectorLevelGrids(const Dune::CpGrid& grid,
                                            int maxLevel,
                                            const std::vector<ScalarType>& leafVector,
                                            const std::vector<std::vector<int>>& toOutput_refinedLevels,
                                            const std::vector<double>& porv_levelZero,
                                            std::vector<std::vector<ScalarType>>& levelVectors)
{
    for (int level = 0; level <= maxLevel; ++level) {
        levelVectors[level].resize(grid.levelGridView(level).size(0));
    }
    // For level cells that appear in the leaf, extract the data value from leafVector
    // and assign it to the equivalent level cell.
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        levelVectors[element.level()][element.getLevelElem().index()] = leafVector[element.index()];
    }
    // Note that all cells from maxLevel have assigned values at this point.
    // Now, assign values for parent cells (for now, average of children values).
    if (maxLevel)  {
        for (int level = maxLevel-1; level >= 0; --level) {
            for (const auto& element : Dune::elements(grid.levelGridView(level))) {
                if (!element.isLeaf()) {
                    levelVectors[level][element.index()] = Opm::Lgr::averageChildrenData(levelVectors,
                                                                                         element,
                                                                                         grid,
                                                                                         porv_levelZero);
                }
            }
        }
    }
    // Use toOutput_levels to reorder in ascending level cartesian indices
    for (int level = 1; level<=maxLevel; ++level) { // exclude level zero (does not need reordering)
        levelVectors[level] = Opm::Lgr::reorderForOutput(levelVectors[level], toOutput_refinedLevels[level-1]);
    }
}

template <typename Grid>
void Opm::Lgr::extractRestartValueLevelGrids(const Grid& grid,
                                             const Opm::RestartValue& leafRestartValue,
                                             const std::vector<double>& porv_levelZero,
                                             std::vector<Opm::RestartValue>& restartValue_levels)
{
    if constexpr (std::is_same_v<Grid, Dune::CpGrid>) {

        int maxLevel = grid.maxLevel();
        restartValue_levels.resize(maxLevel+1); // level 0, 1, ..., max level
        
        // To store leafRestartValue.extra data in the order expected
        // by outout files (increasing level Cartesian indices)
        std::vector<std::vector<int>> toOutput_refinedLevels{};
        toOutput_refinedLevels.resize(maxLevel); // exclude level zero (does not need reordering)

        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (int level = 1; level <= maxLevel; ++level) { // exclude level zero (does not need reordering)
            toOutput_refinedLevels[level-1] = mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
        }

        std::vector<Opm::data::Solution> dataSolutionLevels{};
        extractSolutionLevelGrids(grid,
                                  toOutput_refinedLevels,
                                  leafRestartValue.solution,
                                  porv_levelZero,
                                  dataSolutionLevels);



        for (int level = 0; level <= maxLevel; ++level) {
            restartValue_levels[level] = Opm::RestartValue(std::move(dataSolutionLevels[level]),
                                                           leafRestartValue.wells,
                                                           leafRestartValue.grp_nwrk,
                                                           leafRestartValue.aquifer,
                                                           level);
        }

        for (const auto& [rst_key, leafVector] : leafRestartValue.extra) {

            std::vector<std::vector<double>> levelVectors{};
            levelVectors.resize(maxLevel+1);

            Opm::Lgr::populateDataVectorLevelGrids<double>(grid,
                                                           maxLevel,
                                                           leafVector,
                                                           toOutput_refinedLevels,
                                                           porv_levelZero,
                                                           levelVectors);

            for (int level = 0; level <= maxLevel; ++level) {
                restartValue_levels[level].addExtra(rst_key.key, rst_key.dim, std::move(levelVectors[level]));
            }
        }
    }
}

#endif // OPM_GRID_CPGRID_LGROUTPUTHELPERS_HEADER_INCLUDED
