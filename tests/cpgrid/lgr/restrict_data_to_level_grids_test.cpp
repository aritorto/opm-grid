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
#include "config.h"

#define BOOST_TEST_MODULE RestrictDataToLevelGridsTests
#include <boost/test/unit_test.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>

#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>
#include <opm/grid/cpgrid/LgrHelpers.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>
#include <opm/output/data/Solution.hpp>

#include <array>
#include <numeric>
#include <string>  // for std::iota
#include <utility> // for std::move
#include <vector>

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

void restrictFakeLeafDataToLevelGrids(const Dune::CpGrid& grid)
{
    // Create a fake leaf solution
    Opm::data::Solution leafSolution{};
    BOOST_CHECK(!leafSolution.has("NOTHING")); // dummy check

    std::vector<double> leafVec{};
    leafVec.resize(grid.leafGridView().size(0));
    std::iota(leafVec.begin(), leafVec.end(), 0);  // fills with 0, 1, 2, ..., total amount of leaf cells -1

    leafSolution.insert("FAKEPROP",
                        Opm::UnitSystem::measure::liquid_surface_volume, // just one possible meassure
                        std::move(leafVec),
                        Opm::data::TargetType::RESTART_OPM_EXTENDED);    // just one possible TargetType

    BOOST_CHECK(leafSolution.has("FAKEPROP"));

    const auto& leafFakePropData = leafSolution.data<double>("FAKEPROP");
    BOOST_CHECK_EQUAL(leafFakePropData.size(), grid.leafGridView().size(0));

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    std::vector<Opm::data::Solution> levelSolutions{};
    levelSolutions.resize(grid.maxLevel()+1);

    for (const auto& leaf_data : leafSolution)
    {
        std::vector<std::vector<double>> aux_levels_data{};
        aux_levels_data.resize(grid.maxLevel()+1);

        for (int level = 0; level <=grid.maxLevel(); ++level) {
            aux_levels_data[level].resize(grid.levelGridView(level).size(0));
        }

        // For level cells that appear in the leaf, extract the data value from leaf_data
        for (const auto& element : Dune::elements(grid.leafGridView())) {
            aux_levels_data[element.level()][element.getLevelElem().index()] = leaf_data.second.data<double>()[element.index()];
        }

        // For the level cells that vanished (don't make it to the leaf), compute avarage of children data values
        for (int level = grid.maxLevel()-1; level >= 0; --level) {
            for (const auto& element : Dune::elements(grid.levelGridView(level))) {
                if (element.isLeaf()) // entry already populated
                    continue;

                auto it = element.hbegin(level);
                const auto& endIt = element.hend(level);

                int child_count = 0;
    
                for (; it != endIt; ++it) {
                    aux_levels_data[level][element.index()] += aux_levels_data[it->level()][it->index()];
                    ++child_count;
                }
                aux_levels_data[level][element.index()] /= child_count;
            }

        }

        // By now, all level cells have data assigned. Reorder the containers in the order
        // expected by outout files (increasing level Cartesian indices)
        const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
        for (int level = 0; level <= grid.maxLevel(); ++level) {

            const auto toOutput = Opm::Lgr::mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
            const auto outputContainer = Opm::Lgr::reorderForOutput( aux_levels_data[level],
                                                                     toOutput);

            levelSolutions[level].insert(leaf_data.first,
                                         Opm::UnitSystem::measure::liquid_surface_volume,
                                         std::move(outputContainer),
                                         Opm::data::TargetType::RESTART_OPM_EXTENDED);

            BOOST_CHECK(levelSolutions[level].has("FAKEPROP"));

            const auto& levelFakePropData = levelSolutions[level].data<double>("FAKEPROP");
            BOOST_CHECK_EQUAL(levelFakePropData.size(), grid.levelGridView(level).size(0));
        }
    }
}

BOOST_AUTO_TEST_CASE(restrictDataForNonNestedLgrsSharingEdges)
{
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {4,3,3}, /* cell_sizes = */ {1.,1.,1.});
    //                          LGR1 parent cells          LGR2 parent cells
    // k = 2   32 33 34 35 |
    //         28 29 30 31 |          29 30
    //         24 25 26 27 |          25 26
    // --------------------|    ------------------       ------------------
    // k = 1   20 21 22 23 |
    //         16 17 18 19 |          17 18
    //         12 13 14 15 |          13 14
    // --------------------|    ------------------       ------------------
    // k = 0    8  9 10 11 |                                   9 10
    //          4  5  6  7 |
    //          0  1  2  3 |
    //---------------------|
    // Refine a few cells, each into 2x2x2 children (2 subdivisions per x-,y-,z- direction).
    grid.addLgrsUpdateLeafView(/* cells_per_dim_vec = */ {{2,2,2}, {2,3,2}},
                               /* startIJK_vec = */ {{1,0,1}, {1,2,0}},
                               /* endIJK_vec = */ {{3,2,3}, {3,3,1}},
                               /* lgr_name_vec = */ {"LGR1", "LGR2"});

    restrictFakeLeafDataToLevelGrids(grid);
}


BOOST_AUTO_TEST_CASE(restrictDataForNestedRefinementOnly) {

    const std::vector<std::array<int,3>>  cells_per_dim_vec = {{3,2,2}, {2,2,2}, {2,2,1}, {3,4,2}};
    const std::vector<std::array<int,3>>       startIJK_vec = {{1,1,0}, {1,1,0}, {1,1,0}, {1,1,0}};
    const std::vector<std::array<int,3>>         endIJK_vec = {{2,2,1}, {2,2,1}, {2,2,1}, {2,2,1}};
    const std::vector<std::string>             lgr_name_vec = { "LGR1",  "LGR2",  "LGR3",  "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL", "LGR1",  "LGR2",  "LGR3"};

    // GLOBAL            The grids are stored: GLOBAL,
    //  |                                      LGR1 (child grid from GLOBAL),
    // LGR1                                    LGR2 (child grid from LGR1),
    //  |                                      LGR3 (child grid from LGR2),
    // LGR2                                    LGR4 (child grid from LGR3),
    //  |                                      leaf grid view (without name).
    // LGR3
    //  |
    // LGR4
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                               startIJK_vec,
                               endIJK_vec,
                               lgr_name_vec,
                               lgr_parent_grid_name_vec);

    restrictFakeLeafDataToLevelGrids(grid);
}

BOOST_AUTO_TEST_CASE(restrictDataForMixNameOrderAndNestedRefinement){

    const std::vector<std::array<int,3>>  cells_per_dim_vec = { {3,2,2}, {2,2,2}, {2,2,1}, {3,4,2}};
    const std::vector<std::array<int,3>>       startIJK_vec = { {1,1,0}, {0,0,0}, {0,0,0}, {1,1,0}};
    const std::vector<std::array<int,3>>         endIJK_vec = { {2,2,1}, {1,1,1}, {1,1,1}, {2,2,1}};
    const std::vector<std::string>             lgr_name_vec = {  "LGR1", "LGR2",   "LGR3", "LGR4"};
    const std::vector<std::string> lgr_parent_grid_name_vec = {"GLOBAL", "LGR1", "GLOBAL", "LGR3"};

    //   GLOBAL            The grids are stored: GLOBAL,
    //   |    |                                  LGR1, LGR3 (child grid from GLOBAL),
    // LGR1  LGR3                                LGR2 (child grid from LGR1),
    //   |    |                                  LGR4 (child grid from LGR3),
    // LGR2  LGR4                                leaf grid view (without name).
    Dune::CpGrid grid;
    grid.createCartesian(/* grid_dim = */ {3,3,1}, /* cell_sizes = */ {1.0, 1.0, 1.0});

    grid.addLgrsUpdateLeafView(cells_per_dim_vec,
                               startIJK_vec,
                               endIJK_vec,
                               lgr_name_vec,
                               lgr_parent_grid_name_vec);

    restrictFakeLeafDataToLevelGrids(grid);
}
