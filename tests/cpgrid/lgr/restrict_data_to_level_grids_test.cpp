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

#include <dune/grid/common/mcmgmapper.hh>

#include <array>
#include <string>
#include <tuple>
#include <utility>  // for std::move
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

template <typename Scalar, int numPhases> 
std::vector<Opm::data::Solution> restrictScalarBufferToLevelGrids(const Dune::CpGrid& grid,
                                                                  const Opm::data::Solution& leafSolution,
                                                                  int gasPhaseIdx,
                                                                  int oilPhaseIdx,
                                                                  int waterPhaseIdx)
{
    using DataEntry = std::tuple<std::string, Opm::UnitSystem::measure, std::vector<Scalar>&>;

    using ScalarBuffer = std::vector<Scalar>;

    // if index not specified, we treat it as valid (>= 0)
    auto addEntryIf = [&](std::vector<std::vector<DataEntry>>& container_levels,
                          const std::string& name,
                          Opm::UnitSystem::measure measure,
                          ScalarBuffer& flowArray_level, int level, int index = 1) 
    {
        assert(container_levels.size() == grid.maxLevel()+1);

        if (index >= 0 &&  leafSolution.has(name)) {  // Only add if index is valid

            const auto& leaf_data = leafSolution.template data<Scalar>(name);

            flowArray_level.resize(grid.levelGridView(level).size(0));

            // For level cells that appear in the leaf, extract the data value from leaf_data
            for (const auto& element : Dune::elements(grid.leafGridView())) {
                if (element.level() == level) {
                    flowArray_level[element.getLevelElem().index()] = leaf_data.second.template data<Scalar>()[element.index()];
                }
            }
            // Reorder the containers in the order expected by outout files (increasing level Cartesian indices)
            const Opm::LevelCartesianIndexMapper<Dune::CpGrid> levelCartMapp(grid);
            const auto toOutput = Opm::Lgr::mapLevelIndicesToCartesianOutputOrder(grid, levelCartMapp, level);
            const auto outputContainer = Opm::Lgr::reorderForOutput( flowArray_level,
                                                                     toOutput);

            container_levels[level].emplace_back(name, measure, outputContainer);
        }
    };


    int maxLevel = grid.maxLevel();

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    std::vector<Opm::data::Solution> levelSolutions{};
    levelSolutions.resize(maxLevel+1);

    std::vector<std::vector<DataEntry>> baseSolutionVector_levels{};
    baseSolutionVector_levels.resize(maxLevel+1);

    for (int level = 0; level <= maxLevel; ++level) {

        ScalarBuffer gasFormationVolumeFactor_level{};
        ScalarBuffer hydrocarbonPoreVolume_level{};
        ScalarBuffer pressureTimesPoreVolume_level{};
        ScalarBuffer pressureTimesHydrocarbonVolume_level{};
        ScalarBuffer dynamicPoreVolume_level{};
        ScalarBuffer rPorV_level{};
        ScalarBuffer fluidPressure_level{};
        ScalarBuffer temperature_level{};
        ScalarBuffer rs_level{};
        ScalarBuffer rsw_level{};
        ScalarBuffer rv_level{};
        ScalarBuffer rvw_level{};
        ScalarBuffer overburdenPressure_level{};
        ScalarBuffer oilSaturationPressure_level{};
        ScalarBuffer drsdtcon_level{};
        ScalarBuffer sSol_level{};
        ScalarBuffer rswSol_level{};
        ScalarBuffer cPolymer_level{};
        ScalarBuffer cFoam_level{};
        ScalarBuffer cSalt_level{};
        ScalarBuffer pSalt_level{};
        ScalarBuffer permFact_level{};
        ScalarBuffer soMax_level{};
        ScalarBuffer swMax_level{};
        ScalarBuffer sgmax_level{};
        ScalarBuffer shmax_level{};
        ScalarBuffer somin_level{};
        ScalarBuffer swmin_level{};
        ScalarBuffer ppcw_level{};
        ScalarBuffer gasDissolutionFactor_level{};
        ScalarBuffer oilVaporizationFactor_level{};
        ScalarBuffer gasDissolutionFactorInWater_level{};
        ScalarBuffer waterVaporizationFactor_level{};
        ScalarBuffer bubblePointPressure_level{};
        ScalarBuffer dewPointPressure_level{};
        ScalarBuffer rockCompPorvMultiplier_level{};
        ScalarBuffer minimumOilPressure_level{};
        ScalarBuffer saturatedOilFormationVolumeFactor_level{};
        ScalarBuffer rockCompTransMultiplier_level{};
        ScalarBuffer pcgw_level{};
        ScalarBuffer pcow_level{};
        ScalarBuffer pcog_level{};

        std::array<ScalarBuffer, numPhases> saturation_level{};
        std::array<ScalarBuffer, numPhases> invB_level{};
        std::array<ScalarBuffer, numPhases> density_level{};
        std::array<ScalarBuffer, numPhases> viscosity_level{};
        std::array<ScalarBuffer, numPhases> relativePermeability_level{};

        addEntryIf(baseSolutionVector_levels, "1OVERBG",  Opm::UnitSystem::measure::gas_inverse_formation_volume_factor,   invB_level[gasPhaseIdx], level, gasPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "1OVERBO",  Opm::UnitSystem::measure::oil_inverse_formation_volume_factor,   invB_level[oilPhaseIdx], level, oilPhaseIdx);

        // avoid output with generic fluid system and disabled water phase
        if constexpr (numPhases > 2) {
            addEntryIf(baseSolutionVector_levels, "1OVERBW", Opm::UnitSystem::measure::water_inverse_formation_volume_factor, invB_level[waterPhaseIdx], level, waterPhaseIdx);
        }
        addEntryIf(baseSolutionVector_levels, "FOAM",     Opm::UnitSystem::measure::identity,                              cFoam_level, level);
        addEntryIf(baseSolutionVector_levels, "GASKR",    Opm::UnitSystem::measure::identity,                              relativePermeability_level[gasPhaseIdx], level, gasPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "GAS_DEN",  Opm::UnitSystem::measure::density,                               density_level[gasPhaseIdx], level, gasPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "GAS_VISC", Opm::UnitSystem::measure::viscosity,                             viscosity_level[gasPhaseIdx], level, gasPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "OILKR",    Opm::UnitSystem::measure::identity,                              relativePermeability_level[oilPhaseIdx], level, oilPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "OIL_DEN",  Opm::UnitSystem::measure::density,                               density_level[oilPhaseIdx], level, oilPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "OIL_VISC", Opm::UnitSystem::measure::viscosity,                             viscosity_level[oilPhaseIdx], level, oilPhaseIdx);
        addEntryIf(baseSolutionVector_levels, "PBUB",     Opm::UnitSystem::measure::pressure,                              bubblePointPressure_level, level);
        addEntryIf(baseSolutionVector_levels, "PCGW",     Opm::UnitSystem::measure::pressure,                              pcgw_level, level);
        addEntryIf(baseSolutionVector_levels, "PCOG",     Opm::UnitSystem::measure::pressure,                              pcog_level, level);
        addEntryIf(baseSolutionVector_levels, "PCOW",     Opm::UnitSystem::measure::pressure,                              pcow_level, level);
        addEntryIf(baseSolutionVector_levels, "PDEW",     Opm::UnitSystem::measure::pressure,                              dewPointPressure_level, level);
        addEntryIf(baseSolutionVector_levels, "POLYMER",  Opm::UnitSystem::measure::identity,                              cPolymer_level, level);
        addEntryIf(baseSolutionVector_levels, "PPCW",     Opm::UnitSystem::measure::pressure,                              ppcw_level, level);
        addEntryIf(baseSolutionVector_levels, "PRESROCC", Opm::UnitSystem::measure::pressure,                              minimumOilPressure_level, level);
        addEntryIf(baseSolutionVector_levels, "PRESSURE", Opm::UnitSystem::measure::pressure,                              fluidPressure_level, level);
        addEntryIf(baseSolutionVector_levels, "RPORV",    Opm::UnitSystem::measure::volume,                                rPorV_level, level);
        addEntryIf(baseSolutionVector_levels, "RS",       Opm::UnitSystem::measure::gas_oil_ratio,                         rs_level, level);
        addEntryIf(baseSolutionVector_levels, "RSSAT",    Opm::UnitSystem::measure::gas_oil_ratio,                         gasDissolutionFactor_level, level);
        addEntryIf(baseSolutionVector_levels, "RV",       Opm::UnitSystem::measure::oil_gas_ratio,                         rv_level, level);
        addEntryIf(baseSolutionVector_levels, "RVSAT",    Opm::UnitSystem::measure::oil_gas_ratio,                         oilVaporizationFactor_level, level);
        addEntryIf(baseSolutionVector_levels, "SALT",     Opm::UnitSystem::measure::concentration,                         cSalt_level, level);
        addEntryIf(baseSolutionVector_levels, "SGMAX",    Opm::UnitSystem::measure::identity,                              sgmax_level, level);
        addEntryIf(baseSolutionVector_levels, "SHMAX",    Opm::UnitSystem::measure::identity,                              shmax_level, level);
        addEntryIf(baseSolutionVector_levels, "SOMAX",    Opm::UnitSystem::measure::identity,                              soMax_level, level);
        addEntryIf(baseSolutionVector_levels, "SOMIN",    Opm::UnitSystem::measure::identity,                              somin_level, level);
        addEntryIf(baseSolutionVector_levels, "SSOLVENT", Opm::UnitSystem::measure::identity,                              sSol_level, level);
        addEntryIf(baseSolutionVector_levels, "SWHY1",    Opm::UnitSystem::measure::identity,                              swmin_level, level);
        addEntryIf(baseSolutionVector_levels, "SWMAX",    Opm::UnitSystem::measure::identity,                              swMax_level, level);

        // avoid output with generic fluid system and disabled water phase
        if constexpr (numPhases > 2) {
            addEntryIf(baseSolutionVector_levels, "WATKR",    Opm::UnitSystem::measure::identity,                          relativePermeability_level[waterPhaseIdx], level, waterPhaseIdx);
            addEntryIf(baseSolutionVector_levels, "WAT_DEN",  Opm::UnitSystem::measure::density,                           density_level[waterPhaseIdx], level, waterPhaseIdx);
            addEntryIf(baseSolutionVector_levels, "WAT_VISC", Opm::UnitSystem::measure::viscosity,                         viscosity_level[waterPhaseIdx], level, waterPhaseIdx);
        }

        /* auto extendedSolutionArrays = std::array {
            DataEntry{"DRSDTCON", Opm::UnitSystem::measure::gas_oil_ratio_rate, drsdtcon_level},
            DataEntry{"PERMFACT", Opm::UnitSystem::measure::identity,           permFact_level},
            DataEntry{"PORV_RC",  Opm::UnitSystem::measure::identity,           rockCompPorvMultiplier_level},
            DataEntry{"PRES_OVB", Opm::UnitSystem::measure::pressure,           overburdenPressure_level},
            DataEntry{"RSW",      Opm::UnitSystem::measure::gas_oil_ratio,      rsw_level},
            DataEntry{"RSWSAT",   Opm::UnitSystem::measure::gas_oil_ratio,      gasDissolutionFactorInWater_level},
            DataEntry{"RSWSOL",   Opm::UnitSystem::measure::gas_oil_ratio,      rswSol_level},
            DataEntry{"RVW",      Opm::UnitSystem::measure::oil_gas_ratio,      rvw_level},
            DataEntry{"RVWSAT",   Opm::UnitSystem::measure::oil_gas_ratio,      waterVaporizationFactor_level},
            DataEntry{"SALTP",    Opm::UnitSystem::measure::identity,           pSalt_level},
            DataEntry{"TMULT_RC", Opm::UnitSystem::measure::identity,           rockCompTransMultiplier_level},
            };*/

        auto doInsert = [&levelSolutions, &level](DataEntry&       entry,
                                                  const Opm::data::TargetType target)
        {
            if (std::get<2>(entry).empty()) {
                return;
            }

            levelSolutions[level].insert(std::get<std::string>(entry),
                                         std::get<Opm::UnitSystem::measure>(entry),
                                         std::move(std::get<2>(entry)),
                                         target);
        };

        for (auto& array : baseSolutionVector_levels[level]) {
            doInsert(array, Opm::data::TargetType::RESTART_SOLUTION, level);
        }
    }
    return levelSolutions;
}


// BioeffectsContainer<Scalar> bioeffectsC_;
// ExtboContainer<Scalar> extboC_;


void restrictFakeLeafDataToLevelGrids(const Dune::CpGrid& grid,
                                      const std::vector<std::vector<double>>& expected_data_levels)
{
    // Create a fake leaf solution
    Opm::data::Solution leafSolution{};
    BOOST_CHECK(!leafSolution.has("NOTHING")); // dummy check

    const auto& leafView = grid.leafGridView();

    using LeafMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LeafGridView>;
    const LeafMapper leafMapper(leafView, Dune::mcmgElementLayout());

    std::vector<double> leafVec{};
    leafVec.resize(leafView.size(0));
    for (const auto& element : Dune::elements(leafView)) {
        leafVec[leafMapper.index(element)] = element.hasFather()? 1 : 0;
    }
    // leafVec is a vector of 0's and 1's.
    // Set 1 for those cells that have a father, 0 otherwise.

    leafSolution.insert("FAKEPROP",
                        Opm::UnitSystem::measure::liquid_surface_volume, // just one possible meassure
                        std::move(leafVec),
                        Opm::data::TargetType::RESTART_OPM_EXTENDED);    // just one possible TargetType

    BOOST_CHECK(leafSolution.has("FAKEPROP"));

    const auto& leafFakePropData = leafSolution.data<double>("FAKEPROP");
    BOOST_CHECK_EQUAL(leafFakePropData.size(), leafView.size(0));

    int maxLevel = grid.maxLevel();

    // To restrict/create the level cell data, based on the leaf cells and the hierarchy
    std::vector<Opm::data::Solution> levelSolutions{};
    levelSolutions.resize(maxLevel+1);

    using LevelMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Dune::CpGrid::LevelGridView>;
    for (const auto& leaf_data : leafSolution)
    {
        std::vector<std::vector<double>> aux_levels_data{};
        aux_levels_data.resize(maxLevel+1);

        std::vector<LevelMapper> levelMappers{};
        levelMappers.reserve(maxLevel+1);

        for (int level = 0; level <=grid.maxLevel(); ++level) {
            aux_levels_data[level].resize(grid.levelGridView(level).size(0));
            levelMappers.push_back(LevelMapper{grid.levelGridView(level), Dune::mcmgElementLayout()});
        }

        // For level cells that appear in the leaf, extract the data value from leaf_data
        for (const auto& element : Dune::elements(leafView)) {

            int elemLevel = element.level();
            const auto& equivalentLevelElement = element.getLevelElem();
            /** Entity::getLevelElem() Not DUNE interface.
                Alternative: use map CpGridData leaf_to_level_cells_
                leaf_to_level_cells_[ leaf_cell_idx ] = {level (cell parent grid index), cell index in that level} */

            aux_levels_data[elemLevel][levelMappers[elemLevel].index(equivalentLevelElement)] =
                leaf_data.second.data<double>()[leafMapper.index(element)];
        }

        // For the level cells that vanished (don't make it to the leaf), compute avarage of children data values
        for (int level = grid.maxLevel()-1; level >= 0; --level) {
            for (const auto& element : Dune::elements(grid.levelGridView(level))) {
                if (element.isLeaf()) // entry already populated
                    continue;

                auto it = element.hbegin(level+1);
                const auto& endIt = element.hend(level+1);

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
            BOOST_CHECK_EQUAL_COLLECTIONS(levelFakePropData.begin(), levelFakePropData.end(),
                                          expected_data_levels[level].begin(), expected_data_levels[level].end());
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

    std::vector<std::vector<double>> expected_data_levels{};
    expected_data_levels.resize(grid.maxLevel()+1);

    // Expected data for level grids is 1 for parent and children cells, 0 otherwise.

    // level 0 data vector size 36, with all entries 0 except for parent cells (1)
    expected_data_levels[0] = std::vector<double>(36,0);
    for (const auto& idx : std::vector<int>{9,10,13,14,17,18,25,26,29,30}) {
        expected_data_levels[0][idx] = 1;
    }

    // level 1 data vector of 1's, size 64 (cells_per_dim = {2,2,2}, block parent cells dim = {2,2,2})
    expected_data_levels[1] = std::vector<double>(64,1);

    // level 2 data vector of 1's, size 24 (cells_per_dim = {2,3,2}, block parent cells dim = {2,1,1})
    expected_data_levels[2] = std::vector<double>(24,1);

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
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

    std::vector<std::vector<double>> expected_data_levels{};
    expected_data_levels.resize(grid.maxLevel()+1);

    // Expected data for level grids is 1 for parent and children cells, 0 otherwise.

    // level 0 data vector size 36, with all entries 0 except for parent cells (1)
    expected_data_levels[0] = std::vector<double>(9,0);
    expected_data_levels[0][4] = 1; // parent index LGR1 = {4}

    // level 1 data vector of 1's, size 12 (cells_per_dim = {3,2,2}, block parent cells dim = {1,1,1})
    expected_data_levels[1] = std::vector<double>(12,1);

    // level 2 data vector of 1's, size 8 (cells_per_dim = {2,2,2}, block parent cells dim = {1,1,1})
    expected_data_levels[2] = std::vector<double>(8,1);

    // level 3 data vector of 1's, size 4 (cells_per_dim = {2,2,1}, block parent cells dim = {1,1,1})
    expected_data_levels[3] = std::vector<double>(4,1);

    // level 4 data vector of 1's, size 24 (cells_per_dim = {3,4,2}, block parent cells dim = {1,1,1})
    expected_data_levels[4] = std::vector<double>(24,1);

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
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

    std::vector<std::vector<double>> expected_data_levels{};
    expected_data_levels.resize(grid.maxLevel()+1);

    // Expected data for level grids is 1 for parent and children cells, 0 otherwise.
    
    // level 0 data vector size 36, with all entries 0 except for parent cells (1)
    expected_data_levels[0] = std::vector<double>(9,0);
    expected_data_levels[0][4] = 1; // parent index LGR1 = {4}
    expected_data_levels[0][0] = 1; // parent index LGR3 = {0}
    /** level 1 contains LGR1, whose parent grid is GLOBAL
        level 2 contains LGR3, whose parent grid is GLOBAL
        level 3 contains LGR2, whose parent grid is LGR1
        level 4 contains LGR4, whose parent grid is LGR3   */

    // level 1 (LGR1) data vector of 1's, size 12 (cells_per_dim = {3,2,2}, block parent cells dim = {1,1,1})
    expected_data_levels[1] = std::vector<double>(12,1);

    // level 2 (LGR3) data vector of 1's, size 4 (cells_per_dim = {2,2,1}, block parent cells dim = {1,1,1})
    expected_data_levels[2] = std::vector<double>(4,1);

    // level 3 (LGR2) data vector of 1's, size 8 (cells_per_dim = {2,2,2}, block parent cells dim = {1,1,1})
    expected_data_levels[3] = std::vector<double>(8,1);

    // level 4 (LGR4) data vector of 1's, size 24 (cells_per_dim = {3,4,2}, block parent cells dim = {1,1,1})
    expected_data_levels[4] = std::vector<double>(24,1);

    restrictFakeLeafDataToLevelGrids(grid, expected_data_levels);
}
