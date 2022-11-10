//===========================================================================
//
// File: CpGrid.hpp
//
// Created: Fri May 29 20:26:36 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            B�rd Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010, 2014, 2022 Equinor ASA.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulartion-Software & Services
  Copyright 2015       NTNU

  This file is part of The Open Porous Media project  (OPM).

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

#ifndef OPM_CPGRID_HEADER
#define OPM_CPGRID_HEADER

#include <string>
#include <map>
#include <array>
#include <unordered_set>
#include <opm/grid/utility/ErrorMacros.hpp>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/version.hh>

#if HAVE_MPI
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
#include <dune/common/parallel/variablesizecommunicator.hh>
#else
#include <opm/grid/utility/VariableSizeCommunicator.hpp>
#endif
#endif

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/gridenums.hh>

#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include "cpgrid/Intersection.hpp"
#include "cpgrid/Entity.hpp"
#include "cpgrid/Geometry.hpp"
#include "cpgrid/CpGridData.hpp"
#include "cpgrid/Iterators.hpp"
#include "cpgrid/Indexsets.hpp"
#include "cpgrid/DefaultGeometryPolicy.hpp"
#include "common/GridEnums.hpp"
#include "common/Volumes.hpp"
#include <opm/grid/cpgpreprocess/preprocess.h>

#include <opm/grid/utility/OpmParserIncludes.hpp>

#include <iostream>
#if ! HAVE_MPI
#include <list>
#endif

namespace Dune
{

    class CpGrid;

    namespace cpgrid
    {
        class CpGridData;
    }

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridTraits
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridTraits
    {
        /// \brief The type that implements the grid.
        typedef CpGrid Grid;

        /// \brief The type of the intersection at the leafs of the grid.
        typedef cpgrid::Intersection LeafIntersection;
        /// \brief The type of the intersection at the levels of the grid.
        typedef cpgrid::Intersection LevelIntersection;
        /// \brief The type of the intersection iterator at the leafs of the grid.
        typedef cpgrid::IntersectionIterator LeafIntersectionIterator;
        /// \brief The type of the intersection iterator at the levels of the grid.
        typedef cpgrid::IntersectionIterator LevelIntersectionIterator;

        /// \brief The type of the  hierarchic iterator.
        typedef cpgrid::HierarchicIterator HierarchicIterator;

        /// \brief Traits associated with a specific codim.
        /// \tparam cd The codimension.
        template <int cd>
        struct Codim
        {
            /// \brief The type of the geometry associated with the entity.
            /// IMPORTANT: Codim<codim>::Geometry == Geometry<dim-codim,dimw>
             typedef cpgrid::Geometry<3-cd, 3> Geometry;
            //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> Geometry;
            /// \brief The type of the local geometry associated with the entity.
             typedef cpgrid::Geometry<3-cd, 3> LocalGeometry;
            //typedef Dune::Geometry<3-cd, 3, CpGrid, cpgrid::Geometry> LocalGeometry;
            /// \brief The type of the entity.
            typedef cpgrid::Entity<cd> Entity;

            /// \brief The type of the iterator over all level entities of this codim.
            typedef cpgrid::Iterator<cd, All_Partition> LevelIterator;

            /// \brief The type of the iterator over all leaf entities of this codim.
            typedef cpgrid::Iterator<cd, All_Partition> LeafIterator;

            /// \brief The type of the entity pointer for entities of this codim.
            typedef cpgrid::Entity<cd> EntitySeed;

            /// \brief Traits associated with a specific grid partition type.
            /// \tparam pitype The type of the grid partition.
            template <PartitionIteratorType pitype>
            struct Partition
            {
                /// \brief The type of the iterator over the level entities of this codim on this partition.
                typedef cpgrid::Iterator<cd, pitype> LevelIterator;
                /// \brief The type of the iterator over the leaf entities of this codim on this partition.
                typedef cpgrid::Iterator<cd, pitype> LeafIterator;
            };
        };

        /// \brief Traits associated with a specific grid partition type.
        /// \tparam pitype The type of the grid partition.
        template <PartitionIteratorType pitype>
        struct Partition
        {
            /// \brief The type of the level grid view associated with this partition type.
            typedef Dune::GridView<DefaultLevelGridViewTraits<CpGrid> > LevelGridView;
            /// \brief The type of the leaf grid view associated with this partition type.
            typedef Dune::GridView<DefaultLeafGridViewTraits<CpGrid> > LeafGridView;

        };

        /// \brief The type of the level grid view associated with this partition type.
        typedef Dune::GridView<DefaultLevelGridViewTraits<CpGrid> > LevelGridView;
        /// \brief The type of the leaf grid view associated with this partition type.
        typedef Dune::GridView<DefaultLeafGridViewTraits<CpGrid> > LeafGridView;

        /// \brief The type of the level index set.
        typedef cpgrid::IndexSet LevelIndexSet;
        /// \brief The type of the leaf index set.
        typedef cpgrid::IndexSet LeafIndexSet;
        /// \brief The type of the global id set.
        typedef cpgrid::GlobalIdSet GlobalIdSet;
        /// \brief The type of the local id set.
        typedef GlobalIdSet LocalIdSet;

        /// \brief The type of the collective communication.

    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
    using Communication = Dune::Communication<MPICommunicator>;
    using CollectiveCommunication = Communication;
#else
    using CollectiveCommunication = Dune::CollectiveCommunication<MPICommunicator>;
    using Communication = Dune::CollectiveCommunication<MPICommunicator>;
#endif
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGridFamily
    //
    ////////////////////////////////////////////////////////////////////////

    struct CpGridFamily
    {
        typedef CpGridTraits Traits;
    };

    ////////////////////////////////////////////////////////////////////////
    //
    //   CpGrid
    //
    ////////////////////////////////////////////////////////////////////////

    /// \brief [<em> provides \ref Dune::Grid </em>]
    class CpGrid
        : public GridDefaultImplementation<3, 3, double, CpGridFamily >
    {
        friend class cpgrid::CpGridData;
        template<int dim>
        friend cpgrid::Entity<dim> createEntity(const CpGrid&,int,bool);
        friend
        void ::refine_and_check(const Dune::cpgrid::Geometry<3, 3>&,
                                const std::array<int, 3>&,
                                bool);

    public:

        // --- Typedefs ---


        /// Family typedef, why is this not defined by Grid<>?
        typedef CpGridFamily GridFamily;


        // --- Methods ---


        /// Default constructor
        CpGrid();

        CpGrid(MPIHelper::MPICommunicator comm);

        /// \name IO routines
        //@{
        /// Read the Sintef legacy grid format ('topogeom').
        /// \param grid_prefix the grid name, such that topology is
        /// found in <grid_prefix>-topo.dat etc.
        void readSintefLegacyFormat(const std::string& grid_prefix);


        /// Write the Sintef legacy grid format ('topogeom').
        /// \param grid_prefix the grid name, such that topology will be
        /// found in <grid_prefix>-topo.dat etc.
        void writeSintefLegacyFormat(const std::string& grid_prefix) const;


#if HAVE_ECL_INPUT
        /// Read the Eclipse grid format ('grdecl').
        ///
        /// \return Vector of global indices to the cells which have
        ///         been removed in the grid processing due to small pore volume. Function only returns
        ///         indices on rank 0, the vector is empty of other ranks.
        /// \param ecl_grid the high-level object from opm-parser which represents the simulation's grid
        ///        In a parallel run this may be a nullptr on all ranks but rank zero.
        /// \param ecl_state the object from opm-parser provide information regarding to pore volume, NNC,
        ///        aquifer information when ecl_state is available. NNC and aquifer connection
        ///        information will also be updated during the function call when available and necessary.
        /// \param periodic_extension if true, the grid will be (possibly) refined, so that
        ///        intersections/faces along i and j boundaries will match those on the other
        ///        side. That is, i- faces will match i+ faces etc.
        /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
        /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
        /// \param pinchActive Force specific pinch behaviour. If true a face will connect two vertical cells, that are
        ///           topological connected, even if there are cells with zero volume between them. If false these
        ///           cells will not be connected despite their faces coinciding.
        std::vector<std::size_t> processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                                                      Opm::EclipseState* ecl_state,
                                                      bool periodic_extension, bool turn_normals, bool clip_z,
                                                      bool pinchActive);

        /// Read the Eclipse grid format ('grdecl').
        ///
        /// Pinch behaviour is determind from the parameter ecl_grid. If ecl_grid is a nullptr or PINCH was specified for
        /// the grid, then a face will connect two vertical cells, that are topological connected, even if there are
        /// cells with zero volume between them, Otherwise these cells will not be connected despite their faces coinciding.
        ///
        /// \return Vector of global indices to the cells which have
        ///         been removed in the grid processing due to small pore volume. Function only returns
        ///         indices on rank 0, the vector is empty of other ranks.
        /// \param ecl_grid the high-level object from opm-parser which represents the simulation's grid
        ///        In a parallel run this may be a nullptr on all ranks but rank zero.
        /// \param ecl_state the object from opm-parser provide information regarding to pore volume, NNC,
        ///        aquifer information when ecl_state is available. NNC and aquifer connection
        ///        information will also be updated during the function call when available and necessary.
        /// \param periodic_extension if true, the grid will be (possibly) refined, so that
        ///        intersections/faces along i and j boundaries will match those on the other
        ///        side. That is, i- faces will match i+ faces etc.
        /// \param turn_normals if true, all normals will be turned. This is intended for handling inputs with wrong orientations.
        /// \param clip_z if true, the grid will be clipped so that the top and bottom will be planar.
        std::vector<std::size_t> processEclipseFormat(const Opm::EclipseGrid* ecl_grid,
                                                      Opm::EclipseState* ecl_state,
                                                      bool periodic_extension, bool turn_normals = false, bool clip_z = false);

#endif

        /// Read the Eclipse grid format ('grdecl').
        /// \param input_data the data in grdecl format, declared in preprocess.h.
        /// \param remove_ij_boundary if true, will remove (i, j) boundaries. Used internally.
        void processEclipseFormat(const grdecl& input_data, bool remove_ij_boundary, bool turn_normals = false);

        //@}

        /// \name Cartesian grid extensions.
        ///
        /// A cornerpoint grid can be seen as a degenerated and distorted cartesian
        /// grid. Therefore it provides mappings from cells to the underlying cartesian
        /// index.
        //@{
        /// Create a cartesian grid.
        /// \param dims the number of cells in each cartesian direction.
        /// \param cellsize the size of each cell in each dimension.
        void createCartesian(const std::array<int, 3>& dims,
                             const std::array<double, 3>& cellsize);

        /// The logical cartesian size of the global grid.
        /// This function is not part of the Dune grid interface,
        /// and should be used with caution.
        const std::array<int, 3>& logicalCartesianSize() const
        {
            return current_view_data_->logical_cartesian_size_;
        }

        /// Retrieve mapping from internal ("compressed") active grid
        /// cells to external ("uncompressed") cells.  Specifically,
        /// @code globalCell()[i] @endcode is the linearized Cartesian
        /// index of grid cell @code i @endcode.  This method should
        /// only be used by classes which really need it, such as
        /// those dealing with permeability fields from the input deck
        /// from whence the current CpGrid was constructed.
        const std::vector<int>& globalCell() const
        {
            return current_view_data_->global_cell_;
        }

        /// @brief
        ///    Extract Cartesian index triplet (i,j,k) of an active cell.
        ///
        /// @param [in] c
        ///    Active cell index.
        ///
        /// @param [out] ijk  Cartesian index triplet
        void getIJK(const int c, std::array<int,3>& ijk) const
        {
            current_view_data_->getIJK(c, ijk);
        }
        //@}

        /// Is the grid currently using unique boundary ids?
        /// \return true if each boundary intersection has a unique id
        ///         false if we use the (default) 1-6 ids for i- i+ j- j+ k- k+ boundaries.
        bool uniqueBoundaryIds() const
        {
            return current_view_data_->uniqueBoundaryIds();
        }

        /// Set whether we want to have unique boundary ids.
        /// \param uids if true, each boundary intersection will have a unique boundary id.
        void setUniqueBoundaryIds(bool uids)
        {
            current_view_data_->setUniqueBoundaryIds(uids);
        }

        // --- Dune interface below ---

        /// \name The DUNE grid interface implementation
        // \@{
        /// \brief Get the grid name.
        ///
        /// It's the same as the class name.
        /// What did you expect, something funny?
        std::string name() const
        {
            return "CpGrid";
        }


        /// Return maximum level defined in this grid. Levels are numbered
        /// 0 ... maxlevel with 0 the coarsest level.
        int maxLevel() const
        {
            return 0;
        }


        /// Iterator to first entity of given codim on level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
        }


        /// one past the end on this level
        template<int codim>
        typename Traits::template Codim<codim>::LevelIterator lend (int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,All_Partition>(*current_view_data_, size(codim), true );

        }


        /// Iterator to first entity of given codim on level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType>(*current_view_data_, 0, true );
        }


        /// one past the end on this level
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return cpgrid::Iterator<codim,PiType>(*current_view_data_, size(codim), true);
        }


        /// Iterator to first leaf entity of given codim
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafbegin() const
        {
            return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, 0, true);
        }


        /// one past the end of the sequence of leaf entities
        template<int codim>
        typename Traits::template Codim<codim>::LeafIterator leafend() const
        {
            return cpgrid::Iterator<codim, All_Partition>(*current_view_data_, size(codim), true);
        }


        /// Iterator to first leaf entity of given codim
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const
        {
            return cpgrid::Iterator<codim, PiType>(*current_view_data_, 0, true);
        }


        /// one past the end of the sequence of leaf entities
        template<int codim, PartitionIteratorType PiType>
        typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const
        {
            return cpgrid::Iterator<codim, PiType>(*current_view_data_, size(codim), true);
        }


        /// \brief Number of grid entities per level and codim
        int size (int level, int codim) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(codim);
        }


        /// number of leaf entities per codim in this process
        int size (int codim) const
        {
            return current_view_data_->size(codim);
        }


        /// number of entities per level and geometry type in this process
        int size (int level, GeometryType type) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return size(type);
        }


        /// number of leaf entities per geometry type in this process
        int size (GeometryType type) const
        {
            return current_view_data_->size(type);
        }


        /// \brief Access to the GlobalIdSet
        const Traits::GlobalIdSet& globalIdSet() const
        {
            return global_id_set_;
        }


        /// \brief Access to the LocalIdSet
        const Traits::LocalIdSet& localIdSet() const
        {
            return global_id_set_;
        }


        /// \brief Access to the LevelIndexSets
        const Traits::LevelIndexSet& levelIndexSet(int level) const
        {
            if (level<0 || level>maxLevel())
                DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
            return *current_view_data_->index_set_;
        }


        /// \brief Access to the LeafIndexSet
        const Traits::LeafIndexSet& leafIndexSet() const
        {
            return *current_view_data_->index_set_;
        }


        /// global refinement
        void globalRefine (int)
        {
            std::cout << "Warning: Global refinement not implemented, yet." << std::endl;
        }

        const std::vector< Dune :: GeometryType >& geomTypes( const int codim ) const
        {
          return leafIndexSet().geomTypes( codim );
        }

        /// given an EntitySeed (or EntityPointer) return an entity object
        template <int codim>
        cpgrid::Entity<codim> entity( const cpgrid::Entity< codim >& seed ) const
        {
            return seed;
        }

        // @TODO Given a CpGridData that's level 0 of a data object,
        // create "future_leaf_cells", "future_leaf_corners"
        //  Adding a 0: {0, cell index} , {0, corner index}.

        // ADD LEVEL, 'PREDICT' LEAF CELLS AND CORNERS (void; we add a new entry to 3 exisiting objects)
        // Take as references a vector with shared pointers of type CpGridData ("data") and 2 vectors with the
        // indices of cells/corners, coming possible from different levels, that all together would form the leaf cells/corners,
        // (kind of predicting the "future_leaf_cells/corners"). 
        // Construct an LGR choosing an existing entry of "data", given the amount of
        // children cells in each direction, and the begining and end of the patch to be refined.
        // Add this 'level' to "data" and its corresponding index-access information in "future_leaf_view".
        // @param data                          Vector of shared pointers of type CpGridData (each ptr ~ one level).
        // @param future_leaf_cells             Vector. Each entry looks like: {cell level, cell index}.
        //                                      cell level = which entry of "data" points at the CpGridData object
        //                                                   where the cell belongs.
        //                                      cell index -> to access the Geometry<3,3> via the entry of "data"
        //                                                    that points at the CpGridData object
        //                                                    where the cell belongs.
        // @param future_leaf_corners           Vector. Each entry looks like: {corner level, corner index}.
        //                                      corner level = which entry of "data" points at the CpGridData object
        //                                                   where the corner belongs.
        //                                      corner index -> to access the Geometry<0,3> via the entry of "data"
        //                                                    that points at the CpGridData object
        //                                                    where the corners belongs.
        // @param level_to_refine               Integer (smaller than data.size()) representing the level from where
        //                                      the patch to refine is taken.      
        // @param cells_per_dim                 Number of sub-cells in each direction (for each cell) in the lgr.
        // @param start_ijk                     Minimum values of i,j,k where the patch/lgr 'starts'.
        // @param end_ijk                       Maximum values of i,j,k where the patch/lgr 'ends'.
        //
        // Idea of "future_leaf_cells":
        // We start with a grid that plays the role of level 0, so then "future_leaf_cells" will contain all
        // the (active) cells of the grid, with a 0 to indicate they were born in level 0:
        // {{0, cell index}, {0, next cell index}, ..., {0, last cell index}}
        // For example, {{0,0}, {0,1}, ...,  {0,23}} for a grid with 4,3,2 cells in x,y,z-direction respectively.
        // When constructing an LGR out of the grid from level 0, we call it level 1 and will erase from "future_leaf_cells"
        // all the std::array<int,2> which correspond to parents of this patch. 
        // In the previous example, let's say our patch consists of the cells with indices 0,1,4,5
        // (strat_ijk = {0,0,0}, end_ijk = {2,2,1}). Then we earase from "future_leaf_cells" the entries
        // {0,0}, {0,1}, {0,4}, {0,5}, and 'replace' them [push_back] by (if, for instance, cells_per_dim = {5,6,7})
        // {1,0}, {1,1}, ..., {1, ([2x2x1]x[5x6x7]) -1}.
        // Then the vector "future_leaf_cells", after 'level 1'-refinement, looks like:
        // {{0,2}, {0,3}, {0,6}, {0,7}, ..., {0,23}, {1,0}, {1,1}, ..., {1, ([2x2x1]x[5x6x7]) -1}}.
        // Notice that 'the union of the cells stored in "future_leaf_cells"' would give us the leaf cells  IN THE END.
        // Every time we refine, we delete the parent cells and replace them with their children. 
        // Now, when we build level 2 we choose a patch in level 0 or level 1, we proceed in the same way.
        // First delete from "future_leaf_cells" the parent cells, then push back the new child cells with
        // their correspongind level.
        // To construct a leaf view, it 'will(might) be' enough to iterate on the entries of "future_leaf_cells".
        // Remark: in this code only a block patch belonging to one level can be refined,
        // namely, we cannot refine in the same LGR 'a cell from some_level' and 'a cell from some_other_level'.
        //
        // Idea of "future_leaf_corners":
        // In the same spitit as in "future_leaf_view", we want to delete in each refinement to avoid repetition.
        // Assume that before the first refinement, "future_leaf_corners" looks like:
        // {{0, corner 0}, {0, corner 1}, ..., {0, total_corners -1}}.
        // Let's say we want to refine the 'cell 0'. For a grid with 4x3x2 cells, this means that our 'cell 0'
        // has corners with indices { 0,3,15,18 (bottom) , 1,4,16,19 (top) } (fake{0,1,..,7}).
        // Then, we delete from "future_leaf_corners" the following entries:
        // {0,0}, {0,1}, {0,3}, {0,4}, {0,15}, {0,16}, {0,18}, {0,19}
        // and 'replace them' (push_back) the new corners arising from the first LGR ('level 1').
        void addLevel(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                      std::vector<std::array<int,2>>& future_leaf_cells,
                      std::vector<std::array<int,2>>& future_leaf_corners,
                      // @TODO ADD CONTAINER FOR 'FUTURE LEAF FACES' TO AVOID REPETITION
                      int level_to_refine,
                      const std::array<int,3>& cells_per_dim,
                      std::array<int,3> start_ijk, std::array<int,3> end_ijk)
        {
            // Use *data.back() instead, if we only want to allow refinement based on the last level stored in data
            auto [new_data_entry, parent_to_child, child_to_parent] =
                (*data[level_to_refine]).refine_block_patch(cells_per_dim, start_ijk, end_ijk);
            data.push_back(new_data_entry);
            // @todo Where (and why? is it needed?) to store "parent_to_child" and "chil_to_parent"
            auto [patch_dim, patch_cell_indices] =
                (*data[level_to_refine]).get_patch_dim_and_cellIndices(start_ijk, end_ijk);
            for (auto idx : patch_cell_indices) {
                std::array<int, 2> parent_to_erase = { level_to_refine, idx};
                auto parent_to_erase_it = std::find(future_leaf_cells.begin(),
                                                     future_leaf_cells.end(), parent_to_erase);
                future_leaf_cells.erase(parent_to_erase_it);
                }
            // Add the (all) child cells to "future_leaf_cells". We do not separate them by 'parent'.
            // Recall that the numbering for cells follows the rule:
            // from left to right (increasing i/x-direction),
            // front to back (increasing j/y-direction),
            // from bottom to top (increasing k/z-direction).
            for (int new_cell = 0; new_cell < patch_dim[0]*cells_per_dim[0]
                     *patch_dim[1]*cells_per_dim[1]*patch_dim[2]*cells_per_dim[2]; ++new_cell) {
                future_leaf_cells.push_back({data.size() +1, new_cell});
            }
            // We delete the corners of the coarse level involved in the patch,
            // to avoid repetition.
            // Get size of the coarse level where the pacth is (needed to compute corner indices to be deleted).
            std::array<int,3> coarse_level_dim = (*data[level_to_refine]).logicalCartesianSize();
            for (int j = start_ijk[1]; j < end_ijk[1] +1; ++j) {
                for (int i = start_ijk[0]; i < end_ijk[0] +1; ++i) { 
                    for (int k = start_ijk[2]; k < end_ijk[2] +1; ++k) {
                        // Index (in the coarse level) of the corner we want to delete 
                        std::array<int,2> corner_to_erase = {level_to_refine, (j*(coarse_level_dim[0]+1)*(coarse_level_dim[2]+1))
                            + (i*(coarse_level_dim[2]+1)) +k};
                        auto corner_to_erase_it = std::find(future_leaf_corners.begin(),
                                                            future_leaf_corners.end(),
                                                            corner_to_erase);
                        future_leaf_corners.erase(corner_to_erase_it);     
                    } // end k-for-loop
                } // end i-for-loop
            } // end j-for-loop
            // Get the size of the new corners
            int total_new_corners = ((cells_per_dim[0]*patch_dim[0]) + 1) *
                ((cells_per_dim[1]*patch_dim[1]) + 1) * ((cells_per_dim[2]*patch_dim[2]) + 1);
            // Add the (child) corners to "future_leaf_corners". We do not separate them by 'parent'.
            // Recall that the numbering for corners follows the rule:
            // from bottom to top - from left to right - from front to back.
            for (int new_corner = 0; new_corner < total_new_corners; ++new_corner) {
                future_leaf_corners.push_back({data.size() +1, new_corner});
            }
        }

        // GET LEAF VIEW -------------- IN PROGRESS
        typedef Dune::FieldVector<double,3> PointType;
        std::shared_ptr<Dune::cpgrid::CpGridData> getLeafView(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> data,
                                                              std::vector<std::array<int,2>> future_leaf_view,
                                                              std::vector<std::array<int,2>> future_leaf_corners) const
        {
            // To store the leaf view.
            std::shared_ptr<Dune::cpgrid::CpGridData> leaf_view_ptr = std::make_shared<Dune::cpgrid::CpGridData>(); // ccobj_
            auto& leaf_view = *leaf_view_ptr;
            Dune::cpgrid::DefaultGeometryPolicy& leaf_geometries = leaf_view.geometry_;
            std::vector<std::array<int,8>>& leaf_cell_to_point = leaf_view.cell_to_point_;
            cpgrid::OrientedEntityTable<0,1>& leaf_cell_to_face = leaf_view.cell_to_face_;
            Opm::SparseTable<int>& leaf_face_to_point = leaf_view.face_to_point_;
            cpgrid::OrientedEntityTable<1,0>& leaf_face_to_cell = leaf_view.face_to_cell_;
            cpgrid::EntityVariable<enum face_tag,1>& leaf_face_tags = leaf_view.face_tag_;
            cpgrid::SignedEntityVariable<Dune::FieldVector<double,3>,1>& leaf_face_normals = leaf_view.face_normals_;

            Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<0,3>>& leaf_corners =
                leaf_geometries.geomVector(std::integral_constant<int,3>());
            Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<2,3>>& leaf_faces =
                leaf_geometries.geomVector(std::integral_constant<int,1>());
            Dune::cpgrid::EntityVariableBase<cpgrid::Geometry<3,3>>& leaf_cells =
                leaf_geometries.geomVector(std::integral_constant<int,0>());
            Dune::cpgrid::EntityVariableBase<enum face_tag>& mutable_face_tags = leaf_face_tags;
            Dune::cpgrid::EntityVariableBase<PointType>& mutable_face_normals = leaf_face_normals;

            
//leaf_cells/faces/corners -> resize  and then []
//leaf_cell_to_point -> resize and then []. QUESTION: SHOULD IT BE appendRow and std::vector<cpgrid::EntityRep<0>???
//                                         QUESTION: IS THE TYPE CORRECT? 
//                                       Current type (also in Geometry.hpp): std::vector<std::array<int,8>>
//                                          Should it be  Opm::SparseTable<int> (as in *_face_to_point)?
//leaf_cell_to_face -> appendRow( bla.begin(), bla.end()), bla type std::vector<cpgrid::EntityRep<1>>
//leaf_face_to_point -> appendRow(bla.begin(), bla.end()) bla type std::array<int,4>
//leaf_face_to_cell -> appendRow(bla.begin(), bla.end()) bla typ std::vector<cpgrid::EntityRep<0>>

            // Get the leaf corners
            // Determine the size
            int total_corners = future_leaf_corners.size();
            leaf_corners.reserve(total_corners);
            // The order of the leaf corners follows the rule of
            // 0. all the ('survivor') corners of level 0 (they were never involved in any refinement),
            //   with the order from bottom to top-from left to right- from front to back,
            // 1. all the ('survivor') corners of level 1 (not involved in the rest of the refinement process)
            //   with the order from bottom to top-from left to right- from front to back, and so on...
            for (auto idx : future_leaf_corners) {
                leaf_corners.push_back((*data[idx[0]]).geometry_.geomVector(std::integral_constant<int,3>())
                                       [Dune::cpgrid::EntityRep<3>(idx[1], true)]);
            }

            // Get total cells in the leaf view (possibly coming from different levels).
            int total_cells = future_leaf_view.size();
            int leaf_cell_idx = 0;
            leaf_cells.resize(total_cells);
            leaf_cell_to_point.resize(total_cells);
            // To populate "leaf_cells"(and "leaf_cell_to_*") we follow the rule:
            // 0. all the cells from level 0 that haven't been refined in the entire process.
            // 1. all the cells from level 1 that haven't been refined in the rest of the process. [and so on ...]
            // This is exactly how the cells are ordered in "future_leaf_view"
            // (thanks to the deleting-pushing_back constuction).
            for (auto idx : future_leaf_view) {
                // idx = {level the cell was born in, index cell in that level}.
                // idx[0] = the level where the cell was born = the entry of "data" we need to access to this cell.
                leaf_cells[leaf_cell_idx] = (*data[idx[0]]).geometry_.geomVector(std::integral_constant<int,0>())
                                     [Dune::cpgrid::EntityRep<0>(idx[1], true)];
                leaf_cell_to_point[leaf_cell_idx] =  (*data[idx[0]]).cell_to_point_[idx[1]];
                leaf_cell_to_face.appendRow((*data[idx[0]]).cell_to_face_[idx[1]].begin(), (*data[idx[0]]).cell_to_face_[idx[1]].end());
                leaf_cell_idx +=1;
            }
    
            return leaf_view_ptr;
           }


        
        /*  No refinement implemented. GridDefaultImplementation's methods will be used.

        /// \brief Mark entity for refinement
        ///
        /// This only works for entities of codim 0.
        /// The parameter is currently ignored
        ///
        /// \return <ul>
        /// <li> true, if marking was succesfull </li>
        /// <li> false, if marking was not possible </li>
        /// </ul>

        bool mark(int refCount, const typename Traits::template Codim<0>::EntityPointer & e)
        {
            return hostgrid_->mark(refCount, getHostEntity<0>(*e));
        }

        /// \brief Return refinement mark for entity
        ///
        /// \return refinement mark (1,0,-1)

        int getMark(const typename Traits::template Codim<0>::EntityPointer & e) const
        {
            return hostgrid_->getMark(getHostEntity<0>(*e));
        }

        /// \todo Please doc me !
        bool preAdapt() {
            return hostgrid_->preAdapt();
        }


        /// Triggers the grid refinement process
        bool adapt()
        {
            return hostgrid_->adapt();
        }

        /// \brief Clean up refinement markers
        void postAdapt() {
            return hostgrid_->postAdapt();
        }

        end of refinement section */

        /// \brief Size of the overlap on the leaf level
        unsigned int overlapSize(int) const {
            return 1;
        }


        /// \brief Size of the ghost cell layer on the leaf level
        unsigned int ghostSize(int) const {
            return 0;
        }


        /// \brief Size of the overlap on a given level
        unsigned int overlapSize(int, int) const {
            return 1;
        }


        /// \brief Size of the ghost cell layer on a given level
        unsigned int ghostSize(int, int) const {
            return 0;
        }

        /// \brief returns the number of boundary segments within the macro grid
        unsigned int numBoundarySegments() const
        {
            if( uniqueBoundaryIds() )
            {
                return current_view_data_->unique_boundary_ids_.size();
            }
            else
            {
                unsigned int numBndSegs = 0;
                const int num_faces = numFaces();
                for (int i = 0; i < num_faces; ++i) {
                    cpgrid::EntityRep<1> face(i, true);
                    if (current_view_data_->face_to_cell_[face].size() == 1) {
                        ++numBndSegs;
                    }
                }
                return numBndSegs;
            }
        }

        void setZoltanParams(const std::map<std::string,std::string>& params)
        {
          zoltanParams = params;
        }

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \warning May only be called once.
        bool loadBalance(int overlapLayers=1, bool useZoltan=true)
        {
            using std::get;
            return get<0>(scatterGrid(defaultTransEdgeWgt, false, nullptr, false, nullptr, true, overlapLayers, useZoltan ));
        }

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        ///
        /// This will construct the corresponding graph to the grid and use the transmissibilities
        /// specified as weights associated with its edges. The graph will be passed to the load balancer.
        /// \param wells The wells of the eclipse If null wells will be neglected.
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param transmissibilities The transmissibilities used as the edge weights.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(const std::vector<cpgrid::OpmWellType> * wells,
                    const double* transmissibilities = nullptr,
                    int overlapLayers=1, bool useZoltan=true)
        {
            return scatterGrid(defaultTransEdgeWgt, false, wells, false, transmissibilities, false, overlapLayers, useZoltan);
        }

        // loadbalance is not part of the grid interface therefore we skip it.

        /// \brief Distributes this grid over the available nodes in a distributed machine
        ///
        /// This will construct the corresponding graph to the grid and use the transmissibilities
        /// specified to calculate the  weights associated with its edges. The graph will be passed
        ///  to the load balancer.
        /// \param method The edge-weighting method to be used on the Zoltan partitioner.
        /// \param wells The wells of the eclipse If null wells will be neglected.
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(EdgeWeightMethod method, const std::vector<cpgrid::OpmWellType> * wells,
                    const double* transmissibilities = nullptr, bool ownersFirst=false,
                    bool addCornerCells=false, int overlapLayers=1,
                    bool useZoltan = true)
        {
            return scatterGrid(method, ownersFirst, wells, false, transmissibilities, addCornerCells, overlapLayers, useZoltan);
        }

        /// \brief Distributes this grid and data over the available nodes in a distributed machine.
        /// \param data A data handle describing how to distribute attached data.
        /// \param wells The wells of the eclipse  Default: null
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights.
        /// \param overlapLayers The number of layers of overlap cells to be added
        ///        (default: 1)
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        template<class DataHandle>
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(DataHandle& data,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    const double* transmissibilities = nullptr,
                    int overlapLayers=1, bool useZoltan = true)
        {
            auto ret = loadBalance(wells, transmissibilities, overlapLayers, useZoltan);
            using std::get;
            if (get<0>(ret))
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid over the available nodes in a distributed machine
        ///
        /// This will construct the corresponding graph to the grid and use the transmissibilities
        /// specified to calculate the  weights associated with its edges. The graph will be passed
        ///  to the load balancer.
        /// \param data A data handle describing how to distribute attached data.
        /// \param method The edge-weighting method to be used on the Zoltan partitioner.
        /// \param wells The information about all possible wells. If null then
        ///            the wells will be neglected. Otherwise the wells will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This is done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param serialPartitioning If true, the partitioning will be done on a single process.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \param zoltanImbalanceTol Set the imbalance tolerance used by Zoltan
        /// \param allowDistributedWells Allow the perforation of a well to be distributed to the
        ///        interior region of multiple processes.
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        template<class DataHandle>
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(DataHandle& data, EdgeWeightMethod method,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    bool serialPartitioning,
                    const double* transmissibilities = nullptr, bool ownersFirst=false,
                    bool addCornerCells=false, int overlapLayers=1, bool useZoltan = true,
                    double zoltanImbalanceTol = 1.1,
                    bool allowDistributedWells = false)
        {
            auto ret = scatterGrid(method, ownersFirst, wells, serialPartitioning, transmissibilities,
                                   addCornerCells, overlapLayers, useZoltan, zoltanImbalanceTol, allowDistributedWells);
            using std::get;
            if (get<0>(ret))
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param data A data handle describing how to distribute attached data.
        /// \param parts The partitioning information. For a cell with local index i the entry
        ///              parts[i] is the partion number. Partition numbers need to start with zero
        ///              and need to be consectutive also parts.size()==grid.leafGridView().size()
        ///              and the ranks communicator need to be able to map all parts. Needs to valid
        ///              at rank 0. Number of parts cannot exceed the number of ranks. Parts need to
        ///              numbered consecutively starting from zero.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        template<class DataHandle>
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        loadBalance(DataHandle& data, const std::vector<int>& parts,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    bool ownersFirst=false,
                    bool addCornerCells=false, int overlapLayers=1)
        {
            using std::get;
            auto ret = scatterGrid(defaultTransEdgeWgt,  ownersFirst, wells,
                                   /* serialPartitioning = */ false,
                                   /* transmissibilities = */ {},
                                   addCornerCells, overlapLayers, /* useZoltan =*/ false,
                                   /* zoltanImbalanceTol (ignored) = */ 0.0,
                                   /* allowDistributedWells = */ true, parts);
            using std::get;
            if (get<0>(ret))
            {
                scatterData(data);
            }
            return ret;
        }
        /// \brief Distributes this grid and data over the available nodes in a distributed machine.
        /// \param data A data handle describing how to distribute attached data.
        /// \param overlapLayers The number of layers of overlap cells to be added
        ///        (default: 1)
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \tparam DataHandle The type implementing DUNE's DataHandle interface.
        /// \warning May only be called once.
        template<class DataHandle>
        bool loadBalance(DataHandle& data,
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
                         decltype(data.fixedSize(0,0)) overlapLayers=1, bool useZoltan = true)
#else
                         decltype(data.fixedsize(0,0)) overlapLayers=1, bool useZoltan = true)
#endif
        {
            // decltype usage needed to tell the compiler not to use this function if first
            // argument is std::vector but rather loadbalance by parts
            bool ret = loadBalance(overlapLayers, useZoltan);
            if (ret)
            {
                scatterData(data);
            }
            return ret;
        }

        /// \brief Distributes this grid over the available nodes in a distributed machine
        /// \param parts The partitioning information. For a cell with local index i the entry
        ///              parts[i] is the partion number. Partition numbers need to start with zero
        ///              and need to be consectutive also parts.size()==grid.leafGridView().size()
        ///              and the ranks communicator need to be able to map all parts. Needs to valid
        ///              at rank 0. Number of parts cannot exceed the number of ranks. Parts need to
        ///              numbered consecutively starting from zero.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \warning May only be called once.
        bool loadBalance(const std::vector<int>& parts, bool ownersFirst=false,
                         bool addCornerCells=false, int overlapLayers=1)
        {
            using std::get;
            return get<0>(scatterGrid(defaultTransEdgeWgt,  ownersFirst, /* wells = */ {},
                                      /* serialPartitioning = */ false,
                                      /* trabsmissibilities = */ {},
                                      addCornerCells, overlapLayers, /* useZoltan =*/ false,
                                      /* zoltanImbalanceTol (ignored) = */ 0.0,
                                      /* allowDistributedWells = */ true, parts));
        }

        /// \brief Distributes this grid and data over the available nodes in a distributed machine
        /// \param data A data handle describing how to distribute attached data.
        /// \param parts The partitioning information. For a cell with local index i the entry
        ///              parts[i] is the partion number. Partition numbers need to start with zero
        ///              and need to be consectutive also parts.size()==grid.leafGridView().size()
        ///              and the ranks communicator need to be able to map all parts. Needs to valid
        ///              at rank 0. Number of parts cannot exceed the number of ranks. Parts need to
        ///              numbered consecutively starting from zero.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region (default: 1).
        /// \warning May only be called once.
        template<class DataHandle>
        bool loadBalance(DataHandle& data, const std::vector<int>& parts, bool ownersFirst=false,
                         bool addCornerCells=false, int overlapLayers=1)
        {
            bool ret = loadBalance(parts, ownersFirst, addCornerCells, overlapLayers);
            if (ret)
            {
                scatterData(data);
            }
            return ret;
        }

         /// \brief Partitions the grid using Zoltan without decomposing and distributing it among processes.
         /// \param wells The wells of the eclipse.
         /// \param transmissibilities The transmissibilities used to calculate the edge weights.
         /// \param numParts Number of parts in the partition.
         /// \return An array with the domain index for each cell.
         std::vector<int> zoltanPartitionWithoutScatter(const std::vector<cpgrid::OpmWellType> * wells,
                                                        const double* transmissibilities, int numParts,
                                                        const double zoltanImbalanceTol);

        /// The new communication interface.
        /// \brief communicate objects for all codims on a given level
        /// \param data The data handle describing the data. Has to adhere to the
        /// Dune::DataHandleIF interface.
        /// \param iftype The interface to use for the communication.
        /// \param dir The direction of the communication along the interface (forward or backward).
        /// \param level discarded as CpGrid is not adaptive.
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir, int /*level*/) const
        {
            communicate(data, iftype, dir);
        }

        /// The new communication interface.
        /// \brief communicate objects for all codims on a given level.
        /// \tparam DataHandle The type of the data handle describing the data.
        /// \param data The data handle describing the data. Has to adhere to the Dune::DataHandleIF interface.
        /// \param iftype The interface to use for the communication.
        /// \param dir The direction of the communication along the interface (forward or backward).
        /// \param level discarded as CpGrid is not adaptive.
        template<class DataHandle>
        void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
        {
            current_view_data_->communicate(data, iftype, dir);
        }

        /// \brief Get the collective communication object.
        const typename CpGridTraits::Communication& comm () const
        {
            return current_view_data_->ccobj_;
        }
        //@}

        // ------------ End of Dune interface, start of simplified interface --------------

        /// \name The simplified grid interface.
        ///
        /// It provides additional methods not in the DUNE interface but needed by OPM.
        /// Vertices, faces, and cells are not represented as entities but identified by
        /// indices.
        //@{
        // enum { dimension = 3 }; // already defined

        typedef Dune::FieldVector<double, 3> Vector;


        const std::vector<double>& zcornData() const {
            return current_view_data_->zcornData();
        }


        // Topology
        /// \brief Get the number of cells.
        int numCells() const
        {
            return current_view_data_->cell_to_face_.size();
        }
        /// \brief Get the number of faces.
        int numFaces() const
        {
            return current_view_data_->face_to_cell_.size();
        }
        /// \brief Get The number of vertices.
        int numVertices() const
        {
            return current_view_data_->geomVector<3>().size();
        }
        /// \brief Get the number of faces of a cell.
        ///
        /// Due to faults, and collapsing vertices (along pillars) this
        /// number is quite arbitrary. Its lower bound is 4, but there is
        /// no upper bound.
        /// \parame cell the index identifying the cell.
        int numCellFaces(int cell) const
        {
            return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)].size();
        }
        /// \brief Get a specific face of a cell.
        /// \param cell The index identifying the cell.
        /// \param local_index The local index (in [0,numFaces(cell))) of the face in this cell.
        /// \return The index identifying the face.
        int cellFace(int cell, int local_index) const
        {
            return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)][local_index].index();
        }

        /// \brief Get a list of indices identifying all faces of a cell.
        /// \param cell The index identifying the cell.
        const cpgrid::OrientedEntityTable<0,1>::row_type cellFaceRow(int cell) const
        {
            return current_view_data_->cell_to_face_[cpgrid::EntityRep<0>(cell, true)];
        }
        /// \brief Get the index identifying a cell attached to a face.
        ///
        /// Note that a face here is always oriented. If there are two
        /// neighboring cells then the orientation will be from local_index 0
        /// to local_index 1
        /// \param face The index identifying the face.
        /// \param local_index The local_index of the cell.
        /// \return The index identifying a cell or -1 if there is no such
        /// cell due the face being part of the grid boundary or the
        /// cell being stored on another process.
        int faceCell(int face, int local_index) const
        {
            // In the parallel case we store non-existent cells for faces along
            // the front region. Theses marked with index std::numeric_limits<int>::max(),
            // orientation might be arbitrary, though.
            cpgrid::OrientedEntityTable<1,0>::row_type r
                = current_view_data_->face_to_cell_[cpgrid::EntityRep<1>(face, true)];
            bool a = (local_index == 0);
            bool b = r[0].orientation();
            bool use_first = a ? b : !b;
            // The number of valid cells.
            int r_size = r.size();
            // In the case of only one valid cell, this is the index of it.
            int index = 0;
            if(r[0].index()==std::numeric_limits<int>::max()){
                assert(r_size==2);
                --r_size;
                index=1;
            }
            if(r.size()>1 && r[1].index()==std::numeric_limits<int>::max())
            {
                assert(r_size==2);
                --r_size;
            }
            if (r_size == 2) {
                return use_first ? r[0].index() : r[1].index();
            } else {
                return use_first ? r[index].index() : -1;
            }
        }
        /// \brief Get the sum of all faces attached to all cells.
        ///
        /// Each face identified by a unique index is counted as often
        /// as there are neigboring cells attached to it.
        /// \f$ numCellFaces()=\sum_{c} numCellFaces(c) \f$
        /// \see numCellFaces(int)const
        int numCellFaces() const
        {
            return current_view_data_->cell_to_face_.dataSize();
        }
        int numFaceVertices(int face) const
        {
            return current_view_data_->face_to_point_[face].size();
        }
        /// \brief Get the index identifying a vertex of a face.
        /// \param cell The index identifying the face.
        /// \param local_index The local_index (in [0,numFaceVertices(vertex) - 1]])
        ///  of the vertex.
        int faceVertex(int face, int local_index) const
        {
            return current_view_data_->face_to_point_[face][local_index];
        }
        /// \brief Get vertical position of cell center ("zcorn" average).
        /// \brief cell_index The index of the specific cell.
        double cellCenterDepth(int cell_index) const
        {
            // Here cell center depth is computed as a raw average of cell corner depths.
            // This generally gives slightly different results than using the cell centroid.
            double zz = 0.0;
            const int nv = current_view_data_->cell_to_point_[cell_index].size();
            const int nd = 3;
            for (int i=0; i<nv; ++i) {
                zz += vertexPosition(current_view_data_->cell_to_point_[cell_index][i])[nd-1];
            }
            return zz/nv;
        }

        const Vector faceCenterEcl(int cell_index, int face) const
        {
            // This method is an alternative to the method faceCentroid(...).
            // The face center is computed as a raw average of cell corners.
            // For faulted cells this gives different results then average of face nodes
            // that seems to agree more with eclipse.
            // This assumes the cell nodes are ordered
            // 6---7
            // | T |
            // 4---5
            //   2---3
            //   | B |
            //   0---1

            // this follows the DUNE reference cube
            static const int faceVxMap[ 6 ][ 4 ] = { {0, 2, 4, 6}, // face 0
                                                     {1, 3, 5, 7}, // face 1
                                                     {0, 1, 4, 5}, // face 2
                                                     {2, 3, 6, 7}, // face 3
                                                     {0, 1, 2, 3}, // face 4
                                                     {4, 5, 6, 7}  // face 5
                                                   };


            assert (current_view_data_->cell_to_point_[cell_index].size() == 8);
            Vector center(0.0);
            for( int i=0; i<4; ++i )
            {
               center += vertexPosition(current_view_data_->cell_to_point_[cell_index][ faceVxMap[ face ][ i ] ]);
            }

            for (int i=0; i<3; ++i) {
                center[i] /= 4;
            }
            return center;

        }

        const Vector faceAreaNormalEcl(int face) const
        {
            // same implementation as ResInsight
            const int nd = Vector::dimension;
            const int nv =  numFaceVertices(face);
            switch (nv)
            {
            case 0:
            case 1:
            case 2:
                {
                    return Vector(0.0);
                }
                break;
            case 3:
                {
                Vector a = vertexPosition(current_view_data_->face_to_point_[face][0]) - vertexPosition(current_view_data_->face_to_point_[face][2]);
                Vector b = vertexPosition(current_view_data_->face_to_point_[face][1]) - vertexPosition(current_view_data_->face_to_point_[face][2]);
                Vector areaNormal = cross(a,b);
                for (int i=0; i<nd; ++i) {
                    areaNormal[i] /= 2;
                }
                return areaNormal;
            }
                                break;
            case 4:
                {
                Vector a = vertexPosition(current_view_data_->face_to_point_[face][0]) - vertexPosition(current_view_data_->face_to_point_[face][2]);
                Vector b = vertexPosition(current_view_data_->face_to_point_[face][1]) - vertexPosition(current_view_data_->face_to_point_[face][3]);
                Vector areaNormal = cross(a,b);
                areaNormal *= 0.5;
                return areaNormal;
                }
                break;
            default:
                {
                    int h = (nv - 1)/2;
                    int k = (nv % 2) ? 0 : nv - 1;

                    Vector areaNormal(0.0);
                    // First quads
                    for (int i = 1; i < h; ++i)
                    {
                        Vector a = vertexPosition(current_view_data_->face_to_point_[face][2*i]) - vertexPosition(current_view_data_->face_to_point_[face][0]);
                        Vector b = vertexPosition(current_view_data_->face_to_point_[face][2*i+1]) - vertexPosition(current_view_data_->face_to_point_[face][2*i-1]);
                        areaNormal += cross(a,b);
                    }

                    // Last triangle or quad
                    Vector a = vertexPosition(current_view_data_->face_to_point_[face][2*h]) - vertexPosition(current_view_data_->face_to_point_[face][0]);
                    Vector b = vertexPosition(current_view_data_->face_to_point_[face][k]) - vertexPosition(current_view_data_->face_to_point_[face][2*h-1]);
                    areaNormal += cross(a,b);

                    areaNormal *= 0.5;

                    return areaNormal;
                }

            }
        }

        // Geometry
        /// \brief Get the Position of a vertex.
        /// \param cell The index identifying the cell.
        /// \return The coordinates of the vertex.
        const Vector& vertexPosition(int vertex) const
        {
            return current_view_data_->geomVector<3>()[cpgrid::EntityRep<3>(vertex, true)].center();
        }
        /// \brief Get the area of a face.
        /// \param cell The index identifying the face.
        double faceArea(int face) const
        {
            return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].volume();
        }
        /// \brief Get the coordinates of the center of a face.
        /// \param cell The index identifying the face.
        const Vector& faceCentroid(int face) const
        {
            return current_view_data_->geomVector<1>()[cpgrid::EntityRep<1>(face, true)].center();
        }
        /// \brief Get the unit normal of a face.
        /// \param cell The index identifying the face.
        /// \see faceCell
        const Vector& faceNormal(int face) const
        {
            return current_view_data_->face_normals_.get(face);
        }
        /// \brief Get the volume of the cell.
        /// \param cell The index identifying the cell.
        double cellVolume(int cell) const
        {
            return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].volume();
        }
        /// \brief Get the coordinates of the center of a cell.
        /// \param cell The index identifying the face.
        const Vector& cellCentroid(int cell) const
        {
            return current_view_data_->geomVector<0>()[cpgrid::EntityRep<0>(cell, true)].center();
        }

        /// \brief An iterator over the centroids of the geometry of the entities.
        /// \tparam codim The co-dimension of the entities.
        template<int codim>
        class CentroidIterator
            : public RandomAccessIteratorFacade<CentroidIterator<codim>,
                                                FieldVector<double, 3>,
                                                const FieldVector<double, 3>&, int>
        {
        public:
            /// \brief The type of the iterator over the codim geometries.
            typedef typename std::vector<cpgrid::Geometry<3-codim, 3> >::const_iterator
            GeometryIterator;
            /// \brief Constructs a new iterator from an iterator over the geometries.
            /// \param iter The iterator of the geometry objects.
            CentroidIterator(GeometryIterator iter)
            : iter_(iter)
            {}

            const FieldVector<double, 3>& dereference() const
            {
                return iter_->center();
            }
            void increment()
            {
                ++iter_;
            }
            const FieldVector<double, 3>& elementAt(int n)
            {
                return iter_[n]->center();
            }
            void advance(int n)
            {
                iter_+=n;
            }
            void decrement()
            {
                --iter_;
            }
            int distanceTo(const CentroidIterator& o)
            {
                return o-iter_;
            }
            bool equals(const CentroidIterator& o) const
            {
                return o==iter_;
            }
        private:
            /// \brief The iterator over the underlying geometry objects.
            GeometryIterator iter_;
        };

        /// \brief Get an iterator over the cell centroids positioned at the first one.
        CentroidIterator<0> beginCellCentroids() const
        {
            return CentroidIterator<0>(current_view_data_->geomVector<0>().begin());
        }

        /// \brief Get an iterator over the face centroids positioned at the first one.
        CentroidIterator<1> beginFaceCentroids() const
        {
            return CentroidIterator<1>(current_view_data_->geomVector<1>().begin());
        }

        // Extra
        int boundaryId(int face) const
        {
            // Note that this relies on the following implementation detail:
            // The grid is always construct such that the faces where
            // orientation() returns true are oriented along the positive IJK
            // direction. Oriented means that the first cell attached to face
            // has the lower index.
            int ret = 0;
            cpgrid::EntityRep<1> f(face, true);
            if (current_view_data_->face_to_cell_[f].size() == 1) {
                if (current_view_data_->uniqueBoundaryIds()) {
                    // Use the unique boundary ids.
                    ret = current_view_data_->unique_boundary_ids_[f];
                } else {
                    // Use the face tag based ids, i.e. 1-6 for i-, i+, j-, j+, k-, k+.
                    const bool normal_is_in =
                        !(current_view_data_->face_to_cell_[f][0].orientation());
                    enum face_tag tag = current_view_data_->face_tag_[f];
                    switch (tag) {
                    case I_FACE:
                        //                   LEFT : RIGHT
                        ret = normal_is_in ? 1    : 2; // min(I) : max(I)
                        break;
                    case J_FACE:
                        //                   BACK : FRONT
                        ret = normal_is_in ? 3    : 4; // min(J) : max(J)
                        break;
                    case K_FACE:
                        // Note: TOP at min(K) as 'z' measures *depth*.
                        //                   TOP  : BOTTOM
                        ret = normal_is_in ? 5    : 6; // min(K) : max(K)
                        break;
                    case NNC_FACE:
                        // This should not be possible, as NNC "faces" always
                        // have two cell neighbours.
                        OPM_THROW(std::logic_error, "NNC face at boundary. This should never happen!");
                    }
                }
            }
            return ret;
        }

        /// \brief Get the cartesian tag associated with a face tag.
        ///
        /// The tag tells us in which direction the face would point
        /// in the underlying cartesian grid.
        /// \param An iterator that points to the face and was obtained
        /// by iterating over Opm::UgGridHelpers::cell2Faces(grid).
        template<class Cell2FacesRowIterator>
        int
        faceTag(const Cell2FacesRowIterator& cell_face) const
        {
            // Note that this relies on the following implementation detail:
            // The grid is always constructed such that the interior faces constructed
            // with orientation set to true are
            // oriented along the positive IJK direction. Oriented means that
            // the first cell attached to face has the lower index.
            // For faces along the boundary (only one cell, always  attached at index 0)
            // the orientation has to be determined by the orientation of the cell.
            // If it is true then in UnstructuredGrid it would be stored at index 0,
            // otherwise at index 1.
            const int cell = cell_face.getCellIndex();
            const int face = *cell_face;
            assert (0 <= cell);  assert (cell < numCells());
            assert (0 <= face);  assert (face < numFaces());

            typedef cpgrid::OrientedEntityTable<1,0>::row_type F2C;

            const cpgrid::EntityRep<1> f(face, true);
            const F2C&     f2c = current_view_data_->face_to_cell_[f];
            const face_tag tag = current_view_data_->face_tag_[f];

            assert ((f2c.size() == 1) || (f2c.size() == 2));

            int inside_cell = 0;

            if ( f2c.size() == 2 ) // Two cells => interior
            {
                if ( f2c[1].index() == cell )
                {
                    inside_cell = 1;
                }
            }
            const bool normal_is_in = ! f2c[inside_cell].orientation();

            switch (tag) {
            case I_FACE:
                //                    LEFT : RIGHT
                return normal_is_in ? 0    : 1; // min(I) : max(I)
            case J_FACE:
                //                    BACK : FRONT
                return normal_is_in ? 2    : 3; // min(J) : max(J)
            case K_FACE:
                // Note: TOP at min(K) as 'z' measures *depth*.
                //                    TOP  : BOTTOM
                return normal_is_in ? 4    : 5; // min(K) : max(K)
            case NNC_FACE:
                // For nnc faces we return the otherwise unused value -1.
                return -1;
            default:
                OPM_THROW(std::logic_error, "Unhandled face tag. This should never happen!");
            }
        }

        //@}

        // ------------ End of simplified interface --------------

        //------------- methods not in the DUNE grid interface.

        /// \name Parallel grid extensions.
        /// Methods extending the DUNE's parallel grid interface.
        /// These are basically for scattering/gathering data to/from
        /// distributed views.
        //@{
        ///
        /// \brief Moves data from the global (all data on process) view to the distributed view.
        ///
        /// This method does not do communication but assumes that the global grid
        /// is present on every process and simply copies data to the distributed view.
        /// \tparam DataHandle The type of the data handle describing the data and responsible for
        ///         gathering and scattering the data.
        /// \param handle The data handle describing the data and responsible for
        ///         gathering and scattering the data.
        template<class DataHandle>
        void scatterData(DataHandle& handle) const
        {
#if HAVE_MPI
            if(distributed_data_.empty())
                OPM_THROW(std::runtime_error, "Moving Data only allowed with a load balanced grid!");
            distributed_data_[0]->scatterData(handle, data_[0].get(), distributed_data_[0].get(), cellScatterGatherInterface(),
                                           pointScatterGatherInterface());
#else
            // Suppress warnings for unused argument.
            (void) handle;
#endif
        }

        ///
        /// \brief Moves data from the distributed view to the global (all data on process) view.
        /// \tparam DataHandle The type of the data handle describing the data and responsible for
        ///         gathering and scattering the data.
        /// \param handle The data handle describing the data and responsible for
        ///         gathering and scattering the data.
        template<class DataHandle>
        void gatherData(DataHandle& handle) const
        {
#if HAVE_MPI
            if(distributed_data_.empty())
                OPM_THROW(std::runtime_error, "Moving Data only allowed with a load balance grid!");
            distributed_data_[0]->gatherData(handle, data_[0].get(), distributed_data_[0].get());
#else
            // Suppress warnings for unused argument.
            (void) handle;
#endif
        }
#if HAVE_MPI
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
        /// \brief The type of the map describing communication interfaces.
        using InterfaceMap = VariableSizeCommunicator<>::InterfaceMap;
#else
        /// \brief The type of the map describing communication interfaces.
        using InterfaceMap = Opm::VariableSizeCommunicator<>::InterfaceMap;
#endif
#else
        // bogus definition for the non parallel type. VariableSizeCommunicator not
        // availabe

        /// \brief The type of the map describing communication interfaces.
        typedef std::map<int, std::list<int> > InterfaceMap;
#endif

        /// \brief Get an interface for gathering/scattering data attached to cells with communication.
        ///
        /// Scattering means sending data from the indices of the global grid on
        /// process 0 to the distributed grid on all ranks independent of the grid.
        /// Gathering is the other way around.
        /// The interface can be used with VariableSizeCommunicator and a custom
        /// index based data handle to scatter (forward direction of the communicator)
        /// and gather data (backward direction of the communicator).
        /// Here is a small example that prints the received values when scattering:
        /// \code
        /// struct Handle{
        ///   typedef int DataType;
        ///   const std::vector<int>& vals;
        ///   bool fixedsize() { return true; }
        ///   size_t size(std::size_t) { return 1; }
        ///   void gather(auto& B buf, size_t i)[ buf.write(vals[i]); }
        ///   void scatter(auto& B buf, size_t i, std::size_t) {
        ///     int val;
        ///     buf.read(val);
        ///     cout<<i<<": "<<val<<" "; }
        /// };
        ///
        /// Handle handle;
        /// handle.vals.resize(grid.size(0), -1);
        /// Dune::VariableSizeCommunicator<> comm(grid.comm(),
        ///                                       grid.cellScatterGatherInterface());
        /// comm.forward(handle);
        /// \endcode
        const InterfaceMap& cellScatterGatherInterface() const
        {
            return *cell_scatter_gather_interfaces_;
        }

        /// \brief Get an interface for gathering/scattering data attached to points with communication.
        /// \see cellScatterGatherInterface
        const InterfaceMap& pointScatterGatherInterface() const
        {
            return *point_scatter_gather_interfaces_;
        }

        /// \brief Switch to the global view.
        void switchToGlobalView()
        {
            current_view_data_=data_[0].get();
        }

        /// \brief Switch to the distributed view.
        void switchToDistributedView()
        {
            if (distributed_data_.empty())
                OPM_THROW(std::logic_error, "No distributed view available in grid");
            current_view_data_=distributed_data_[0].get();
        }
        //@}

#if HAVE_MPI
        /// \brief The type of the parallel index set
        typedef cpgrid::CpGridData::ParallelIndexSet ParallelIndexSet;
        /// \brief The type of the remote indices information
        typedef cpgrid::CpGridData::RemoteIndices RemoteIndices;

        /// \brief The type of the owner-overlap-copy communication
        using CommunicationType = cpgrid::CpGridData::CommunicationType;

        /// \brief Get the owner-overlap-copy communication for cells
        ///
        /// Suitable e.g. for parallel linear algebra used by CCFV
        const CommunicationType& cellCommunication() const
        {
            return current_view_data_->cellCommunication();
        }

        ParallelIndexSet& getCellIndexSet()
        {
            return current_view_data_->cellIndexSet();
        }

        RemoteIndices& getCellRemoteIndices()
        {
            return current_view_data_->cellRemoteIndices();
        }

        const ParallelIndexSet& getCellIndexSet() const
        {
            return current_view_data_->cellIndexSet();
        }

        const RemoteIndices& getCellRemoteIndices() const
        {
            return current_view_data_->cellRemoteIndices();
        }

#endif

        /// \brief Get sorted active cell indices of numerical aquifer
       const std::vector<int>& sortedNumAquiferCells() const
       {
           return current_view_data_->sortedNumAquiferCells();
       }

    private:
        /// \brief Scatter a global grid to all processors.
        /// \param method The edge-weighting method to be used on the Zoltan partitioner.
        /// \param ownersFirst Order owner cells before copy/overlap cells.
        /// \param wells The wells of the eclipse If null wells will be neglected.
        ///            If this is not null then complete well information of
        ///            of the last scheduler step of the eclipse state will be
        ///            used to make sure that all the possible completion cells
        ///            of each well are stored on one process. This done by
        ///            adding an edge with a very high edge weight for all
        ///            possible pairs of cells in the completion set of a well.
        /// \param transmissibilities The transmissibilities used to calculate the edge weights in
        ///                           the Zoltan partitioner. This is done to improve the numerical
        ///                           performance of the parallel preconditioner.
        /// \param addCornerCells Add corner cells to the overlap layer.
        /// \param overlapLayers The number of layers of cells of the overlap region.
        /// \param useZoltan Whether to use Zoltan for partitioning or our simple approach based on
        ///        rectangular partitioning the underlying cartesian grid.
        /// \param zoltanImbalanceTol Set the imbalance tolerance used by Zoltan
        /// \param allowDistributedWells Allow the perforation of a well to be distributed to the
        ///        interior region of multiple processes.
        /// \param cell_part When using an external loadbalancer the partition number for each cell.
        ///                  If empty or not specified we use internal load balancing.
        /// \return A pair consisting of a boolean indicating whether loadbalancing actually happened and
        ///         a vector containing a pair of name and a boolean, indicating whether this well has
        ///         perforated cells local to the process, for all wells (sorted by name)
        std::pair<bool, std::vector<std::pair<std::string,bool> > >
        scatterGrid(EdgeWeightMethod method,
                    bool ownersFirst,
                    const std::vector<cpgrid::OpmWellType> * wells,
                    bool serialPartitioning,
                    const double* transmissibilities,
                    bool addCornerCells,
                    int overlapLayers,
                    bool useZoltan = true,
                    double zoltanImbalanceTol = 1.1,
                    bool allowDistributedWells = true,
                    const std::vector<int>& input_cell_part = {});

        /** @brief The data stored in the grid.
         *
         * All the data of all grids are stored there and
         * calls are forwarded to relevant grid.*/
        std::vector<std::shared_ptr<cpgrid::CpGridData>> data_;
        /** @brief A pointer to data of the current View. */
        cpgrid::CpGridData* current_view_data_;
        /** @brief The data stored for the distributed grid. */
        std::vector<std::shared_ptr<cpgrid::CpGridData>> distributed_data_;
        /**
         * @brief Interface for scattering and gathering cell data.
         *
         * @warning Will only update owner cells
         */
        std::shared_ptr<InterfaceMap> cell_scatter_gather_interfaces_;
        /*
         * @brief Interface for scattering and gathering point data.
         *
         * @warning Will only update owner cells
         */
        std::shared_ptr<InterfaceMap> point_scatter_gather_interfaces_;
        /**
         * @brief The global id set (also used as local one).
         */
        cpgrid::GlobalIdSet global_id_set_;

        /**
         * @brief Zoltan partitioning parameters
         */
        std::map<std::string,std::string> zoltanParams;

    }; // end Class CpGrid



    namespace Capabilities
    {
        /// \todo Please doc me !
        template <>
        struct hasEntity<CpGrid, 0>
        {
            static const bool v = true;
        };

        /// \todo Please doc me !
        template <>
        struct hasEntity<CpGrid, 3>
        {
            static const bool v = true;
        };

        template<>
        struct canCommunicate<CpGrid,0>
        {
            static const bool v = true;
        };

        template<>
        struct canCommunicate<CpGrid,3>
        {
            static const bool v = true;
        };

        /// \todo Please doc me !
        template <>
        struct hasBackupRestoreFacilities<CpGrid>
        {
            static const bool v = false;
        };

    }


    template<int dim>
    cpgrid::Entity<dim> createEntity(const CpGrid& grid,int index,bool orientation)
    {
        return cpgrid::Entity<dim>(*grid.current_view_data_, index, orientation);
    }

} // namespace Dune

#include <opm/grid/cpgrid/PersistentContainer.hpp>
#include <opm/grid/cpgrid/CartesianIndexMapper.hpp>
#endif // OPM_CPGRID_HEADER
