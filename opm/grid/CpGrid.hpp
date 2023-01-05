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

        // --------------- GET LEAF VIEW FROM 2 LEVELS------SINGLE CELL------
        // Assume we have Level 0 (a CpGridData object) stored in the entry 0
        // of a vector "data" of shared pointers of CpGridData objects.
        // We choose a cell from level 0, and refine it. Store this information
        // in entry 1 of "data". Create the entry 2 of data with the leaf view
        // built with the entities from level 0 that weren't involded in the
        // refinenment together with the new born entities created in level 1.
        //
        // @param data            Vector of shared pointers of CpGridData type.
        //                        data[0] points at the CpGridData where a cell will be chosen, to be refined.
        // @param cells_per_dim   Amount of new born cells per dimension ({in x-direction, y-direction, z-direction})
        // @param parent_idx      Index of the chosen cell from level 0 to be refined.
        void getLeafView2Levels(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                                const std::array<int,3>& cells_per_dim,
                                const int& parent_idx)
        {
            // Build level 1 from the selected cell from level 0 (level 0 = data[0]).
            const auto& [level1_ptr, parent_to_refined_corners,
                         parent_to_children_faces, parent_to_children_cells, child_to_parent_faces, child_to_parent_cell,
                         isParent_faces, isParent_cells]
                = (*data[0]).refineSingleCell(cells_per_dim, parent_idx);
            // Add level 1 to "data".
            data.push_back(level1_ptr);
            // Get some information about the parent cell.
            // Parent cell corners.
            std::vector<int> parent_corners;
            parent_corners.resize(parent_to_refined_corners.size());
            for (int c = 0; c < parent_to_refined_corners.size(); ++c){
                parent_corners[c] = parent_to_refined_corners[c][0];
            }
            // Parent cell faces.
            const auto& parent_faces =  (*data[0]).cell_to_face_[Dune::cpgrid::EntityRep<0>(parent_idx, true)];

            // To store the leaf view.
            typedef Dune::FieldVector<double,3> PointType;
            std::shared_ptr<Dune::cpgrid::CpGridData> leaf_view_ptr =
                std::make_shared<Dune::cpgrid::CpGridData>((*data[0]).ccobj_);
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

            // LEAF CORNER MAP
            // This map will generate a consecutive index-numbering, associating this leaf index with
            // the corresponding level (0 or 1) index.
            int corner_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_corners;
            // Corners coming from the level 0, EXCLUDING parent_corners.
            // Check all the corners from level 0.
            for (int corner = 0; corner < data[0]->size(3); ++corner) {
                // Auxiliary bool to check if a corner from level 0 is a corner of the parent cell
                bool isThere_corn = false;
                // Check if the corner does not belong to the patch. In that case, store it.
                // If the face does not belong to the cell to be refined, store it.
                for(auto& parent_cell_corn : parent_corners) {
                    isThere_corn = isThere_corn || (corner == parent_cell_corn); //true->corn coincides with one parent cell corns
                }
                if(!isThere_corn) {
                    level_to_leaf_corners[{0, corner}] = corner_count;
                    corner_count +=1;
                }
            }
            // Corners coming from level 1, i.e. refined corners.
            // Notice that new born corners start at the entry "total level 0 corners - parent_corners.size()"
            for (int corner = 0; corner < data[1]->size(3); ++corner) {
                level_to_leaf_corners[{1, corner}] = corner_count;
                corner_count +=1;
            }
            // Resize the container of the leaf corners (size: total level 0 - parent_corners.size() + total level 1).
            leaf_corners.resize(corner_count);
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_corners) {
                leaf_corners[leaf_idx] =
                    (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,3>()).get(level_levelIdx[1]);
            }
            // Map to replace old parent cell corners with their equivalent new refined corners.
            // Recall that parent_to_refined_corners = {{parent '0', refined '0'}, {parent 1, refined '1'}, ...}
            std::map<std::array<int,2>, std::array<int,2>> old_to_new_corners;
            for (int corner = 0; corner < parent_to_refined_corners.size(); ++corner) {
                old_to_new_corners[{0,  parent_to_refined_corners[corner][0]}] = {1, parent_to_refined_corners[corner][1]};
            }
            // LEAF FACE MAP
            // This map will generate a consecutive index-numbering, associating this leaf index with
            // the corresponding level (0 or 1) index.
            int face_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_faces;
            // Faces coming from the level 0, that do not belong to the parent cell (that got refined).
            // Check all the faces from level 0.
            for (int face = 0; face < data[0]->face_to_cell_.size(); ++face) {
                // Auxiliary bool to check if a face coincides with a parent cell face.
                bool isThere_face = false;
                // If the face does not belong to the cell to be refined, store it.
                for(auto& parent_cell_face : parent_faces) {
                    isThere_face = isThere_face || (face == parent_cell_face.index()); //true->face coincides with one parent cell face
                }
                if(!isThere_face) {
                    level_to_leaf_faces[{0, face}] = face_count;
                    face_count +=1;
                }
            }
            // Faces coming from level 1, i.e. refined faces.
            // Notice that new born faces start at entry "total level 0 faces - parent_faces.size()"
            for (int face = 0; face < data[1]-> face_to_cell_.size(); ++face) {
                level_to_leaf_faces[{1, face}] = face_count;
                face_count +=1;
            }
            // Resize leaf_faces, mutable_face_tags/normals.
            leaf_faces.resize(face_count);
            mutable_face_tags.resize(face_count);
            mutable_face_normals.resize(face_count);
            // Create the leaf faces, their tags, normals, and corners (face_to_point).
            // Auxiliary vector to store face_to_point "with non consecutive indices". Later on,
            // combine it with the leaf face map, to store the topogolical information with the
            // consecutive index numbering of entities in the leaf view.
            std::vector<std::vector<int>> aux_face_to_point;
            // Auxiliary integer to count the total points in leaf_face_to_point.
            int num_points = 0;
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_faces) {
                // Compute face, face tag, normal, and corners, for faces coming from level 0.
                // Get the face geometry.
                leaf_faces[leaf_idx] = (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,1>())
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the face tag.
                mutable_face_tags[leaf_idx] = (*data[level_levelIdx[0]]).face_tag_
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the face normal.
                mutable_face_normals[leaf_idx] = (*data[level_levelIdx[0]]).face_normals_
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the  leaf indices of corners of the face. Add this info to aux_face_to_point
                // in the right entry to later on populate leaf_face_to_point using appendRow().
                auto old_face_to_point = (*data[level_levelIdx[0]]).face_to_point_[level_levelIdx[1]];
                aux_face_to_point[leaf_idx].reserve(old_face_to_point.size());
                num_points += old_face_to_point.size();
                if (level_levelIdx[0] == 0) {
                    for (int corn = 0; corn < old_face_to_point.size(); ++corn) {
                        // Check if the corner is one of the parent cell corners that got replaced by a refined one.
                        // In that case, we use the refined corner, via the map connecting old parent corners and
                        // new refined ones.
                        // If the corner was not involved in the refinement:
                        bool isThere_corn = false;
                        for(auto& parent_cell_corn : parent_corners) {
                            isThere_corn = isThere_corn || (corn == parent_cell_corn); //true->corn coincides with one parent cell corns
                        }
                        if(!isThere_corn) {
                            aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners[{0, old_face_to_point[corn]}]);
                        }
                        // If the corner was involved in the refinement (one of the 8 corners of the parent cell):
                        else {
                            aux_face_to_point[leaf_idx].push_back(
                                                                  level_to_leaf_corners
                                                                  [old_to_new_corners[{0, old_face_to_point[corn]}]]);
                        }
                    }
                }
                else {
                    for (int corn = 0; corn < old_face_to_point.size(); ++corn) {
                        aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners[{1, old_face_to_point[corn]}]);
                    }
                }
            }
            // LEAF FACE_TO_POINT
            leaf_face_to_point.reserve(face_count, num_points);
            for (int face = 0; face < face_count; ++face) {
                leaf_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
            }
            // Map to connect each of the old parent faces with its children faces.
            std::map<std::array<int,2>,std::vector<std::array<int,2>>> old_to_new_faces;
            for (auto& face : parent_faces) {
                std::vector<int> children = std::get<1>(parent_to_children_faces[face.index()]);
                for (int child = 0; child < children.size(); ++child) {
                    old_to_new_faces[{0,face.index()}].push_back({1, children[child]});
                }
            }
            // CELLS
            int cell_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_cells;
            // Cells coming from the level 0, that are not the parent cell (that got refined).
            // All (cell) indices from level 0 smaller than parent_idx.
            for (int cell = 0; cell < parent_idx; ++cell) {
                level_to_leaf_cells[{0, cell}] = cell_count;
                cell_count +=1;
            }
            // All (cell) indices from level 0 larger than parent_idx.
            for (int cell = parent_idx +1; cell < data[0]-> size(0); ++cell) {
                level_to_leaf_cells[{0, cell}] = cell_count;
                cell_count +=1;
            }
            // Cells coming from level 1, i.e. refined cells.
            // Notice that new born cells start at entry "total level 0 cells - 1"
            for (int cell = 0; cell < data[1]-> size(0); ++cell) {
                level_to_leaf_cells[{1, cell}] = cell_count;
                cell_count +=1;
            }
            leaf_cells.resize(cell_count);
            leaf_cell_to_point.resize(cell_count);
            // As it was done for faces, we store in an additional vector the topological
            // information cell_to_face and later on, we will use a map to store with the
            // consecutive indices the leaf cell_to_face information.
            std::map<int,std::vector<cpgrid::EntityRep<1>>> aux_cell_to_face;
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_cells) {
                leaf_cells[leaf_idx] = (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,0>())
                    [Dune::cpgrid::EntityRep<0>(level_levelIdx[1], true)];
                auto old_cell_to_face = (*data[level_levelIdx[0]]).cell_to_face_[Dune::cpgrid::EntityRep<0>(level_levelIdx[1], true)];
                auto old_cell_to_point = (*data[level_levelIdx[0]]).cell_to_point_[level_levelIdx[1]];
                if (level_levelIdx[0] == 0) {
                    // CELL TO POINT
                    for (int corn = 0; corn < old_cell_to_point.size(); ++corn) {
                        // Auxiliary bool to identify corners of the parent cell.
                        bool isThere_corn = false;
                        for(auto& parent_cell_corn : parent_corners) {
                            isThere_corn = isThere_corn || (corn == parent_cell_corn); //true->corn coincides with one parent cell corns
                        }
                        if(isThere_corn) {
                            leaf_cell_to_point[leaf_idx][corn] =
                                level_to_leaf_corners[old_to_new_corners[{0, old_cell_to_point[corn]}]];
                        }
                        else {
                            leaf_cell_to_point[leaf_idx][corn] = level_to_leaf_corners[{0, old_cell_to_point[corn]}];
                        }
                    }
                    // CELL TO FACE
                    for (auto& face : old_cell_to_face) {
                        // Auxiliary bool to identify faces of the parent cell.
                        bool isThere_face = false;
                        for(auto& parent_cell_face : parent_faces) {
                            isThere_face = isThere_face || (face == parent_cell_face);
                            // true->face coincides with one of the parent cell faces
                        }
                        if(isThere_face) {
                            for (auto& level_newFace : old_to_new_faces[{0, face.index()}]) {
                                aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[level_newFace],
                                        face.orientation()}); // orientation
                            }
                        }
                        // Otherwise, keep the face but with its new leaf index.
                        else {
                            aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[{0, face.index()}], // neigh cell
                                    face.orientation()}); // orientation
                        }
                    }

                }
                else {
                    // CELL TO POINT
                    for (int corn = 0; corn < old_cell_to_point.size(); ++corn) {
                        leaf_cell_to_point[leaf_idx][corn] = level_to_leaf_corners[{1, old_cell_to_point[corn]}];
                    }
                    // CELL TO FACE
                    for (int face = 0; face < old_cell_to_face.size(); ++face) {
                        aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[{1, old_cell_to_face[face].index()}], // neigh cell
                                old_cell_to_face[face].orientation()}); // orientation
                    }
                }
            }
            // LEAF CELL TO FACE
            for (int cell = 0; cell < cell_count; ++cell) {
                leaf_cell_to_face.appendRow(aux_cell_to_face[cell].begin(), aux_cell_to_face[cell].end());
            }
            // LEAF FACE TO CELL
            leaf_cell_to_face.makeInverseRelation(leaf_face_to_cell);
            //  Add level 2 (leafview) to "data".
            data.push_back(leaf_view_ptr);
        }


        // --------------- GET LEAF VIEW FROM 2 LEVELS------BLOCK PATCH------
        // Assume we have Level 0 (a CpGridData object) stored in the entry 0
        // of a vector "data" of type shared pointers of CpGridData objects.
        // We choose a (block) patch from level 0, and refine it. Store this information
        // in entry 1 of "data". Create the entry 2 of data with the leaf view
        // construct with the entities from level 0 that weren't involded in the
        // refinenment together with the new born entities created in level 1.
        // @param data            Vector of shared pointers of CpGridData type.
        //                        data[0] points at the CpGridData where a patch of cells will be chosen, to be refined.
        // @param cells_per_dim   Amount of new born cells per dimension ({in x-direction, y-direction, z-direction})
        // @param start/end_ijk   To build the patch to be refined.
        void getLeafView2LevelsPatch(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                                     const std::array<int,3>& cells_per_dim,
                                     const std::array<int,3>& start_ijk,
                                     const std::array<int,3>& end_ijk)
        {
            // Build level 1 from the selected patch from level 0 (level 0 = data[0]).
            const auto& [level1_ptr, boundary_old_to_new_corners, boundary_old_to_new_faces,
                         parent_to_children_faces, parent_to_children_cells,
                         child_to_parent_faces, child_to_parent_cells,
                         isParent_faces, isParent_cells]
                = (*data[0]).refineBlockPatch(cells_per_dim, start_ijk, end_ijk);
            // "level1_ptr" is a CpGridData object. New born corners, faces, cells, with
            //              their topological information, e.g., cell_to_face_.
            // Add level 1 to "data".
            data.push_back(level1_ptr);
            // To store the leaf view.
            typedef Dune::FieldVector<double,3> PointType;
            std::shared_ptr<Dune::cpgrid::CpGridData> leaf_view_ptr =
                std::make_shared<Dune::cpgrid::CpGridData>((*data[0]).ccobj_);
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

            // Get patch corners, faces, and cells indices
            const auto& [patch_corners, patch_faces, patch_cells] = (*data[0]).getPatchGeomIndices(start_ijk, end_ijk);

            // LEAF CORNER MAP
            int corner_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_corners;
            // Corners coming from the level 0, EXCLUDING patch_corners.
            for (int corner = 0; corner < data[0]->size(3); ++corner) {
                // Auxiliary bool to discard patch corners.
                bool isThere_corn = false;
                for(auto& patch_corn : patch_corners) {
                    isThere_corn = isThere_corn || (corner == patch_corn); //true->corn coincides with one patch_corners
                }
                if(!isThere_corn) {
                    level_to_leaf_corners[{0, corner}] = corner_count;
                    corner_count +=1;
                }
            }
            // Corners coming from level 1, i.e. refined corners.
            // Notice that new born corners start at entry "total level 0 corners - patch_corners.size()"
            for (int corner = 0; corner < (data[1]->size(3)); ++corner) {
                level_to_leaf_corners[{1, corner}] = corner_count;
                corner_count +=1;
            }
            // Resize the container of the leaf corners (size: total level 0 - patch_corners.size() + total level 1).
            leaf_corners.resize(corner_count);
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_corners) {
                leaf_corners[leaf_idx] =
                    (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,3>()).get(level_levelIdx[1]);
            }
            // Assume each parent cell has 8 corners (maybe allowing repetition of corners?).
            // parent_cells_to_new_corners[fixed cell index] = {{parent corner '0', equivalent refined corner '0'},
            //                                                  {parent corner '1', equivalent refined corner '1'}, ...}
            // For all parent cells, we connect old corners to new ones, and store them in a map:
            std::map<std::array<int,2>, std::array<int,2>> old_to_new_corners;
            // Auxiliary vector with boundary corners of the patch
            std::vector<int> boundary_patch_corners;
            boundary_patch_corners.reserve(boundary_old_to_new_corners.size());
            for (int corner = 0; corner < boundary_old_to_new_corners.size(); ++corner) {
                old_to_new_corners[{0, boundary_old_to_new_corners[corner][0]}] = {1, boundary_old_to_new_corners[corner][1]};
                boundary_patch_corners.push_back(boundary_old_to_new_corners[corner][0]);
            }
            // FACES
            int face_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_faces;
            // Faces coming from the level 0, that do not belong to the patch (that got refined).
            for (int face = 0; face < data[0]->face_to_cell_.size(); ++face) {
                // Auxiliary bool to discard patch faces.
                bool isThere_face = false;
                for(auto& patch_face : patch_faces) {
                    isThere_face = isThere_face || (face == patch_face); //true->face coincides with one patch faces
                }
                if(!isThere_face) {
                    //  if (std::find(patch_faces.begin(), patch_faces.end(), face) ==  patch_faces.end()) {
                    level_to_leaf_faces[{0, face}] = face_count;
                    face_count +=1;
                }
            }
            // Faces coming from level 1, i.e. refined faces.
            // Notice that new born faces start at entry "total level 0 faces - patch_faces.size()"
            // int start_new_faces = (data[0]-> face_to_cell_.size()) - patch_faces.size();
            for (int face = 0; face < (data[1]-> face_to_cell_.size()); ++face) {
                level_to_leaf_faces[{1, face}] = face_count;
                face_count +=1;
            }
            // Resize leaf_faces, mutable_face_tags/normals.
            leaf_faces.resize(face_count);
            mutable_face_tags.resize(face_count);
            mutable_face_normals.resize(face_count);
            // Auxiliary integer to count all the points in leaf_face_to_point.
            int num_points = 0;
            // Create the leaf faces, their tags, normals, and 4 corners.
            // Auxiliary vector to store face_to_point "with non consecutive indices". Later on,
            // we use this and a map, to be able to store the topogolical information with the
            // consecutive index numbering of entities in the lead view.
            std::vector<std::vector<int>> aux_face_to_point;
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_faces) {
                // Compute face, face tag, normal, and corners, for faces coming from level 0.
                // Get the face geometry.
                leaf_faces[leaf_idx] = (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,1>())
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the face tag.
                mutable_face_tags[leaf_idx] = (*data[level_levelIdx[0]]).face_tag_
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the face normal.
                mutable_face_normals[leaf_idx] = (*data[level_levelIdx[0]]).face_normals_
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the  leaf indices of corners of the face. Add this info to aux_face_to_point
                // in the right entry to later on populate leaf_face_to_point using appendRow().
                auto old_face_to_point = (*data[level_levelIdx[0]]).face_to_point_[level_levelIdx[1]];
                aux_face_to_point[leaf_idx].reserve(old_face_to_point.size()); // we can put 4 here
                // Add the amount of points to the count num_points.
                num_points += old_face_to_point.size();
                if (level_levelIdx[0] == 0) {
                    for (int corn = 0; corn < old_face_to_point.size(); ++corn) { // we can put 4 here
                        // Check if the corner is one of the parent cell corners that got replaced by a refined one.
                        // In that case, we use the refined corner, via the map connecting old parent corners and
                        // new refined ones.
                        // If the corner was not involved in the refinement:
                        // Auxiliary bool to identify boundary patch corners
                        bool isThere_bound_corn = false;
                        for(auto& bound_corn : boundary_patch_corners) {
                            isThere_bound_corn = isThere_bound_corn || (corn == bound_corn); //true->coincides with one boundary corn
                        }
                        // If it does not belong to the boundary of the patch:
                        if(!isThere_bound_corn) {
                            aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners[{0, old_face_to_point[corn]}]);
                        }
                        // If the corner was involved in the refinement (one of the corners of a cell on the boundary of the patch):
                        else {
                            aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners
                                                                  [old_to_new_corners[{0, old_face_to_point[corn]}]]);
                        }
                    }
                }
                else {
                    for (int corn = 0; corn < old_face_to_point.size(); ++corn) {
                        aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners[{1, old_face_to_point[corn]}]);
                    }
                }
            }
            // LEAF FACE_TO_POINT
            leaf_face_to_point.reserve(face_count, num_points);
            for (int face = 0; face < face_count; ++face) {
                leaf_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
            }

            // Map to connect each of the old parent faces with its children faces.
            // Recall const std::vector<std::tuple<int, std::vector<int>>> boundary_old_to_new_faces
            std::map<std::array<int,2>,std::vector<std::array<int,2>>> old_to_new_faces;
            std::vector<int> boundary_patch_faces;
            boundary_patch_faces.reserve(boundary_old_to_new_faces.size());
            for (int face = 0; face < boundary_old_to_new_faces.size(); ++face) {
                for (auto& child : std::get<1>(boundary_old_to_new_faces[face])) {
                    old_to_new_faces[{0, std::get<0>(boundary_old_to_new_faces[face])}].push_back({1, child});
                }
                boundary_patch_faces.push_back(std::get<0>(boundary_old_to_new_faces[face]));
            }

            // CELLS
            int cell_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_cells;
            // Cells coming from the level 0, that are not the patch cells (that got refined).
            for (int cell = 0; cell < data[0]-> size(0); ++cell) {
                // Auxiliary bool to identify cells of the patch.
                bool isThere_cell = false;
                for(auto& patch_cell : patch_cells) {
                    isThere_cell = isThere_cell || (cell == patch_cell); //true-> coincides with one patch cells
                }
                // If it does not belong to the patch:
                if(!isThere_cell) {
                    level_to_leaf_cells[{0, cell}] = cell_count;
                    cell_count +=1;
                }
            }
            // Cells coming from level 1, i.e. refined cells.
            // Notice that new born cells start at entry "total level 0 cells - patch_cells.size()"
            for (int cell = 0; cell < data[1]-> size(0); ++cell) {
                level_to_leaf_cells[{1, cell}] = cell_count;
                cell_count +=1;
            }
            leaf_cells.resize(cell_count);
            leaf_cell_to_point.resize(cell_count);
            // As it was done for faces, we store in an additional vector the topological
            // information cell_to_face and later on, we will use a map to store with the
            // consecutive indices the leaf cell_to_face information.
            std::map<int,std::vector<cpgrid::EntityRep<1>>> aux_cell_to_face;
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_cells) {
                // Get the cell geometry.
                leaf_cells[leaf_idx] = (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,0>())
                    [Dune::cpgrid::EntityRep<0>(level_levelIdx[1], true)];
                // Get old faces of the cell that will be replaced with leaf faces
                auto old_cell_to_face = (*data[level_levelIdx[0]]).cell_to_face_[Dune::cpgrid::EntityRep<0>(level_levelIdx[1], true)];
                auto old_cell_to_point = (*data[level_levelIdx[0]]).cell_to_point_[level_levelIdx[1]];
                if (level_levelIdx[0] == 0) {
                    // CELL TO POINT
                    // Get the  leaf indices of corners of the cell.
                    for (int corn = 0; corn < old_cell_to_point.size(); ++corn) {
                        // Auxiliary bool to identity patch corners
                        bool isThere_corn = false;
                        for(auto& patch_corn : patch_corners) {
                            isThere_corn = isThere_corn || (corn == patch_corn); //true-> coincides with one patch corns
                        }
                        // If it DOES belong to the patch:
                        if(isThere_corn) {
                            leaf_cell_to_point[leaf_idx][corn] =
                                level_to_leaf_corners[old_to_new_corners[{0, old_cell_to_point[corn]}]];
                        }
                        else {
                            leaf_cell_to_point[leaf_idx][corn] = level_to_leaf_corners[{0, old_cell_to_point[corn]}];
                        }
                    }
                    // CELL TO FACE
                    for (auto& face : old_cell_to_face)
                    {
                        // Auxiliary bool to identity boundary patch faces
                        bool isThere_face = false;
                        // Check if the face is one of those faces that was involved in the refinement.
                        // If it was, replace it with its children faces.
                        for(auto& bound_face : boundary_patch_faces) {
                            isThere_face = isThere_face || (face.index() == bound_face); //true-> coincides with one bound faces
                        }
                        // If it DOES belong to the boundary:
                        if(isThere_face) {
                            for (auto& level_newFace : old_to_new_faces[{0, face.index()}]) {
                                aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[level_newFace], // (new) refined face
                                        face.orientation()}); // orientation
                            }
                        }
                        // Otherwise, keep the face but with its new leaf index.
                        else {
                            aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[{0, face.index()}], // (new) refined face
                                    face.orientation()}); // orientation
                        }
                    }
                }
                else {
                    // CELL TO POINT
                    for (int corn = 0; corn < old_cell_to_point.size(); ++corn) {
                        // Get the  leaf indices of corners of the cell.
                        leaf_cell_to_point[leaf_idx][corn] = level_to_leaf_corners[{1, old_cell_to_point[corn]}];
                    }
                    // CELL TO FACE
                    for (auto& face : old_cell_to_face) {
                        // Get the leaf index and orientation of faces of the cell.
                        aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[{1, face.index()}], // neigh cell
                                face.orientation()}); // orientation
                    }
                }
            }
            // LEAF CELL TO FACE
            for (int cell = 0; cell < cell_count; ++cell) {
                leaf_cell_to_face.appendRow(aux_cell_to_face[cell].begin(), aux_cell_to_face[cell].end());
            }
            // LEAF FACE TO CELL
            leaf_cell_to_face.makeInverseRelation(leaf_face_to_cell);
            //  Add level 2 (leafview) to "data".
            data.push_back(leaf_view_ptr);
        }

        //////////////////////////////////////
        // --------------- GET LEAF VIEW FROM 2 LEVELS----- SET OF CELLS------
        // Assume we have Level 0 (a CpGridData object) stored in the entry 0
        // of a vector "data" of shared pointers of CpGridData objects.
        // We choose some cells from level 0 (not necessary a block of cells), and refine them.
        // Store this information in entry 1 of "data". Create the entry 2 of data with the leaf view
        // built with the entities from level 0 that weren't involded in the
        // refinenment together with the new born entities created in level 1.
        //
        // @param data            Vector of shared pointers of CpGridData type.
        //                        data[0] points at the CpGridData where a cell will be chosen, to be refined.
        // @param cells_per_dim   Amount of new born cells per dimension ({in x-direction, y-direction, z-direction})
        // @param cells2refine    Indices of the chosen cells from level 0 to be refined.
        void getLeafView2LevelsSet(std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>>& data,
                                const std::array<int,3>& cells_per_dim,
                                std::vector<int> cells2refine) 
        {
            // Sort the vector of cell indices to be refined.
            //std::vector<int> sorted_cells2refine =
                std::sort(cells2refine.begin(), cells2refine.end());
            // Amount of cells to be refined.
            int num_cells2refine = cells2refine.size();
            // Container to collect auxiliary CpGridData objects for each parent cell.
            // All CpGridData objects, one per parent cell
            std::vector<std::shared_ptr<Dune::cpgrid::CpGridData>> aux_grids;
            aux_grids.reserve(num_cells2refine);
            // All parent_to_refined_corners (the corners of a parent cell get "replaced by" the new born (coinciding) ones.)
            std::vector<std::vector<std::array<int,2>>> all_parent_to_refined_corners;
            all_parent_to_refined_corners.reserve(num_cells2refine);
            // All parent_to_children_faces
            std::vector<std::vector<std::tuple<int,std::vector<int>>>> all_parent_to_children_faces;
            // all_parent_to_children_faces.reserve(???);
            // All parent_children_cells = {parent cell idx, {child 0, child 1, ...}}
            std::vector<std::tuple<int, std::vector<int>>> all_parent_to_children_cells;
            all_parent_to_children_cells.reserve(num_cells2refine);
            // All child_to_parent_faces [POSTPONED]
            // std::vector<std::vector<std::array<int,2>>> all_child_to_parent_faces;
            // all_child_to_parent_faces.reserve(???)
            // All child_to_parent_cell [POSTPONED]
            // std::vector<std::vector<std::array<int,2>>> all)child_to_parent_cell;
            // All isParent_faces. Default: false. It'll be rewritten for actual parents.
            std::map<int,bool> isParent_faces;
            for (int face = 0; face < (data[0]-> face_to_cell_.size()); ++face) {
                isParent_faces[face] = false;
            }
            // All isParent_cells. Default: false. It'll be rewritten for actual parents.
            std::map<int,bool> isParent_cells;
            for (int c = 0; c < (data[0]-> size(0)); ++c) {
            isParent_cells[c] = false;
            }
            // All REFINED corners, with repetition. We will make a set out of it to 'remove' repeated corners.
            std::vector<Dune::cpgrid::Geometry<0,3>> all_corners_repetition;
            // Map to relate a corner {x,y,z} with its (if any) parent cell, index in its level.
            // The value will be {int, std::vector<int>} where the first int takes value -1 when the corner
            // does not come from refinenment, or parent cell index (!= -1) when the corner comes from
            // refinement. The second integer is the corner index (when first entry is -1, it corresponds
            // to the corner index from level 0, when first entry is > -1, it corresponds to the corner index
            // from the LRG of the parent cell).
            // MODIFY DOCUMENTATION, ONLY REFINED CORNERS HERE
            std::map<Dune::cpgrid::Geometry<0,3>, std::array<int,2>> globalCoord_to_cornId;
            // Reverse map
            std::map<std::array<int,2>,Dune::cpgrid::Geometry<0,3>> cornId_to_globalCoord;
            
            
            // Build one LRG for each cell that get refined. We'll store only one final CpGridData object as a new entry in "data".
            for (auto& cell : cells2refine) {
            // Build level from the selected cell from level 0 (level 0 = data[0]).
            const auto& [level_ptr, parent_to_refined_corners,
                         parent_to_children_faces, parent_to_children_cells, child_to_parent_faces, child_to_parent_cell,
                         level_isParent_faces, level_isParent_cells]
                = (*data[0]).refineSingleCell(cells_per_dim, cell);
            // Add the level to the auxiliary vector with grids. 
            aux_grids.push_back(level_ptr);
            // Add parent_to_refined_corners.
            all_parent_to_refined_corners.push_back(parent_to_refined_corners);
            // Add parent_to_children_faces.
            all_parent_to_children_faces.push_back(parent_to_children_faces);
            // Add parent_to_children_cells.
            all_parent_to_children_cells.push_back(parent_to_children_cells);
            // Add child_to_parent_faces [POSTPONED]
            // Add child_to_parent_cell [POSTPONED]
            // Rewrite the keys of isParent_faces to incorporate the new parents.
            // Recall parent_to_children_faces contain the parent faces and it's type
            // std::vector<std::tuple<int,std::vector<int>>. (tuple: {face idx, {child 0, child 1, ...}})
            for (auto& parent_face : parent_to_children_faces) {
                isParent_faces[std::get<0>(parent_face)] = true;
            }
            // Rewrite the keys of isParent_cells to incorporate the new parents.
            isParent_cells[cell] = true;
               auto level_corners =  (*level_ptr).geometry_.geomVector(std::integral_constant<int,3>());
               
                 int corn_count = 0;
             for (auto& corn : level_corners){
                auto corn_center = level_corners.get(corn_count).center();
                  // Add refined corners
                all_corners_repetition.push_back(corn_center);
                // Create entries on the maps to relate refined corners with their actual global coordiantes.
                // Notice that maps re-write a value when the key is the same, so we would not repeat them.
                //    globalCoord_to_cornId[ corn_center ] = {cell, corn_count};
                cornId_to_globalCoord[{cell, corn_count}] = corn_center;
                // globalCoord_to_cornId[corn_center] = {cell, corn_count};
                corn_count += 1;
            }
             // All refined corners WITHOUT repetition
            // std::set<Dune::cpgrid::Geometry<0,3>> refined_corners_set;
            std::vector<Dune::cpgrid::Geometry<0,3>> all_different_corners;
            all_different_corners.push_back(all_corners_repetition[0].center());
            for (auto& corn : all_corners_repetition) {
                for (auto& one_diff_corn : all_different_corners) {
                    if ( corn.center() != one_diff_corn.center()) {
                        all_different_corners.push_back(corn);
                    }   
                }
            }
           
            
            
         
            
            /* // Compute boundary/inner refined corners, per parent cell
             const auto& [b_front, b_back, b_left, b_right, b_bottom, b_top, inner]
                    = (*data[0]).geometry_.geomVector(std::integral_constant<int,0>())
                    [Dune::cpgrid::EntityRep<0>(cell, true)].getBoundaryInnerRefinedCorners(cells_per_dim);
             // Add boundary/inner refined corners in their containers.
             all_boundary_front_refined_cornes.push_back(b_front);
             all_boundary_back_refined_cornes.push_back(b_back);
             all_boundary_left_refined_cornes.push_back(b_left);
             all_boundary_right_refined_cornes.push_back(b_right);
             all_boundary_bottom_refined_cornes.push_back(b_bottom);
             all_boundary_top_refined_cornes.push_back(b_top);
             all_boundary_inner_refined_cornes.push_back(inner); */
             
            }
            
            // Two different cells can be completely disjoint, or can intersect in only one corner (they do not
            // share any edge or face), more than one corner. In the case they share more than one corner, they may share
            // an edge-and-no-face or, at least, a face.

            // Get some information about the parent cells.
            // Parent corners. May repeat indices. 
            std::vector<int> parent_corners;
            // Parent faces. May repeat indices. 
            // std::vector<int> parent_faces;
            // To store the refined corners of the set of cells that got refined, without repetition.
            // CONTAINER 
            for (int cell = 0; cell < cells2refine.size(); ++cell) {
            // Parent cell corners.
                //     std::array<int,8> cell2point  = (*data[0]).cell_to_point_[cell];
                //for (int corn = 0; corn < 8; ++corn) {
                //  parent_corners.push_back(cell2point[corn]);
                //   }
                // Boundary refined corners.
                // Boundary front,back,left,right,bottom,top
                // Inner refined corners.
                const auto& [b_front, b_back, b_left, b_right, b_bottom, b_top, inner]
                    = (*data[0]).geometry_.geomVector(std::integral_constant<int,0>())
                    [Dune::cpgrid::EntityRep<0>(cell, true)].getBoundaryInnerRefinedCorners(cells_per_dim);
                // Store the inner corners (they are unique no matter possible neighboring cells)
                // STORE THEM...
                /*     for (int other_cell = 0; other_cell < cells2refine.size(); ++other_cell) {
                    if (cell != other_cell)
                        // for loop in cell faces
                        // check if one of the other_cell faces coincides 
                        }   */
            }
            
            // Parent cell faces.
                //     cpgrid::OrientedEntityTable<0,1> cell2face =  (*data[0]).cell_to_face_[Dune::cpgrid::EntityRep<0>(cell, true)];
                //for (auto& face : cell2face) {
                //  parent_faces.push_back(face.index());
                // }
            //}
            // Sort parent_corners, may contain repeated indices.
            // std::vector<int> sorted_parent_corners =
                 std::sort(parent_corners.begin(), parent_corners.end());
            // Sort parent_faces, may contain repeated indices.
            //std::vector<int> sorted_parent_faces =
                //        std::sort(parent_faces.begin(), parent_faces.end());
            // Remove repeated corners and faces from sorted_parent_corners/faces.
                //     parent_corners.erase(std::unique( parent_corners.begin(), parent_corners.end() ),
            //                            parent_corners.end());
            //   parent_faces.erase(std::unique( parent_faces.begin(), parent_faces.end() ),
            //                             parent_faces.end());
            
            /*   std::map<int,std::vector<std::array<int,8>>> parents_cell2point;
            // Parent faces. face type: cpgrid::OrientedEntityTable<0,1>
            std::map<int,std::vect> parents_cell2face;
            for (auto& cell : sorted_cells2refine){
            // Parent cell corners.
            parents_cell2point[cell] =  (*data[0]).cell_to_point_[cell];
            // Parent cell faces.
            parents_cell2face[cell] =  (*data[0]).cell_to_face_[Dune::cpgrid::EntityRep<0>(cell, true)];
            }*/
            

            // To store the leaf view.
            typedef Dune::FieldVector<double,3> PointType;
            std::shared_ptr<Dune::cpgrid::CpGridData> leaf_view_ptr =
                std::make_shared<Dune::cpgrid::CpGridData>((*data[0]).ccobj_);
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

            // LEAF CORNER MAP
            // This map will generate a consecutive index-numbering, associating this leaf index with
            // the corresponding level (0 or 1) index.
            int corner_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_corners;
            // Corners coming from the level 0, EXCLUDING parents_corners.
            // Check all the corners from level 0.
            for (int corner = 0; corner < data[0]->size(3); ++corner) {
                // Auxiliary bool to check if a corner from level 0 is a corner of a parent cell
                bool isThere_corn = false;
                // Check if the corner does not belong to the patch. In that case, store it.
                for(auto& parent_corn : parent_corners) { // parent_corners does not repeat indices.
                    isThere_corn = isThere_corn || (corner == parent_corn);
                    //true-> corn coincides with one parent cell corns
                    //false-> corn does not belong to any cell that got refined.
                }
                if(!isThere_corn) { // corner was not involved in refinement.
                    level_to_leaf_corners[{0, corner}] = corner_count;
                    corner_count +=1;
                }
                else { // corner was involved in refinement.
                    
                }
                
            }
            // Corners coming from level l, i.e. refined corners. We have one level per refined-cell.
            // Check if the corner belongs 
            // Notice that new born corners start at the entry "total level 0 corners - parent_corners.size()"
            for (int l = 1; l < num_cells2refine +1; ++l){ // levels are stored in data[1], ..., data[num_cells2refine].
            for (int corner = 0; corner < data[l]->size(3); ++corner) {
                level_to_leaf_corners[{l, corner}] = corner_count;
                corner_count +=1;
            }
            }
            
            // Resize the container of the leaf corners (size: total level 0 - parent_corners.size() + total level 1).
            leaf_corners.resize(corner_count);
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_corners) {
                leaf_corners[leaf_idx] =
                    (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,3>()).get(level_levelIdx[1]);
            }
            // Map to replace old parent cell corners with their equivalent new refined corners.
            // Recall that parent_to_refined_corners = {{parent '0', refined '0'}, {parent 1, refined '1'}, ...}
            //   std::map<std::array<int,2>, std::array<int,2>> old_to_new_corners;
            //for (int corner = 0; corner < parent_to_refined_corners.size(); ++corner) {
            //  old_to_new_corners[{0,  parent_to_refined_corners[corner][0]}] = {1, parent_to_refined_corners[corner][1]};
            // }
            /* // LEAF FACE MAP
            // This map will generate a consecutive index-numbering, associating this leaf index with
            // the corresponding level (0 or 1) index.
            int face_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_faces;
            // Faces coming from the level 0, that do not belong to the parent cell (that got refined).
            // Check all the faces from level 0.
            for (int face = 0; face < data[0]->face_to_cell_.size(); ++face) {
                // Auxiliary bool to check if a face coincides with a parent cell face.
                bool isThere_face = false;
                // If the face does not belong to the cell to be refined, store it.
                for(auto& parent_cell_face : parent_faces) {
                    isThere_face = isThere_face || (face == parent_cell_face.index()); //true->face coincides with one parent cell face
                }
                if(!isThere_face) {
                    level_to_leaf_faces[{0, face}] = face_count;
                    face_count +=1;
                }
            }
            // Faces coming from level 1, i.e. refined faces.
            // Notice that new born faces start at entry "total level 0 faces - parent_faces.size()"
            for (int face = 0; face < data[1]-> face_to_cell_.size(); ++face) {
                level_to_leaf_faces[{1, face}] = face_count;
                face_count +=1;
            }
            // Resize leaf_faces, mutable_face_tags/normals.
            leaf_faces.resize(face_count);
            mutable_face_tags.resize(face_count);
            mutable_face_normals.resize(face_count);
            // Create the leaf faces, their tags, normals, and corners (face_to_point).
            // Auxiliary vector to store face_to_point "with non consecutive indices". Later on,
            // combine it with the leaf face map, to store the topogolical information with the
            // consecutive index numbering of entities in the leaf view.
            std::vector<std::vector<int>> aux_face_to_point;
            // Auxiliary integer to count the total points in leaf_face_to_point.
            int num_points = 0;
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_faces) {
                // Compute face, face tag, normal, and corners, for faces coming from level 0.
                // Get the face geometry.
                leaf_faces[leaf_idx] = (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,1>())
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the face tag.
                mutable_face_tags[leaf_idx] = (*data[level_levelIdx[0]]).face_tag_
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the face normal.
                mutable_face_normals[leaf_idx] = (*data[level_levelIdx[0]]).face_normals_
                    [Dune::cpgrid::EntityRep<1>(level_levelIdx[1], true)];
                // Get the  leaf indices of corners of the face. Add this info to aux_face_to_point
                // in the right entry to later on populate leaf_face_to_point using appendRow().
                auto old_face_to_point = (*data[level_levelIdx[0]]).face_to_point_[level_levelIdx[1]];
                aux_face_to_point[leaf_idx].reserve(old_face_to_point.size());
                num_points += old_face_to_point.size();
                if (level_levelIdx[0] == 0) {
                    for (int corn = 0; corn < old_face_to_point.size(); ++corn) {
                        // Check if the corner is one of the parent cell corners that got replaced by a refined one.
                        // In that case, we use the refined corner, via the map connecting old parent corners and
                        // new refined ones.
                        // If the corner was not involved in the refinement:
                        bool isThere_corn = false;
                        for(auto& parent_cell_corn : parent_corners) {
                            isThere_corn = isThere_corn || (corn == parent_cell_corn); //true->corn coincides with one parent cell corns
                        }
                        if(!isThere_corn) {
                            aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners[{0, old_face_to_point[corn]}]);
                        }
                        // If the corner was involved in the refinement (one of the 8 corners of the parent cell):
                        else {
                            aux_face_to_point[leaf_idx].push_back(
                                                                  level_to_leaf_corners
                                                                  [old_to_new_corners[{0, old_face_to_point[corn]}]]);
                        }
                    }
                }
                else {
                    for (int corn = 0; corn < old_face_to_point.size(); ++corn) {
                        aux_face_to_point[leaf_idx].push_back(level_to_leaf_corners[{1, old_face_to_point[corn]}]);
                    }
                }
            }
            // LEAF FACE_TO_POINT
            leaf_face_to_point.reserve(face_count, num_points);
            for (int face = 0; face < face_count; ++face) {
                leaf_face_to_point.appendRow(aux_face_to_point[face].begin(), aux_face_to_point[face].end());
            }
            // Map to connect each of the old parent faces with its children faces.
            std::map<std::array<int,2>,std::vector<std::array<int,2>>> old_to_new_faces;
            for (auto& face : parent_faces) {
                std::vector<int> children = std::get<1>(parent_to_children_faces[face.index()]);
                for (int child = 0; child < children.size(); ++child) {
                    old_to_new_faces[{0,face.index()}].push_back({1, children[child]});
                }
            }
            // CELLS
            int cell_count = 0;
            std::map<std::array<int,2>, int> level_to_leaf_cells;
            // Cells coming from the level 0, that are not the parent cell (that got refined).
            // All (cell) indices from level 0 smaller than parent_idx.
            for (int cell = 0; cell < parent_idx; ++cell) {
                level_to_leaf_cells[{0, cell}] = cell_count;
                cell_count +=1;
            }
            // All (cell) indices from level 0 larger than parent_idx.
            for (int cell = parent_idx +1; cell < data[0]-> size(0); ++cell) {
                level_to_leaf_cells[{0, cell}] = cell_count;
                cell_count +=1;
            }
            // Cells coming from level 1, i.e. refined cells.
            // Notice that new born cells start at entry "total level 0 cells - 1"
            for (int cell = 0; cell < data[1]-> size(0); ++cell) {
                level_to_leaf_cells[{1, cell}] = cell_count;
                cell_count +=1;
            }
            leaf_cells.resize(cell_count);
            leaf_cell_to_point.resize(cell_count);
            // As it was done for faces, we store in an additional vector the topological
            // information cell_to_face and later on, we will use a map to store with the
            // consecutive indices the leaf cell_to_face information.
            std::map<int,std::vector<cpgrid::EntityRep<1>>> aux_cell_to_face;
            for (auto& [level_levelIdx, leaf_idx] : level_to_leaf_cells) {
                leaf_cells[leaf_idx] = (*data[level_levelIdx[0]]).geometry_.geomVector(std::integral_constant<int,0>())
                    [Dune::cpgrid::EntityRep<0>(level_levelIdx[1], true)];
                auto old_cell_to_face = (*data[level_levelIdx[0]]).cell_to_face_[Dune::cpgrid::EntityRep<0>(level_levelIdx[1], true)];
                auto old_cell_to_point = (*data[level_levelIdx[0]]).cell_to_point_[level_levelIdx[1]];
                if (level_levelIdx[0] == 0) {
                    // CELL TO POINT
                    for (int corn = 0; corn < old_cell_to_point.size(); ++corn) {
                        // Auxiliary bool to identify corners of the parent cell.
                        bool isThere_corn = false;
                        for(auto& parent_cell_corn : parent_corners) {
                            isThere_corn = isThere_corn || (corn == parent_cell_corn); //true->corn coincides with one parent cell corns
                        }
                        if(isThere_corn) {
                            leaf_cell_to_point[leaf_idx][corn] =
                                level_to_leaf_corners[old_to_new_corners[{0, old_cell_to_point[corn]}]];
                        }
                        else {
                            leaf_cell_to_point[leaf_idx][corn] = level_to_leaf_corners[{0, old_cell_to_point[corn]}];
                        }
                    }
                    // CELL TO FACE
                    for (auto& face : old_cell_to_face) {
                        // Auxiliary bool to identify faces of the parent cell.
                        bool isThere_face = false;
                        for(auto& parent_cell_face : parent_faces) {
                            isThere_face = isThere_face || (face == parent_cell_face);
                            // true->face coincides with one of the parent cell faces
                        }
                        if(isThere_face) {
                            for (auto& level_newFace : old_to_new_faces[{0, face.index()}]) {
                                aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[level_newFace],
                                        face.orientation()}); // orientation
                            }
                        }
                        // Otherwise, keep the face but with its new leaf index.
                        else {
                            aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[{0, face.index()}], // neigh cell
                                    face.orientation()}); // orientation
                        }
                    }

                }
                else {
                    // CELL TO POINT
                    for (int corn = 0; corn < old_cell_to_point.size(); ++corn) {
                        leaf_cell_to_point[leaf_idx][corn] = level_to_leaf_corners[{1, old_cell_to_point[corn]}];
                    }
                    // CELL TO FACE
                    for (int face = 0; face < old_cell_to_face.size(); ++face) {
                        aux_cell_to_face[leaf_idx].push_back({level_to_leaf_faces[{1, old_cell_to_face[face].index()}], // neigh cell
                                old_cell_to_face[face].orientation()}); // orientation
                    }
                }
            }
            // LEAF CELL TO FACE
            for (int cell = 0; cell < cell_count; ++cell) {
                leaf_cell_to_face.appendRow(aux_cell_to_face[cell].begin(), aux_cell_to_face[cell].end());
            }
            // LEAF FACE TO CELL
            leaf_cell_to_face.makeInverseRelation(leaf_face_to_cell);
            //  Add level 2 (leafview) to "data".
            */
            data.push_back(leaf_view_ptr);
        }
        //////////////////////////////////////
        
        
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
