#ifndef OPM_CPGRIDCARTESIANINDEXMAPPER_HEADER
#define OPM_CPGRIDCARTESIANINDEXMAPPER_HEADER

#include <array>
#include <cassert>

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>

void lookup_check(const Dune::CpGrid&);

namespace Dune
{
template<>
class CartesianIndexMapper< CpGrid >
{
    friend
    void ::lookup_check(const Dune::CpGrid&);
public:
    static const int dimension = 3 ;
protected:
    typedef CpGrid Grid;
    const Grid& grid_;
    const int cartesianSize_;

    int computeCartesianSize() const
    {
        int size = cartesianDimensions()[ 0 ];
        for( int d=1; d<dimension; ++d )
            size *= cartesianDimensions()[ d ];
        return size ;
    }

public:
    explicit CartesianIndexMapper( const Grid& grid )
        : grid_( grid ),
          cartesianSize_( computeCartesianSize() )
    {
    }

    const std::array<int, dimension>& cartesianDimensions() const
    {
        return grid_.logicalCartesianSize();
    }

    int cartesianSize() const
    {
        return cartesianSize_;
    }

    int compressedSize() const
    {
        return  grid_.globalCell().size();
    }

    int cartesianIndex( const int compressedElementIndex ) const
    {
        assert(  compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
        return grid_.globalCell()[compressedElementIndex];
    }

    void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords) const
    {
        // Build entity so we can get its origin
        const auto& entity = Dune::cpgrid::Entity<0>( *grid_.current_view_data_, compressedElementIndex, true);
        // Get origin-cell index (parent cell if existent, or equivalent cell in level 0), and compute its IJK [in level 0]
        (*(grid_.data_[0])).getIJK(entity.getOrigin().index(), coords);
    }

    /* void cartesianCoordinateInLevel(const int compressedElementIndex, std::array<int,dimension>& coords) const
    {
        // Get level and index of the cell in that level
        const auto& level_levelIdx = (*(grid_.data_.back())).leaf_to_level_cells_[compressedElementIndex];
        // Compute IJK in the corresponding LGR.
        (*(grid_.data_[level_levelIdx[0]])).getIJK(level_levelIdx[1], coords);
        }*/
};

} // end namespace Opm
#endif
