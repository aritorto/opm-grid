#ifndef OPM_CPGRIDCARTESIANINDEXMAPPER_HEADER
#define OPM_CPGRIDCARTESIANINDEXMAPPER_HEADER

#include <array>
#include <cassert>
#include <stdexcept>

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/CpGrid.hpp>

namespace Dune
{
    template<>
    class CartesianIndexMapper< CpGrid >
    {
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
            return size;
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
            return grid_.globalCell().size();
        }

        int compressedLevelZeroSize() const
        {
            
            return (*grid_.currentData()[0]).size(0);
        }

        int cartesianIndex( const int compressedElementIndex ) const
        {
            assert(  compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
            return grid_.globalCell()[ compressedElementIndex ];
        }

        void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords) const
        {
            grid_.getIJK( compressedElementIndex, coords );
        }

        void cartesianCoordinateLevel(const int compressedElementIndexOnLevel, std::array<int,dimension>& coordsOnLevel, int level) const
        {
            if ((level < 0) || (level > grid_.maxLevel())) {
                throw std::invalid_argument("Invalid level.\n");
            }
            (*grid_.currentData()[level]).getIJK( compressedElementIndexOnLevel, coordsOnLevel);
        }
    };

} // end namespace Opm
#endif
