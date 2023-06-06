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
        std::vector<int> l0_global_cell( grid_.data_[0] ->size(0), 0); 
        std::iota(l0_global_cell.begin()+1, l0_global_cell.end(), 1); // from entry[1], adds +1 per entry: {0,1,2,3,...}
        (*(grid_.data_[0])).global_cell_ = l0_global_cell;
    }

    const std::array<int, dimension>& cartesianDimensions() const
    {
        return (*(grid_.data_[0])).logical_cartesian_size_;
    }

    int cartesianSize() const
    {
        return cartesianSize_;
    }

    int compressedSize() const
    {
        return  grid_.size(0,0); // (*(grid_.data_[0])).global_cell_.size();
    }

    int cartesianIndex( const int compressedElementIndex ) const
    {
        assert(  compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
        return (*(grid_.data_[0])).global_cell_[compressedElementIndex];
    }

    void cartesianCoordinate(const int compressedElementIndex, std::array<int,dimension>& coords) const
    {
        (*(grid_.data_[0])).getIJK(compressedElementIndex, coords);
    }
};

} // end namespace Opm
#endif
