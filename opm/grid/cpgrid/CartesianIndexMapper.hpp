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
        static const long long dimension = 3 ;
    protected:
        typedef CpGrid Grid;
        const Grid& grid_;
        const long long cartesianSize_;

        long long computeCartesianSize() const
        {
            long long size = cartesianDimensions()[ 0 ];
            for( long long d=1; d<dimension; ++d )
                size *= cartesianDimensions()[ d ];
            return size;
        }

    public:
        explicit CartesianIndexMapper( const Grid& grid )
            : grid_( grid ),
              cartesianSize_( computeCartesianSize() )
        {
        }

        const std::array<long long, dimension>& cartesianDimensions() const
        {
            return grid_.logicalCartesianSize();
        }

        long long cartesianSize() const
        {
            return cartesianSize_;
        }

        long long compressedSize() const
        {
            return grid_.globalCell().size();
        }

        long long cartesianIndex( const long long compressedElementIndex ) const
        {
            assert(  compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
            return grid_.globalCell()[ compressedElementIndex ];
        }

        void cartesianCoordinate(const long long compressedElementIndex, std::array<long long,dimension>& coords) const
        {
            grid_.getIJK( compressedElementIndex, coords );
        }
    };

} // end namespace Opm
#endif
