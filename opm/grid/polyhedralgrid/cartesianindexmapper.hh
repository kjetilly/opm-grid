#ifndef OPM_POLYHEDRALCARTESIANINDEXMAPPER_HEADER
#define OPM_POLYHEDRALCARTESIANINDEXMAPPER_HEADER

#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/polyhedralgrid.hh>

namespace Dune
{
    template< long long dim, long long dimworld, typename coord_t >
    class CartesianIndexMapper< PolyhedralGrid< dim, dimworld, coord_t > >
    {
        typedef PolyhedralGrid< dim, dimworld, coord_t >  Grid;

        const Grid& grid_;
        const long long cartesianSize_;

        long long computeCartesianSize() const
        {
            long long size = cartesianDimensions()[ 0 ];
            for( long long d=1; d<dim; ++d )
                size *= cartesianDimensions()[ d ];
            return size ;
        }
    public:
        static const long long dimension = Grid :: dimension ;

        explicit CartesianIndexMapper( const Grid& grid )
          : grid_( grid ),
            cartesianSize_( computeCartesianSize() )
        {}

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
            return grid_.size( 0 );
        }

        long long cartesianIndex( const long long compressedElementIndex ) const
        {
            assert( compressedElementIndex >= 0 && compressedElementIndex < compressedSize() );
            return grid_.globalCell()[ compressedElementIndex ];
        }

        void cartesianCoordinate(const long long compressedElementIndex, std::array<long long,dimension>& coords) const
        {
          long long gc = cartesianIndex( compressedElementIndex );
          if( dimension >=2 )
          {
              for( long long d=0; d<dimension-2; ++d )
              {
                coords[d] = gc % cartesianDimensions()[d];  gc /= cartesianDimensions()[d];
              }

              coords[dimension-2] = gc % cartesianDimensions()[dimension-2];
              coords[dimension-1] = gc / cartesianDimensions()[dimension-1];
          }
          else
              coords[ 0 ] = gc ;
        }
    };

} // end namespace Opm
#endif
