#ifndef OPM_CARTESIANINDEXMAPPER_HEADER
#define OPM_CARTESIANINDEXMAPPER_HEADER

#include <array>
#include <cassert>

#include <dune/common/exceptions.hh>

namespace Dune
{
    /** \brief Interface class to access the logical Cartesian grid as used in industry
               standard simulator decks.
               */
    template< class Grid >
    class CartesianIndexMapper
    {
    public:
        /** \brief dimension of the grid */
        static const long long dimension = Grid :: dimension ;

        /** \brief constructor taking grid */
        explicit CartesianIndexMapper( const Grid& )
        {
            DUNE_THROW(InvalidStateException,"CartesianIndexMapper not specialized for given grid");
        }

        /** \brief return Cartesian dimensions, i.e. number of cells in each direction  */
        const std::array<long long, dimension>& cartesianDimensions() const
        {
            static std::array<long long, dimension> a;
            return a;
        }

        /** \brief return total number of cells in the logical Cartesian grid */
        long long cartesianSize() const
        {
            return 0;
        }

        /** \brief return number of cells in the active grid */
        long long compressedSize() const
        {
            return 0;
        }

        /** \brief return index of the cells in the logical Cartesian grid */
        long long cartesianIndex( const long long /* compressedElementIndex */) const
        {
            return 0;
        }

        /** \brief return Cartesian coordinate, i.e. IJK, for a given cell */
        void cartesianCoordinate(const long long /* compressedElementIndex */, std::array<long long,dimension>& /* coords */) const
        {
        }
    };

} // end namespace Opm
#endif
