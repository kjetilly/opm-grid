//===========================================================================
//
// File: LevelCartesianIndexMapper.hpp
//
// Created: Tue October 01  11:44:00 2024
//
// Author(s): Antonella Ritorto <antonella.ritorto@opm-op.com>
//
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2024 Equinor ASA.

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

#ifndef OPM_LEVELCARTESIANINDEXMAPPER_HEADER
#define OPM_LEVELCARTESIANINDEXMAPPER_HEADER

#include <array>

namespace Opm
{
// Interface class to access the local Cartesian grid of each level grid (when refinement).
template< class Grid >
class LevelCartesianIndexMapper
{
public:
    // Dimension of the grid.
    static constexpr long long dimension = Grid :: dimension ;

    // Constructor taking a grid.
    explicit LevelCartesianIndexMapper( const Grid& )
    {}

    // Return the number of cells in each direction (Cartesian dimensions) of a local Cartesian grid with level "level"
    const std::array<long long, dimension>& cartesianDimensions(long long /*level*/) const
    {
        static std::array<long long, dimension> a;
        return a;
    }

    // Return total number of cells in a local Cartesian grid with level "level".
    long long cartesianSize(long long /*level*/) const
    {
        return 0;
    }

    // Return number of cells in the active local Cartesian grid with level "level".
    long long compressedSize(long long /*level*/) const
    {
        return 0;
    }

    // Return index of a cell in the local Cartesian grid with level "level".
    long long cartesianIndex( const long long /* compressedElementIndex */ , const long long /*level*/) const
    {
        return 0;
    }

    // Compute Cartesian coordinate, i.e. IJK, for a given cell, on a given local Cartesian grid with level "level".
    void cartesianCoordinate(const long long /* compressedElementIndexOnLevel */,
                             std::array<long long,dimension>& /* coordsOnLevel */,
                             long long /*level*/) const
    {
    }
};

} // end namespace Opm
#endif
