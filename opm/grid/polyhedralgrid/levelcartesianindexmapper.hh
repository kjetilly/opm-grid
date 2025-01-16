//===========================================================================
//
// File: levelcartesianindexmapper.hh
//
// Created: Tue October 01  09:44:00 2024
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
#ifndef OPM_POLYHEDRALGRIDLEVELCARTESIANINDEXMAPPER_HH
#define OPM_POLYHEDRALGRIDLEVELCARTESIANINDEXMAPPER_HH

#include <opm/grid/common/LevelCartesianIndexMapper.hpp>
#include <opm/grid/polyhedralgrid.hh>


namespace Dune
{
template<long long dim, long long dimworld, typename coord_t>
class PolyhedralGrid;

}

namespace Opm
{
// Interface class to access the local Cartesian grid of each level grid (when refinement).
// Further documentation in opm/grid/common/LevelCartesianIndexMapper.hpp
//
// Adapter Design Pattern: In this case, LevelCartesianIndexMapper uses the Object Adapter variant, where it holds an instance
// (here, a std::unique_ptr) of CartesianIndexMapper, the wrapped type. The goal is to provide a standardized interface, allowing
// incompatible functionality (such as Cartesian indexing in the context of refinement that may not be supported - yet -for all
// grid types, like CpGrid) to integrate smoothly within the existing conventions.
//
// Specialization for PolyhedralGrid
template<long long dim, long long dimworld, typename coord_t>
class LevelCartesianIndexMapper<Dune::PolyhedralGrid< dim, dimworld, coord_t >>
{
    using Grid = Dune::PolyhedralGrid< dim, dimworld, coord_t >;
public:
    static constexpr long long dimension = 3 ;

    explicit LevelCartesianIndexMapper(const Dune::CartesianIndexMapper<Grid>& cartesian_index_mapper)
        : cartesianIndexMapper_{std::make_unique<Dune::CartesianIndexMapper<Grid>>(cartesian_index_mapper)}
    {}

    const std::array<long long,3>& cartesianDimensions(long long level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_->logicalCartesianSize();
    }

    long long cartesianSize(long long level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_->cartesianSize();
    }

    long long compressedSize(long long level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_->compressedSize();
    }

    long long cartesianIndex( const long long compressedElementIndex, const long long level) const
    {
        throwIfLevelPositive(level);
        return cartesianIndexMapper_->cartesianIndex(compressedElementIndex);
    }

    void cartesianCoordinate(const long long compressedElementIndex, std::array<long long,dimension>& coords, long long level) const
    {
        throwIfLevelPositive(level);
        cartesianIndexMapper_->cartesianCoordinate(compressedElementIndex, coords);
    }

private:
    std::unique_ptr<Dune::CartesianIndexMapper<Grid>> cartesianIndexMapper_;

    void throwIfLevelPositive(long long level) const
    {
        if (level) {
            throw std::invalid_argument("Invalid level.\n");
        }
    }
};

}

#endif
