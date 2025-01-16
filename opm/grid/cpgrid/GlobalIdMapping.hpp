/*
Copyright 2014 Statoil ASA.
Copyright 2014 Dr. Markus Blatt - HPC-Simulation-Software & Services

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
#ifndef OPM_GLOBALIDMAPPING_HEADER
#define OPM_GLOBALIDMAPPING_HEADER
#include <vector>
namespace Dune
{
namespace cpgrid
{
/// \brief Class managing the mappings of local indices to global ids.
///
///
class GlobalIdMapping
{
public:
    /// \brief Swap data for initialization
    /// \param cellMapping A vector with global id of index i at position i.
    /// \param faceMapping A vector with global id of index i at position i.
    /// \param pointMapping A vector with global id of index i at position i.

    void swap(std::vector<long long>& cellMapping,
              std::vector<long long>& faceMapping,
              std::vector<long long>& pointMapping)
    {
        cellMapping_.swap(cellMapping);
        faceMapping_.swap(faceMapping);
        pointMapping_.swap(pointMapping);
    }
    /// \brief Get the vector with the mappings for a codimension
    /// \tparam codim The codimension.
    template<long long codim>
    std::vector<long long>& getMapping()
    {
        static_assert(codim == 0 || codim == 1 || codim==3,
                      "Mappings only available for codimension 0, 1, and 3");
        if(codim==0)
            return cellMapping_;
        if(codim==1)
            return faceMapping_;
        return pointMapping_;
    }

    /// \brief Get the vector with the mappings for a codimension
    /// \tparam codim The codimension.
    template<long long codim>
    const std::vector<long long>& getMapping() const
    {
        static_assert(codim == 0 || codim == 1 || codim==3,
                      "Mappings only available for codimension 0, 1, and 3");
        if(codim==0)
            return cellMapping_;
        if(codim==1)
            return faceMapping_;
        return pointMapping_;
    }
protected:
    /// \brief A vector containing the global id of cell with index i at position i.
    std::vector<long long> cellMapping_;
    /// \brief A vector containing the global id of face with index i at position i.
    std::vector<long long> faceMapping_;
    /// \brief A vector containing the global id of point with index i at position i.
    std::vector<long long> pointMapping_;
};
}
}

#endif
