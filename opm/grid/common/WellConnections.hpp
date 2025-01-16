/*
  Copyright 2016 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2016 Statoil AS

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
#ifndef DUNE_CPGRID_WELL_CONNECTIONS_HEADER_INCLUDED
#define DUNE_CPGRID_WELL_CONNECTIONS_HEADER_INCLUDED

#include <set>
#include <unordered_set>
#include <vector>
#include <functional>

#ifdef HAVE_MPI
#include <opm/grid/utility/platform_dependent/disable_warnings.h>
#include "mpi.h"
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>
#endif

#include <dune/common/parallel/communication.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/utility/OpmWellType.hpp>

namespace Dune
{
namespace cpgrid
{

/// \brief A class calculating and representing all connections of wells.
///
/// Wells are identified by their position as exported by the wells method
/// of the eclipse parser.
/// For each well the container stores at the well index all indices of cells
/// that the well perforates.
class WellConnections
{

public:
    /// \brief The const iterator type.
    typedef std::vector<std::set<long long> >::const_iterator const_iterator;

    /// \brief The iterator type (always const).
    typedef const_iterator iterator;

    WellConnections() = default;

    /// \brief Constructor
    /// \param wells The eclipse information about the wells
    /// \param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
    ///                                  The grid will then be partitioned such that these connections are on the same
    ///                                  partition. If NULL, they will be neglected.
    /// \param cartesianSize The logical cartesian size of the grid.
    /// \param cartesian_to_compressed Mapping of cartesian index
    ///        compressed cell index. The compressed index is used
    ///        to represent the well conditions.
    WellConnections(const std::vector<OpmWellType>& wells,
                    const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                    const std::array<long long, 3>& cartesianSize,
                    const std::vector<long long>& cartesian_to_compressed);

    /// \brief Constructor
    /// \param wells The eclipse information about the wells
    /// \param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
    ///                                  The grid will then be partitioned such that these connections are on the same
    ///                                  partition. If NULL, they will be neglected.
    /// \param cpGrid The corner point grid
    WellConnections(const std::vector<OpmWellType>& wells,
                    const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                    const Dune::CpGrid& cpGrid);

    /// \brief Initialze the data of the container
    /// \param wells The eclipse information about the wells
    /// \param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
    ///                                  The grid will then be partitioned such that these connections are on the same
    ///                                  partition. If NULL, they will be neglected.
    /// \param cartesianSize The logical cartesian size of the grid.
    /// \param cartesian_to_compressed Mapping of cartesian index
    ///        compressed cell index. The compressed index is used
    ///        to represent the well conditions.
    void init(const std::vector<OpmWellType>& wells,
              const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
              const std::array<long long, 3>& cartesianSize,
              const std::vector<long long>& cartesian_to_compressed);

    /// \brief Access all connections of a well
    /// \param i The index of the well (position of the well in the
    ///          eclipse schedule.
    /// \return The set of compressed indices of cells perforated by the well.
    const std::set<long long>& operator[](std::size_t i) const
    {
        return well_indices_[i];
    }

    /// \brief Get a begin iterator
    const_iterator begin() const
    {
        return well_indices_.begin();
    }

    /// \brief Get the end iterator
    const_iterator end() const
    {
        return well_indices_.end();
    }

    /// \breif Get the number of wells
    std::size_t size() const
    {
        return well_indices_.size();
    }
private:
    /// Stores at index i all cells that are perforated by
    /// the well at position i of the eclipse schedule.
    std::vector<std::set<long long> > well_indices_;
};


#ifdef HAVE_MPI
/// \brief Determines the wells that have perforate cells for each process.
///
/// On the root process omputes for all processes all indices of wells that
/// will perforate local cells.
/// Note that a well might perforate local cells of multiple processes
///
/// \param parts The partition number for each cell
/// \param wells The eclipse information about the wells.
/// \param possibleFutureConnections Possible future connections of wells that might get added through an ACTIONX.
///                                  The grid will then be partitioned such that these connections are on the same
///                                  partition. If NULL, they will be neglected.
/// \param cpGrid The unbalanced grid we compute on.
/// \return On the rank that has the global grid a vector with the well
///         indices for process i at index i.
std::vector<std::vector<long long> >
perforatingWellIndicesOnProc(const std::vector<long long>& parts,
                  const std::vector<Dune::cpgrid::OpmWellType>& wells,
                  const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                  const CpGrid& cpgrid);

/// \brief Computes wells assigned to processes.
///
/// Computes for all processes all indices of wells that
/// will be assigned to this process.
/// \param parts The partition number for each cell
/// \param gid Functor that turns cell index to global id.
/// \param eclipseState The eclipse information
/// \param well_connecton The information about the perforations of each well.
/// \param exportList List of cells to be exported. Each entry is a tuple of the
///                   global cell index, the process rank to export to, and the
///                   attribute on that rank (assumed to be owner)
/// \param exportList List of cells to be imported. Each entry is a tuple of the
///                   global cell index, the process rank that exports it, and the
///                   attribute on this rank (assumed to be owner)
/// \param cc Information about the parallelism together with the decomposition.
/// \return On rank 0 a vector with the well indices for process i
///         at index i.
std::vector<std::vector<long long> >
postProcessPartitioningForWells(std::vector<long long>& parts,
                                std::function<(long long)(long long)> gid,
                                const std::vector<OpmWellType>&  wells,
                                const WellConnections& well_connections,
                                const std::vector<std::set<long long> >& wellGraph,
                                std::vector<std::tuple<long long,long long,char>>& exportList,
                                std::vector<std::tuple<long long,long long,char,long long>>& importList,
                                const Communication<MPI_Comm>& cc);

/// \brief Computes whether wells are perforating cells on this process.
/// \param wells_on_proc well indices assigned to each process
/// \param eclipseState The eclipse information
/// \param cc The communicator
/// \param root The rank of the process that has the complete partitioning
///             information.
/// \return Vector of pairs of well name and a boolean indicating whether the
///         well with this name perforates cells here. Sorted by well name!
std::vector<std::pair<std::string,bool>>
computeParallelWells(const std::vector<std::vector<long long> >& wells_on_proc,
                     const std::vector<OpmWellType>&  wells,
                     const Communication<MPI_Comm>& cc,
                     long long root);
#endif
} // end namespace cpgrid
} // end namespace Dune

#endif
