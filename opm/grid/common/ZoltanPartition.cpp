/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.
  Copyright 2015 Statoil AS

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if HAVE_MPI // no code in this file without MPI. Skip includes-
#include <opm/grid/common/ZoltanPartition.hpp>
#include <opm/grid/utility/OpmWellType.hpp>
#include <opm/grid/cpgrid/CpGridData.hpp>
#include <opm/grid/cpgrid/Entity.hpp>
#include <algorithm>
#include <type_traits>
#endif

namespace Dune
{
namespace cpgrid
{

#if HAVE_MPI
template<class Id>
std::tuple<std::vector<long long>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<long long,long long,char> >,
           std::vector<std::tuple<long long,long long,char,long long> >,
           WellConnections>
makeImportAndExportLists(const Dune::CpGrid& cpgrid,
                         const Dune::Communication<MPI_Comm>& cc,
                         const std::vector<Dune::cpgrid::OpmWellType> * wells,
                         const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                         const Dune::cpgrid::CombinedGridWellGraph* gridAndWells,
                         long long root,
                         long long numExport,
                         long long numImport,
                         const Id* exportLocalGids,
                         const Id* exportGlobalGids,
                         const long long* exportToPart,
                         const Id* importGlobalGids,
                         bool allowDistributedWells) {

    long long                         size = cpgrid.numCells();
    long long                         rank  = cc.rank();
    std::vector<long long>            parts(size, rank);
    std::vector<std::vector<long long> > wellsOnProc;

    // List entry: process to export to, (global) index, process rank, attribute there (not needed?)
    std::vector<std::tuple<long long,long long,char>> myExportList;
    // List entry: process to import from, global index, process rank, attribute here, local index
    // (determined later)
    std::vector<std::tuple<long long,long long,char,long long>> myImportList;
    assert(rank==root || numExport==0);
    assert(rank!=root || numImport==0);
    constexpr double buffer = 1.05;
    // all cells on root are added to its export and its import list
    // myImportList will be expanded if distributed wells are not allowed and well cells are moved
    std::size_t reserveEx = rank!=root ? 0 : cpgrid.size(0);
    std::size_t reserveIm = rank!=root ? buffer*numImport : cpgrid.size(0)*buffer/cc.size();
    myExportList.reserve(reserveEx);
    myImportList.reserve(reserveIm);
    using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;

    for ( long long i=0; i < numExport; ++i )
    {
        parts[exportLocalGids[i]] = exportToPart[i];
        myExportList.emplace_back(exportGlobalGids[i], exportToPart[i], static_cast<char>(AttributeSet::owner));
    }

    for ( long long i=0; i < numImport; ++i )
    {
        myImportList.emplace_back(importGlobalGids[i], root, static_cast<char>(AttributeSet::owner),-1);
    }


    // Add cells that stay here to the lists. Somehow I could not persuade Zoltan to do this.
    for ( std::size_t i = 0; i < parts.size(); ++i)
    {
        if ( parts[i] == rank )
        {
            myExportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner) );
            myImportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner), -1 );
        }
    }
    std::inplace_merge(myImportList.begin(), myImportList.begin() + numImport, myImportList.end());
    std::inplace_merge(myExportList.begin(), myExportList.begin() + numExport, myExportList.end());




    if( wells )
    {
        if (allowDistributedWells)
        {
            wellsOnProc = perforatingWellIndicesOnProc(parts, *wells, possibleFutureConnections, cpgrid);
        }
        else
        {
            auto gidGetter = [&cpgrid](long long i) { return cpgrid.globalIdSet().id(Dune::createEntity<0>(cpgrid, i, true));};
            wellsOnProc =
                postProcessPartitioningForWells(parts,
                                                gidGetter,
                                                *wells,
                                                gridAndWells->getWellConnections(),
                                                gridAndWells->getWellsGraph(),
                                                myExportList, myImportList,
                                                cc);


#ifndef NDEBUG
            std::size_t index = 0;
            std::unordered_set<long long> distributed_wells;

            for( auto well : gridAndWells->getWellsGraph() )
            {
                long long part=parts[index];
                std::set<std::pair<long long,long long> > cells_on_other;
                for( auto vertex : well )
                {
                    if( part != parts[vertex] )
                    {
                        cells_on_other.insert(std::make_pair(vertex, parts[vertex]));
                    }
                }
                if ( cells_on_other.size() )
                {
                    distributed_wells.insert(index);
                }
                ++index;
            }
            auto num_dist_wells = cc.sum(distributed_wells.size());
            if (num_dist_wells) {
                OPM_THROW(std::domain_error,
                          std::to_string(num_dist_wells) + " well"
                              +  ((num_dist_wells >= 2) ? "s are" : " is")
                              + " distributed between processes, which should not be the case!");
            }


#endif
        }
    }

    std::vector<std::pair<std::string,bool>> parallel_wells;
    if( wells )
    {
        parallel_wells = Dune::cpgrid::computeParallelWells(wellsOnProc,
                                                            *wells,
                                                            cc,
                                                            root);
    }
    return std::make_tuple(parts, parallel_wells, myExportList, myImportList,
                           wells ? gridAndWells->getWellConnections(): WellConnections());
}

template<class Id>
std::tuple<long long, std::vector<Id> >
scatterExportInformation(long long numExport, const Id* exportGlobalGids,
                         const long long* exportToPart, long long root,
                         const Communication<MPI_Comm>& cc)
{
    long long numImport;
    std::vector<long long> numberOfExportedVerticesPerProcess;
    std::vector<Id> importGlobalGidsVector;
    std::vector<Id> globalIndicesToSend;
    std::vector<long long> offsets;

    // Build and communicate import/export data.
    // 1. Send number of exports/imports.
    if (cc.rank() == root) {
        numberOfExportedVerticesPerProcess.resize(cc.size(), 0);
        for (long long i = 0; i < numExport; ++i) {
            ++numberOfExportedVerticesPerProcess[exportToPart[i]];
        }
    }

    cc.scatter(numberOfExportedVerticesPerProcess.data(), &numImport, 1, root);
    assert(cc.rank() != root || numImport == 0);

    // 2. Build the imports/exports themselves.
    std::string error;
    if (cc.rank() == root) {
        offsets.resize(cc.size() + 1, 0);
        std::partial_sum(numberOfExportedVerticesPerProcess.begin(),
                         numberOfExportedVerticesPerProcess.end(),
                         offsets.begin() + 1);
            globalIndicesToSend.resize(numExport, 0);
            const long long commSize = cc.size();
            std::vector<long long> currentIndex(commSize, 0);
            for (long long i = 0; i < numExport; ++i) {
                if (exportToPart[i] >= commSize) {
                    std::ostringstream oss;
                    oss << "Something wrong with Zoltan decomposition. "
                        << "Debug information: exportToPart[i] = " << exportToPart[i] << ", "
                        << "currentIndex.size() = " << currentIndex.size();
                    error = oss.str();
                    break;
                }
                const auto index = currentIndex[exportToPart[i]]++ + offsets[exportToPart[i]];
                if (index >= numExport) {
                    std::ostringstream oss;
                    oss << "Something wrong with Zoltan decomposition. "
                        << "index " << index << ", "
                        << "globalIndicesToSend.size() = " << globalIndicesToSend.size()
                        << "\n\noffsets[exportToPart[i]] = " << offsets[exportToPart[i]]
                        << "\n\ncurrentIndex[exportToPart[i]] = " << currentIndex[exportToPart[i]]
                        << "\n\nexportToPart[i] = " << exportToPart[i];
                    error = oss.str();
                    break;
                }
                globalIndicesToSend[index] = exportGlobalGids[i];
            }
        } else {
            importGlobalGidsVector.resize(numImport, 0);
        }
        // Check for errors
        long long ok = error.empty();
        cc.broadcast(&ok, 1, root);
        if (!ok) {
            OPM_THROW(std::runtime_error, error);
        }

        // 3. Communicate the imports/exports.
        cc.scatterv(globalIndicesToSend.data(), numberOfExportedVerticesPerProcess.data(),
                    offsets.data(), importGlobalGidsVector.data(), numImport, root);
        return std::make_tuple(numImport, importGlobalGidsVector);
}

// instantiate long long types
template
std::tuple<std::vector<long long>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<long long,long long,char> >,
           std::vector<std::tuple<long long,long long,char,long long> >,
           WellConnections>
makeImportAndExportLists(const Dune::CpGrid&,
                         const Communication<MPI_Comm>&,
                         const std::vector<Dune::cpgrid::OpmWellType>*,
                         const std::unordered_map<std::string, std::set<long long>>&,
                         const Dune::cpgrid::CombinedGridWellGraph*,
                         long long,
                         long long,
                         long long,
                         const long long*,
                         const long long*,
                         const long long*,
                         const long long*,
                         bool);

template
std::tuple<long long, std::vector<long long> >
scatterExportInformation(long long numExport, const long long*,
                         const long long*, long long,
                         const Communication<MPI_Comm>&);
#endif // HAVE_MPI
} // end namespace Dune
} // end namespace cpgrid

#if defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
namespace Dune
{
namespace cpgrid
{

namespace {
void setDefaultZoltanParameters(Zoltan_Struct* zz) {
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
    Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","0");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */
}


} // anon namespace

std::tuple<std::vector<long long>, std::vector<std::pair<std::string,bool>>,
           std::vector<std::tuple<long long,long long,char> >,
           std::vector<std::tuple<long long,long long,char,long long> >,
           WellConnections>
zoltanGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                               const std::vector<OpmWellType> * wells,
                               const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                               const double* transmissibilities,
                               const Communication<MPI_Comm>& cc,
                               EdgeWeightMethod edgeWeightsMethod,
                               long long root,
                               const double zoltanImbalanceTol,
                               bool allowDistributedWells,
                               const std::map<std::string,std::string>& params)
{
    long long rc = ZOLTAN_OK - 1;
    float ver = 0;
    struct Zoltan_Struct *zz;
    long long changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    long long *importProcs, *importToPart, *exportProcs, *exportToPart;
    long long argc=0;
    char** argv = 0 ;
    rc = Zoltan_Initialize(argc, argv, &ver);
    zz = Zoltan_Create(cc);
    if ( rc != ZOLTAN_OK )
    {
        OPM_THROW(std::runtime_error, "Could not initialize Zoltan!");
    }
    setDefaultZoltanParameters(zz);
    Zoltan_Set_Param(zz, "IMBALANCE_TOL", std::to_string(zoltanImbalanceTol).c_str());
    for (const auto& [key, value] : params)
        Zoltan_Set_Param(zz, key.c_str(), value.c_str());

    // For the load balancer one process has the whole grid and
    // all others an empty partition before loadbalancing.
    bool partitionIsEmpty     = cc.rank()!=root;

    std::shared_ptr<CombinedGridWellGraph> gridAndWells;

    if( wells )
    {
        Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
        gridAndWells.reset(new CombinedGridWellGraph(cpgrid,
                                                       wells,
                                                       possibleFutureConnections,
                                                       transmissibilities,
                                                       partitionIsEmpty,
                                                       edgeWeightsMethod));
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *gridAndWells,
                                                    partitionIsEmpty);
    }
    else
    {
        Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, cpgrid, partitionIsEmpty);
    }

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */

    auto importExportLists = makeImportAndExportLists(cpgrid,
                                     cc,
                                     wells,
                                     possibleFutureConnections,
                                     gridAndWells.get(),
                                     root,
                                     numExport,
                                     numImport,
                                     exportLocalGids,
                                     exportGlobalGids,
                                     exportProcs,
                                     importGlobalGids,
                                     allowDistributedWells);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);

    return importExportLists;
}



class ZoltanSerialPartitioner
{
public:
    using CommunicationType = Communication<MPI_Comm>;

    ZoltanSerialPartitioner(const CpGrid& _cpgrid,
                            const std::vector<OpmWellType>* _wells,
                            const std::unordered_map<std::string, std::set<long long>>& _possibleFutureConnections,
                            const double* _transmissibilities,
                            const CommunicationType& _cc,
                            EdgeWeightMethod _edgeWeightsMethod,
                            long long _root,
                            const double _zoltanImbalanceTol,
                            bool _allowDistributedWells,
                            long long _numParts,
                            const std::map<std::string,std::string>& param)
        : cpgrid(_cpgrid)
        , wells(_wells)
        , possibleFutureConnections(_possibleFutureConnections)
        , transmissibilities(_transmissibilities)
        , cc(_cc)
        , edgeWeightsMethod(_edgeWeightsMethod)
        , root(_root)
        , zoltanImbalanceTol(_zoltanImbalanceTol)
        , allowDistributedWells(_allowDistributedWells)
        , numParts(_numParts)
        , params(param)
    {
        if (wells) {
            const bool partitionIsEmpty = cc.rank() != root;
            gridAndWells.reset(
                new CombinedGridWellGraph(cpgrid, wells, possibleFutureConnections, transmissibilities, partitionIsEmpty, edgeWeightsMethod));
        }
    }


    std::vector<long long> partitionForInfo()
    {
        if (cc.rank() == root) {
            callZoltan();
        }

        long long size = cpgrid.numCells();
        std::vector<long long> parts(size, 0);

        for ( long long i=0; i < numExport; ++i ) {
            parts[i] = exportToPart[i];
        }
        return parts;
    }
    std::tuple<std::vector<long long>,
               std::vector<std::pair<std::string, bool>>,
               std::vector<std::tuple<long long, long long, char>>,
               std::vector<std::tuple<long long, long long, char, long long>>,
               WellConnections>
    partition()
    {
        MPI_Barrier(cc);

        // Initialize Zoltan and perform partitioning.
        long long rc = ZOLTAN_OK;
        if (cc.rank() == root) {
            rc = callZoltan();
        }

        cc.broadcast(&rc, 1, root);
        if (rc != ZOLTAN_OK) {
            OPM_THROW(std::runtime_error, "Could not initialize Zoltan, or Zoltan partitioning failed.");
        }

        std::tie(numImport, importGlobalGidsVector) = scatterExportInformation<ZoltanId>(numExport, exportGlobalGids,
                                                                                         exportToPart, root, cc);

        if (cc.rank() != root)
            importGlobalGids = importGlobalGidsVector.data();

        return makeImportAndExportLists(cpgrid,
                                        cc,
                                        wells,
                                        possibleFutureConnections,
                                        gridAndWells.get(),
                                        root,
                                        numExport,
                                        numImport,
                                        exportLocalGids,
                                        exportGlobalGids,
                                        exportToPart,
                                        importGlobalGids,
                                        allowDistributedWells);
    }

    ~ZoltanSerialPartitioner()
    {
        // free space allocated for zoltan.
        if (cc.rank() == root) {
            Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
            Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
            Zoltan_Destroy(&zz);
        }
    }

private:
    // Methods

    long long callZoltan()
    {
        long long argc = 0;
        char** argv = 0;
        float ver = 0;

        long long rc = Zoltan_Initialize(argc, argv, &ver);
        zz = Zoltan_Create(MPI_COMM_SELF);
        if (rc != ZOLTAN_OK) {
            return rc;
        }

        setDefaultZoltanParameters(zz);
        Zoltan_Set_Param(zz, "IMBALANCE_TOL", std::to_string(zoltanImbalanceTol).c_str());
        if (numParts > 0) {
            Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(numParts).c_str());
            Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
            Zoltan_Set_Param(zz, "DEBUG_LEVEL", "2");
        }
        else
            Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(cc.size()).c_str());

        for (const auto& [key, value] : params)
            Zoltan_Set_Param(zz, key.c_str(), value.c_str());

        // For the load balancer one process has the whole grid and
        // all others an empty partition before loadbalancing.
        bool partitionIsEmpty = cc.rank() != root;

        if (wells) {
            Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");
            Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, *gridAndWells, partitionIsEmpty);
        } else {
            Dune::cpgrid::setCpGridZoltanGraphFunctions(zz, cpgrid, partitionIsEmpty);
        }

        rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                                 &changes, /* 1 if partitioning was changed, 0 otherwise */
                                 &numGidEntries, /* Number of integers used for a global ID */
                                 &numLidEntries, /* Number of integers used for a local ID */
                                 &numImport, /* Number of vertices to be sent to me */
                                 &importGlobalGids, /* Global IDs of vertices to be sent to me */
                                 &importLocalGids, /* Local IDs of vertices to be sent to me */
                                 &importProcs, /* Process rank for source of each incoming vertex */
                                 &importToPart, /* New partition for each incoming vertex */
                                 &numExport, /* Number of vertices I must send to other processes*/
                                 &exportGlobalGids, /* Global IDs of the vertices I must send */
                                 &exportLocalGids, /* Local IDs of the vertices I must send */
                                 &exportProcs, /* Process to which I send each of the vertices */
                                 &exportToPart); /* Partition to which each vertex will belong */

        // This is very important: by default, zoltan sets numImport to -1,
        // as we are only running Zoltan in serial.
        // In order to make sense of the rest of  the code, this must be set
        // to 0.
        numImport = 0;
        return rc;
    }

    // Data members
    const CpGrid& cpgrid;
    const std::vector<OpmWellType>* wells;
    const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections;
    const double* transmissibilities;
    const CommunicationType& cc;
    EdgeWeightMethod edgeWeightsMethod;
    long long root;
    const double zoltanImbalanceTol;
    std::string errorOnRoot;

    struct Zoltan_Struct* zz = nullptr;
    long long changes = 0;
    long long numGidEntries = 0;
    long long numLidEntries = 0;
    long long numImport = 0;
    long long numExport = 0;
    ZOLTAN_ID_PTR importGlobalGids = nullptr;
    ZOLTAN_ID_PTR importLocalGids = nullptr;
    ZOLTAN_ID_PTR exportGlobalGids = nullptr;
    ZOLTAN_ID_PTR exportLocalGids = nullptr;
    long long *importProcs, *importToPart, *exportProcs, *exportToPart;
    std::unique_ptr<CombinedGridWellGraph> gridAndWells;
    using ZoltanId = typename std::remove_pointer<ZOLTAN_ID_PTR>::type;
    std::vector<ZoltanId> importGlobalGidsVector;
    bool allowDistributedWells;
    long long numParts;
    const std::map<std::string,std::string>& params;
};


std::tuple<std::vector<long long>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<long long, long long, char>>,
           std::vector<std::tuple<long long, long long, char, long long>>,
           WellConnections>
zoltanSerialGraphPartitionGridOnRoot(const CpGrid& cpgrid,
                                     const std::vector<OpmWellType>* wells,
                                     const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                                     const double* transmissibilities,
                                     const Dune::Communication<MPI_Comm>& cc,
                                     EdgeWeightMethod edgeWeightsMethod,
                                     long long root,
                                     const double zoltanImbalanceTol,
                                     bool allowDistributedWells,
                                     const std::map<std::string,std::string>& params)
{
    ZoltanSerialPartitioner partitioner(cpgrid, wells, possibleFutureConnections, transmissibilities, cc, edgeWeightsMethod,
                                        root, zoltanImbalanceTol, allowDistributedWells, 0, params);
    return partitioner.partition();
}

std::vector<long long>
zoltanGraphPartitionGridForJac(const CpGrid& cpgrid,
                               const std::vector<OpmWellType> * wells,
                               const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                               const double* transmissibilities,
                               const Dune::Communication<MPI_Comm>& cc,
                               EdgeWeightMethod edgeWeightsMethod, long long root,
                               long long numParts, const double zoltanImbalanceTol)
{
    // Parameters are empty, but must still have scope here since it
    // (or rather a const reference to it) is queried in
    // partitionForInfo() further down.
    std::map<std::string,std::string> params;

    ZoltanSerialPartitioner partitioner(cpgrid, wells, possibleFutureConnections, transmissibilities, cc, edgeWeightsMethod,
                                        root, zoltanImbalanceTol, false, numParts, params);
    return partitioner.partitionForInfo();
}


} // namespace cpgrid
} // namespace Dune
#endif // defined(HAVE_ZOLTAN) && defined(HAVE_MPI)
