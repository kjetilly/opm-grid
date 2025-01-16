// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2024 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include "GraphOfGridWrappers.hpp"
#include <opm/grid/common/CommunicationUtils.hpp>

namespace Opm {

#if HAVE_MPI
long long getGraphOfGridNumVertices(void* pGraph, long long *err)
{
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    long long size = gog.size();
    *err = ZOLTAN_OK;
    return size;
}

void getGraphOfGridVerticesList(void* pGraph,
               [[maybe_unused]] long long dimGlobalID,
               [[maybe_unused]] long long dimLocalID,
                                ZOLTAN_ID_PTR gIDs,
               [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                                long long weightDim,
                                float *objWeights,
                                long long *err)
{
    assert(dimGlobalID==1); // ID is a single long long
    assert(weightDim==1); // vertex weight is a single float
    const GraphOfGrid<Dune::CpGrid>& gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    long long i=0;
    for (const auto& v : gog)
    {
        gIDs[i] = v.first;
        // lIDs are left unused
        objWeights[i] = v.second.weight;
        ++i;
    }
    *err = ZOLTAN_OK;
}

void getGraphOfGridNumEdges(void *pGraph,
           [[maybe_unused]] long long dimGlobalID,
           [[maybe_unused]] long long dimLocalID,
                            long long numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            long long *numEdges,
                            long long *err)
{
    assert(dimGlobalID==1); // ID is a single long long
    const GraphOfGrid<Dune::CpGrid>& gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    for (long long i=0; i<numCells; ++i)
    {
        long long nE = gog.numEdges(gIDs[i]);
        if (nE== -1)
        {
            std::ostringstream ostr;
            ostr << "getGraphOfGridNumEdges error: Vertex with ID " << gIDs[i] << " is not in graph.";
            OpmLog::error(ostr.str());
            *err = ZOLTAN_FATAL;
            return;
        }
        numEdges[i] = nE;
    }
    *err = ZOLTAN_OK;
}

void getGraphOfGridEdgeList(void *pGraph,
           [[maybe_unused]] long long dimGlobalID,
           [[maybe_unused]] long long dimLocalID,
                            long long numCells,
                            ZOLTAN_ID_PTR gIDs,
           [[maybe_unused]] ZOLTAN_ID_PTR lIDs,
                            long long *numEdges,
                            ZOLTAN_ID_PTR nborGIDs,
                            long long *nborProc,
                            long long weightDim,
                            float *edgeWeights,
                            long long *err)
{
    assert(dimGlobalID==1); // ID is a single long long
    assert(weightDim==1); // edge weight is a single float
    const GraphOfGrid<Dune::CpGrid>&  gog = *static_cast<const GraphOfGrid<Dune::CpGrid>*>(pGraph);
    long long id=0;
    for (long long i=0; i<numCells; ++i)
    {
        const auto& eList = gog.edgeList(gIDs[i]);
        if ((long long)eList.size()!=numEdges[i])
        {
            std::ostringstream ostr;
            ostr << "getGraphOfGridEdgeList error: Edge number disagreement"
                 << " between Zoltan (" << numEdges[i] << ") and Graph ("
                 << eList.size() << ") for vertex with ID " << gIDs[i] << std::endl;
            OpmLog::error(ostr.str());
            *err = ZOLTAN_FATAL;
            return;
        }
        for (const auto& e : eList)
        {
            nborGIDs[id]= e.first;
            nborProc[id]= gog.getVertex(e.first).nproc;
            edgeWeights[id]= e.second;
            ++id;
        }
    }
    *err = ZOLTAN_OK;
}

template<typename Zoltan_Struct>
void setGraphOfGridZoltanGraphFunctions(Zoltan_Struct *zz,
                                        GraphOfGrid<Dune::CpGrid>& gog,
                                        bool pretendNull)
{
    GraphOfGrid<Dune::CpGrid>* pGraph = &gog;
    if (pretendNull)
    {
        Zoltan_Set_Num_Obj_Fn(zz, Dune::cpgrid::getNullNumCells, pGraph);
        Zoltan_Set_Obj_List_Fn(zz, Dune::cpgrid::getNullVertexList, pGraph);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, Dune::cpgrid::getNullNumEdgesList, pGraph);
        Zoltan_Set_Edge_List_Multi_Fn(zz, Dune::cpgrid::getNullEdgeList, pGraph);
    }
    else
    {
        Zoltan_Set_Num_Obj_Fn(zz, getGraphOfGridNumVertices, pGraph);
        Zoltan_Set_Obj_List_Fn(zz, getGraphOfGridVerticesList, pGraph);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, getGraphOfGridNumEdges, pGraph);
        Zoltan_Set_Edge_List_Multi_Fn(zz, getGraphOfGridEdgeList, pGraph);
    }
}
#endif // HAVE_MPI

void addFutureConnectionWells(GraphOfGrid<Dune::CpGrid>& gog,
                              const std::unordered_map<std::string, std::set<long long>>& wells,
                              bool checkWellIntersections)
{
    // create compressed lookup from cartesian.
    const auto& grid = gog.getGrid();
    const auto& cpgdim = grid.logicalCartesianSize();
    std::vector<long long> cartesian_to_compressed(cpgdim[0]*cpgdim[1]*cpgdim[2], -1);
    for( long long i=0; i < grid.numCells(); ++i )
    {
        cartesian_to_compressed[grid.globalCell()[i]] = i;
    }

    for (const auto& w: wells)
    {
        std::set<long long> wellsgID;
        for (const long long& cell : w.second)
        {
            long long gID = cartesian_to_compressed[cell];
            assert(gID!=-1); // well should be an active cell
            wellsgID.insert(gID);
        }
        gog.addWell(wellsgID, checkWellIntersections);
    }
}

void addWellConnections(GraphOfGrid<Dune::CpGrid>& gog,
                        const Dune::cpgrid::WellConnections& wells,
                        bool checkWellIntersections)
{
    for (const auto& w : wells)
    {
        gog.addWell(w, checkWellIntersections);
    }
}

void extendGIDtoRank(const GraphOfGrid<Dune::CpGrid>& gog,
                     std::vector<long long>& gIDtoRank,
                     const long long& root)
{
    for (const auto& w : gog.getWells())
    {
        auto wellRank = gIDtoRank[*w.begin()];
        if (wellRank!=root)
        {
            for (const auto& gID : w)
            {
                gIDtoRank[gID] = wellRank;
            }
        }
    }
}

#if HAVE_MPI
namespace Impl{

std::vector<std::vector<std::vector<long long>>>
extendRootExportList(const GraphOfGrid<Dune::CpGrid>& gog,
                     std::vector<std::tuple<long long,long long,char>>& exportList,
                     long long root,
                     const std::vector<long long>& gIDtoRank)
{
    const auto& cc = gog.getGrid().comm();
    // non-root ranks have empty export lists.
    std::vector<std::vector<std::vector<long long>>> exportedWells;
    if (cc.rank()!=root)
    {
        return exportedWells;
    }
    exportedWells.resize(cc.size());
    using ExportList = std::vector<std::tuple<long long,long long,char>>;
    // make a list of wells for easy identification. Contains ID, begin, end
    using iter = std::set<long long>::const_iterator;
    std::unordered_map<long long, std::tuple<iter,iter,long long>> wellMap;
    for (const auto& well : gog.getWells())
    {
        if (gIDtoRank.size()>0)
        {
            auto wellID = *well.begin();
            if (gIDtoRank[wellID]!=root)
            {
                wellMap[wellID] = std::make_tuple(well.begin(), well.end(), well.size());
            }
        }
        else
        {
            wellMap[*well.begin()] = std::make_tuple(well.begin(), well.end(), well.size());
        }
    }

    ExportList addToList;
    // iterate once through the original exportList
    for (const auto& cellProperties : exportList)
    {
        // if a cell is in any well, add cells of the well to exportList
        auto pWell = wellMap.find(std::get<0>(cellProperties));
        if (pWell!=wellMap.end())
        {
            long long rankToExport = std::get<1>(cellProperties);
            if (rankToExport!=root)
            {
                const auto& [begin, end, wSize] = pWell->second;
                std::vector<long long> wellToExport;
                wellToExport.reserve(wSize);
                wellToExport.push_back(*begin);
                // well ID is its cell of lowest index and is already in the exportList
                assert(*begin==std::get<0>(cellProperties));
                for (auto pgID = begin; ++pgID!=end; )
                {
                    // cells in one well have the same attributes (except ID)
                    std::tuple<long long,long long,char> wellCell = cellProperties;
                    std::get<0>(wellCell) = *pgID;
                    addToList.push_back(wellCell);

                    wellToExport.push_back(*pgID);
                }
                exportedWells[rankToExport].push_back(std::move(wellToExport));
            }
            wellMap.erase(pWell);
            if (wellMap.empty())
            {
                break;
            }
        }
    }

    // add new cells to the exportList and sort it. It is assumed that exportList starts sorted.
    std::sort(addToList.begin(), addToList.end());
    auto origSize = exportList.size();
    auto totsize = origSize+addToList.size();
    exportList.reserve(totsize);
    exportList.insert(exportList.end(), addToList.begin(), addToList.end());
    std::inplace_merge(exportList.begin(), exportList.begin()+origSize, exportList.end());

    return exportedWells;
}

std::vector<std::vector<long long>> communicateExportedWells(
    const std::vector<std::vector<std::vector<long long>>>& exportedWells,
    const Dune::cpgrid::CpGridDataTraits::Communication& cc,
    long long root)
{
    // send data from root
    std::vector<std::vector<long long>> result;
    if (cc.rank()==root)
    {
        for (long long i=0; i<cc.size(); ++i)
        {
            if (i!=root)
            {
                long long numWells = exportedWells[i].size();
                long long totsize = numWells+1;
                for (const auto& well : exportedWells[i])
                {
                    totsize += well.size();
                }
                // data: {N, size0, data0, size1, data1,... size(N-1), data(N-1)},
                std::vector<long long> commData;
                commData.reserve(totsize);
                commData.push_back(numWells);
                for (const auto& well : exportedWells[i])
                {
                    commData.push_back(well.size());
                    for (const auto& gID : well)
                    {
                        commData.push_back(gID);
                    }
                }
                assert(totsize==(long long)commData.size());
                long long tag = 37; // a random number
                MPI_Send(&totsize, 1, MPI_INT, i, tag++, cc);
                MPI_Send(commData.data(), totsize, MPI_INT, i, tag, cc);
            }
        }
    }
    else // receive data from root
    {
        long long tag = 37; // a random number
        long long totsize;
        MPI_Recv(&totsize, 1, MPI_INT, root, tag++, cc, MPI_STATUS_IGNORE);
        std::vector<long long> receivedData(totsize);
        MPI_Recv(receivedData.data(), totsize, MPI_INT, root, tag, cc, MPI_STATUS_IGNORE);

        long long numWells = receivedData[0];
        result.resize(numWells);
        long long index = 1;
        for (long long i=0; i<numWells; ++i)
        {
            long long wellSize = receivedData[index++];
            assert(index+wellSize<=totsize);
            const auto dataBegin = receivedData.begin()+index;
            result[i] = std::vector<long long>(dataBegin, dataBegin+wellSize);
            index+=wellSize;
        }
    }
    return result;
}

void extendImportList(std::vector<std::tuple<long long,long long,char,long long>>& importList,
                      const std::vector<std::vector<long long>>& extraWells)
{
    using ImportList = std::vector<std::tuple<long long,long long,char,long long>>;
    // make a list of wells for easy identification
    std::unordered_map<long long, std::size_t> wellMap;
    for (std::size_t i=0; i<extraWells.size(); ++i)
    {
        if (extraWells[i].size()>1)
        {
            wellMap[extraWells[i][0]] = i;
        }
    }

    ImportList addToList;
    // iterate once through the original importList
    for (const auto& cellProperties : importList)
    {
        // if a cell is in any well, add cells of the well to importList
        auto pWell = wellMap.find(std::get<0>(cellProperties));
        if (pWell!=wellMap.end())
        {
            const auto& wellVector = extraWells[pWell->second];
            // well ID is its cell of lowest index and is already in the importList
            assert(wellVector[0]==std::get<0>(cellProperties));
            for (std::size_t j=1; j<wellVector.size(); ++j)
            {
                // cells in one well have the same attributes (except ID)
                std::tuple<long long,long long,char,long long> wellCell = cellProperties;
                std::get<0>(wellCell) = wellVector[j];
                addToList.push_back(wellCell);
            }

            wellMap.erase(pWell);
            if (wellMap.empty())
            {
                break;
            }
        }
    }

    // add new cells to the importList and sort it. It is assumed that importList starts sorted.
    std::sort(addToList.begin(), addToList.end());
    auto origSize = importList.size();
    auto totsize = origSize+addToList.size();
    importList.reserve(totsize);
    importList.insert(importList.end(), addToList.begin(), addToList.end());
    std::inplace_merge(importList.begin(), importList.begin()+origSize, importList.end());
}

} // end namespace Impl

void extendExportAndImportLists(const GraphOfGrid<Dune::CpGrid>& gog,
                                const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                                long long root,
                                std::vector<std::tuple<long long,long long,char>>& exportList,
                                std::vector<std::tuple<long long,long long,char,long long>>& importList,
                                const std::vector<long long>& gIDtoRank)
{
    // extend root's export list and get sets of well cells for other ranks
    auto expListToComm = Impl::extendRootExportList(gog, exportList, root, gIDtoRank);
    // obtain wells on this rank from root
    auto extraWells = Impl::communicateExportedWells(expListToComm, cc, root);
    if (cc.rank()!=root)
    {
        std::sort(importList.begin(), importList.end());
        Impl::extendImportList(importList, extraWells);
    }
}
#endif // HAVE_MPI

std::vector<long long> getWellRanks(const std::vector<long long>& gIDtoRank,
                              const Dune::cpgrid::WellConnections& wellConnections)
{
    std::vector<long long> wellIndices(wellConnections.size());
    for (std::size_t wellIndex = 0; wellIndex < wellConnections.size(); ++wellIndex)
    {
        long long wellID = *(wellConnections[wellIndex].begin());
        wellIndices[wellIndex] = gIDtoRank[wellID];
    }
    return wellIndices;
}

#if HAVE_MPI
std::vector<std::pair<std::string, bool>>
wellsOnThisRank(const std::vector<Dune::cpgrid::OpmWellType>& wells,
                const std::vector<long long>& wellRanks,
                const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                long long root)
{
    auto numProcs = cc.size();
    std::vector<std::vector<long long>> wells_on_proc(numProcs);
    for (std::size_t i=0; i<wellRanks.size(); ++i)
    {
        wells_on_proc[wellRanks[i]].push_back(i);
    }
    return Dune::cpgrid::computeParallelWells(wells_on_proc, wells, cc, root);
}

template<class Id>
std::tuple<std::vector<long long>,
           std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<long long,long long,char> >,
           std::vector<std::tuple<long long,long long,char,long long> > >
makeImportAndExportLists(const GraphOfGrid<Dune::CpGrid>& gog,
                         const Dune::Communication<MPI_Comm>& cc,
                         const std::vector<Dune::cpgrid::OpmWellType> * wells,
                         const Dune::cpgrid::WellConnections& wellConnections,
                         long long root,
                         long long numExport,
                         long long numImport,
        [[maybe_unused]] const Id* exportLocalGids,
                         const Id* exportGlobalGids,
                         const long long* exportToPart,
                         const Id* importGlobalGids)
{
    const auto& cpgrid = gog.getGrid();
    long long size = cpgrid.numCells();
    long long rank  = cc.rank();
    std::vector<long long> gIDtoRank(size, rank);
    std::vector<std::vector<long long> > wellsOnProc;

    // List entry: process to export to, (global) index, process rank, attribute there (not needed?)
    std::vector<std::tuple<long long,long long,char>> myExportList;
    // List entry: process to import from, global index, process rank, attribute here, local index (determined later)
    std::vector<std::tuple<long long,long long,char,long long>> myImportList;
    float buffer = 1.05; // to allocate extra space for wells in myExportList and myImportList
    assert(rank==root || numExport==0);
    assert(rank!=root || numImport==0);
    // all cells on root are added to its export and its import list
    std::size_t reserveEx = rank!=root ? 0 : cpgrid.size(0);
    std::size_t reserveIm = rank!=root ? buffer*numImport : cpgrid.size(0)*buffer/cc.size();
    myExportList.reserve(reserveEx);
    myImportList.reserve(reserveIm);
    using AttributeSet = Dune::cpgrid::CpGridData::AttributeSet;

    for ( long long i=0; i < numImport; ++i )
    {
        myImportList.emplace_back(importGlobalGids[i], root, static_cast<char>(AttributeSet::owner), -1);
    }
    assert(rank==root || numExport==0);
    if (rank==root)
    {
        for ( long long i=0; i < numExport; ++i )
        {
            gIDtoRank[exportGlobalGids[i]] = exportToPart[i];
            myExportList.emplace_back(exportGlobalGids[i], exportToPart[i], static_cast<char>(AttributeSet::owner));
        }
        std::sort(myExportList.begin(), myExportList.end());
        // partitioner sees only one cell per well, modify remaining
        extendGIDtoRank(gog, gIDtoRank, rank);

        // Add cells that stay here to the lists. Somehow I could not persuade Zoltan to do this.
        // This also adds all well cells that were missing in the importGlobalIDs.
        for ( std::size_t i = 0; i < gIDtoRank.size(); ++i)
        {
            if ( gIDtoRank[i] == rank )
            {
                myExportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner) );
                myImportList.emplace_back(i, rank, static_cast<char>(AttributeSet::owner), -1 );
            }
        }
        std::inplace_merge(myImportList.begin(), myImportList.begin() + numImport, myImportList.end());
        std::inplace_merge(myExportList.begin(), myExportList.begin() + numExport, myExportList.end());
    }


    std::vector<std::pair<std::string, bool>> parallel_wells;
    if( wells )
    {
        // complete root's export and other's import list by adding remaining well cells
        extendExportAndImportLists(gog, cc, root, myExportList, myImportList, gIDtoRank);

        auto wellRanks = getWellRanks(gIDtoRank, wellConnections);
        parallel_wells = wellsOnThisRank(*wells, wellRanks, cc, root);
    }
    return std::make_tuple( std::move(gIDtoRank),
                            std::move(parallel_wells),
                            std::move(myExportList),
                            std::move(myImportList) );
}

namespace {
void setDefaultZoltanParameters(Zoltan_Struct* zz)
{
    Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "PHG_EDGE_SIZE_THRESHOLD", ".35");  /* 0-remove all, 1-remove none */
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
#ifndef NDEBUG
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "2");
#else
    Zoltan_Set_Param(zz, "CHECK_GRAPH", "0");
#endif
}

} // anon namespace

std::tuple<std::vector<long long>, std::vector<std::pair<std::string, bool>>,
           std::vector<std::tuple<long long,long long,char> >,
           std::vector<std::tuple<long long,long long,char,long long> >,
           Dune::cpgrid::WellConnections>
zoltanPartitioningWithGraphOfGrid(const Dune::CpGrid& grid,
                                  const std::vector<Dune::cpgrid::OpmWellType> * wells,
                                  const std::unordered_map<std::string, std::set<long long>>& possibleFutureConnections,
                 [[maybe_unused]] const double* transmissibilities,
                                  const Dune::cpgrid::CpGridDataTraits::Communication& cc,
                 [[maybe_unused]] Dune::EdgeWeightMethod edgeWeightsMethod,
                                  long long root,
                                  const double zoltanImbalanceTol,
                                  const std::map<std::string, std::string>& params)
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

    // root process has the whole grid, other ranks nothing
    bool partitionIsEmpty = cc.rank()!=root;

    // prepare graph and contract well cells
    // non-root processes have empty grid and no wells
    GraphOfGrid gog(grid, transmissibilities);
    assert(gog.size()==0 || !partitionIsEmpty);
    auto wellConnections=partitionIsEmpty ? Dune::cpgrid::WellConnections()
                                          : Dune::cpgrid::WellConnections(*wells, possibleFutureConnections, grid);
    addWellConnections(gog, wellConnections);

    // call partitioner
    setGraphOfGridZoltanGraphFunctions(zz, gog, partitionIsEmpty);
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

    // arrange output into tuples and add well cells
    auto importExportLists = makeImportAndExportLists(gog,
                                                      cc,
                                                      wells,
                                                      wellConnections,
                                                      root,
                                                      numExport,
                                                      numImport,
                                                      exportLocalGids,
                                                      exportGlobalGids,
                                                      exportProcs,
                                                      importGlobalGids);

    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_Destroy(&zz);

    // add wellConnections to the importExportLists and return it
    auto result = std::tuple(std::move(std::get<0>(importExportLists)),
                             std::move(std::get<1>(importExportLists)),
                             std::move(std::get<2>(importExportLists)),
                             std::move(std::get<3>(importExportLists)),
                             std::move(wellConnections));
    return result;
}
#endif // HAVE_MPI

// explicit template instantiations

} // end namespace Opm
