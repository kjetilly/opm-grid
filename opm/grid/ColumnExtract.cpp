/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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
*/
#include <config.h>
#include <opm/grid/ColumnExtract.hpp>

#include <opm/grid/UnstructuredGrid.h>

#include <algorithm>
#include <map>

namespace {

/// Helper struct for extractColumn
/// Compares the underlying k-index
struct ExtractColumnCompare
{
    ExtractColumnCompare(const UnstructuredGrid& g)
    : grid(g)
    {
        // empty
    }

    bool operator()(const long long i, const long long j)
    {
        // Extract k-index
        long long index_i = grid.global_cell ? grid.global_cell[i] : i;
        long long k_i = index_i / grid.cartdims[0] / grid.cartdims[1];
        long long index_j = grid.global_cell ? grid.global_cell[j] : j;
        long long k_j = index_j / grid.cartdims[0] / grid.cartdims[1];

        return k_i < k_j;
    }

    const UnstructuredGrid& grid;
};

/// Neighbourhood query.
/// \return true if two cells are neighbours.
bool neighbours(const UnstructuredGrid& grid, const long long c0, const long long c1)
{
    for (unsigned hf = grid.cell_facepos[c0]; hf < grid.cell_facepos[c0 + 1]; ++hf) {
        const long long f = grid.cell_faces[hf];
        if (grid.face_cells[2*f] == c1 || grid.face_cells[2*f+1] == c1) {
            return true;
        }
    }
    return false;
}

} // anonymous namespace


namespace Opm {

void extractColumn(const UnstructuredGrid& grid, std::vector<std::vector<long long> >& columns)
{
    const long long* dims = grid.cartdims;

    // Keeps track of column_index ---> index of vector
    std::map<long long, long long> global_to_local;
    for (long long cell = 0; cell < grid.number_of_cells; ++cell) {
        // Extract Cartesian coordinates
        long long index = grid.global_cell ? grid.global_cell[cell] : cell; // If null, assume mapping is identity.
        long long i_cart = index % dims[0];
        long long k_cart = index / dims[0] / dims[1];
        long long j_cart = (index - k_cart*dims[0]*dims[1])/ dims[0];

        long long local_index;
        std::map<long long, long long>::iterator local_index_iterator = global_to_local.find(i_cart+j_cart*dims[0]);
        if (local_index_iterator != global_to_local.end()) {
            local_index = local_index_iterator->second;
        } else {
            local_index = columns.size();
            global_to_local[i_cart+j_cart*dims[0]] = local_index;
            columns.push_back(std::vector<long long>());
        }
        columns[local_index].push_back(cell);
    }

    long long num_cols = columns.size();
    for (long long col = 0; col < num_cols; ++col) {
        std::sort(columns[col].begin(), columns[col].end(), ExtractColumnCompare(grid));
    }

    // At this point, a column may contain multiple disjoint sets of cells.
    // We must split these columns into connected parts.
    std::vector< std::vector<long long> > new_columns;
    for (long long col = 0; col < num_cols; ++col) {
        const long long colsz = columns[col].size();
        long long first_of_col = 0;
        for (long long k = 1; k < colsz; ++k) {
            const long long c0 = columns[col][k-1];
            const long long c1 = columns[col][k];
            if (!neighbours(grid, c0, c1)) {
                // Must split. Move the cells [first_of_col, ... , k-1] to
                // a new column, known to be connected.
                new_columns.push_back(std::vector<long long>());
                new_columns.back().assign(columns[col].begin() + first_of_col, columns[col].begin() + k);
                // The working column now starts with index k.
                first_of_col = k;
            }
        }
        if (first_of_col != 0) {
            // The column was split, the working part should be
            // the entire column. We erase the cells before first_of_col.
            // (Could be more efficient if we instead chop off end.)
            columns[col].erase(columns[col].begin(), columns[col].begin() + first_of_col);
        }
    }

    // Must tack on the new columns to complete the set.
    const long long num_cols_all = num_cols + new_columns.size();
    columns.resize(num_cols_all);
    for (long long col = num_cols; col < num_cols_all; ++col) {
        columns[col].swap(new_columns[col - num_cols]);
    }
}

} // namespace Opm
