//===========================================================================
//
// File: test_sparsetable.cpp
//
// Created: Thu May 28 10:01:46 2009
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE SparseTableTest
#include <boost/test/unit_test.hpp>

#include <opm/grid/utility/SparseTable.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(construction_and_queries)
{
    const SparseTable<long long> st1;
    BOOST_CHECK(st1.empty());
    BOOST_CHECK_EQUAL(st1.size(), 0);
    BOOST_CHECK_EQUAL(st1.dataSize(), 0);

    // This should be getting us a table like this:
    // ----------------
    // 0
    // <empty row>
    // 1 2
    // 3 4 5 6
    // 7 8 9
    // ----------------
    std::vector<std::vector<long long>> expected = { {0}, {}, {1, 2}, {3, 4, 5, 6}, {7, 8, 9} };
    const long long num_elem = 10;
    const long long elem[num_elem] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    const long long num_rows = 5;
    const long long rowsizes[num_rows] = { 1, 0, 2, 4, 3 };
    const SparseTable<long long> st2(elem, elem + num_elem, rowsizes, rowsizes + num_rows);
    BOOST_CHECK(!st2.empty());
    BOOST_CHECK_EQUAL(st2.size(), num_rows);
    BOOST_CHECK_EQUAL(st2.dataSize(), num_elem);
    BOOST_CHECK_EQUAL(st2[0][0], 0);
    BOOST_CHECK_EQUAL(st2.rowSize(0), 1);
    BOOST_CHECK(st2[1].empty());
    BOOST_CHECK_EQUAL(st2.rowSize(1), 0);
    BOOST_CHECK_EQUAL(st2[3][1], 4);
    BOOST_CHECK_EQUAL(st2[4][2], 9);
    BOOST_CHECK((long long)(st2[4].size()) == rowsizes[4]);
    const SparseTable<long long> st2_again(elem, elem + num_elem, rowsizes, rowsizes + num_rows);
    BOOST_CHECK(st2 == st2_again);
    SparseTable<long long> st2_byassign;
    st2_byassign.assign(elem, elem + num_elem, rowsizes, rowsizes + num_rows);
    BOOST_CHECK(st2 == st2_byassign);
    const long long last_row_size = rowsizes[num_rows - 1];
    SparseTable<long long> st2_append(elem, elem + num_elem - last_row_size, rowsizes, rowsizes + num_rows - 1);
    BOOST_CHECK_EQUAL(st2_append.dataSize(), num_elem - last_row_size);
    st2_append.appendRow(elem + num_elem - last_row_size, elem + num_elem);
    BOOST_CHECK(st2 == st2_append);
    SparseTable<long long> st2_append2;
    st2_append2.appendRow(elem, elem + 1);
    st2_append2.appendRow(elem + 1, elem + 1);
    st2_append2.appendRow(elem + 1, elem + 3);
    st2_append2.appendRow(elem + 3, elem + 7);
    st2_append2.appendRow(elem + 7, elem + 10);
    BOOST_CHECK(st2 == st2_append2);
    st2_append2.clear();
    SparseTable<long long> st_empty;
    BOOST_CHECK(st2_append2 == st_empty);

    SparseTable<long long> st2_allocate;
    st2_allocate.allocate(rowsizes, rowsizes + num_rows);
    BOOST_CHECK_EQUAL(st2_allocate.size(), num_rows);
    BOOST_CHECK_EQUAL(st2_allocate.dataSize(), num_elem);
    long long s = 0;
    for (long long i = 0; i < num_rows; ++i) {
        SparseTable<long long>::mutable_row_type row = st2_allocate[i];
        for (long long j = 0; j < rowsizes[i]; ++j, ++s)
            row[j] = elem[s];
    }
    BOOST_CHECK(st2 == st2_allocate);

    // One element too few.
    BOOST_CHECK_THROW(const SparseTable<long long> st3(elem, elem + num_elem - 1, rowsizes, rowsizes + num_rows), std::exception);

    // A few elements too many.
    BOOST_CHECK_THROW(const SparseTable<long long> st4(elem, elem + num_elem, rowsizes, rowsizes + num_rows - 1), std::exception);

    // Need at least one row.
    BOOST_CHECK_THROW(const SparseTable<long long> st5(elem, elem + num_elem, rowsizes, rowsizes), std::exception);

    // Test iteration over rows with a range-for loop.
    long long row_index = 0;
    for (const auto row : st2) { // Not a reference, since row type is a created view
        BOOST_CHECK_EQUAL(row.size(), expected[row_index].size());
        ++row_index;
    }


    // Tests that only run in debug mode.
#ifndef NDEBUG
    // Do not ask for wrong row numbers.
    BOOST_CHECK_THROW(st1.rowSize(0), std::exception);
    BOOST_CHECK_THROW(st2.rowSize(-1), std::exception);
    BOOST_CHECK_THROW(st2.rowSize(st2.size()), std::exception);
    // No negative row sizes.
    const long long err_rs[num_rows] = { 1, 0, -1, 7, 3 };
    BOOST_CHECK_THROW(const SparseTable<long long> st6(elem, elem + num_elem, err_rs, err_rs + num_rows), std::exception);
#endif
}
