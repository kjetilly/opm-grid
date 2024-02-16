#include <config.h>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <dune/common/unused.hh>
#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <opm/grid/cpgrid/dgfparser.hh>


#define DISABLE_DEPRECATED_METHOD_CHECK 1
using Dune::referenceElement; // grid check assume usage of Dune::Geometry
#include <dune/grid/test/gridcheck.hh>


// Re-enable warnings.
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <iostream>


#include <opm/grid/cpgpreprocess/preprocess.h>

#define BOOST_TEST_MODULE CPGridIntOverflowTest
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(TestLargeGrid)
{
    const int nx = 128;
    const int ny = 128;
    const int nz = 128;

    std::vector<double> coord((nx + 1) * (ny + 1) * 6, 1.0);
    std::vector<double> zcorn(8 * nx * ny * nz, 0.0);
    std::vector<int> actnum(nx * ny * nz, 1);
    std::vector<int> isAquaCell(nx * ny * nz, 0);
    const std::array<int, 3> dims {{nz, ny, nz}};
    grdecl g {.dims = {nx, ny, nz}, .coord = coord.data(), .zcorn = zcorn.data(), .actnum = actnum.data()};

    processed_grid processedGrid;
    int pinchActive = 0;

    process_grdecl(&g, 1e-16, isAquaCell.data(), &processedGrid, pinchActive);
}
