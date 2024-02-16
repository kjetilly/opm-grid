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

// Note: This is not exported from preprocess, hence we need to redefine it here
int
linearindex(const std::array<int, 3>& dims, int i, int j, int k)
{
    assert(0 <= i);
    assert(0 <= j);
    assert(0 <= k);

    assert(i < dims[0]);
    assert(j < dims[1]);
    assert(k < dims[2]);

    return i + dims[0] * (j + dims[1] * k);
}

std::vector<double>
makeSumIdirAtK(const int nx, const int ny, const int k, const std::vector<double>& dx)
{
    std::vector<double> s(nx * ny, 0.0);
    for (int j = 0; j < ny; ++j) {
        double sum = 0.0;
        for (int i = 0; i < nx; ++i) {
            sum += dx[i + j * nx + k * nx * ny];
            s[i + j * nx] = sum;
        }
    }
    return s;
}

std::vector<double>
makeSumJdirAtK(const int nx, const int ny, const int k, const std::vector<double>& dy)
{
    std::vector<double> s(nx * ny, 0.0);
    for (int i = 0; i < nx; ++i) {
        double sum = 0.0;
        for (int j = 0; j < ny; ++j) {
            sum += dy[i + j * nx + k * nx * ny];
            s[i + j * nx] = sum;
        }
    }
    return s;
}

std::vector<double>
makeSumKdir(const int nx, const int ny, const int nz, const std::vector<double>& dz)
{
    std::vector<double> s(nx * ny, 0.0);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nz; ++k) {
                sum += dz[i + j * nx + k * nx * ny];
            }
            s[i + j * nx] = sum;
        }
    }
    return s;
}

std::vector<double>
makeCoordDxDyDzTops(const std::array<int, 3>& dims,
                    const std::vector<double>& dx,
                    const std::vector<double>& dy,
                    const std::vector<double>& dz,
                    const std::vector<double>& tops)
{
    const size_t nx = size_t(dims[0]);
    const size_t ny = size_t(dims[1]);
    const size_t nz = size_t(dims[2]);



    std::vector<double> coord;
    coord.reserve((nx + 1) * (ny + 1) * 6);

    std::vector<double> sum_idir_top = makeSumIdirAtK(nx, ny, 0, dx);
    std::vector<double> sum_idir_bot = makeSumIdirAtK(nx, ny, nz - 1, dx);
    std::vector<double> sum_jdir_top = makeSumJdirAtK(nx, ny, 0, dy);
    std::vector<double> sum_jdir_bot = makeSumJdirAtK(nx, ny, nz - 1, dy);
    std::vector<double> sum_kdir = makeSumKdir(nx, ny, nz, dz);

    for (std::size_t j = 0; j < ny; j++) {

        double y0 = 0;
        double zt = tops[0];
        double zb = zt + sum_kdir[0 + 0 * nx];

        if (j == 0) {
            double x0 = 0.0;

            coord.push_back(x0);
            coord.push_back(y0);
            coord.push_back(zt);
            coord.push_back(x0);
            coord.push_back(y0);
            coord.push_back(zb);

            for (std::size_t i = 0; i < nx; i++) {

                size_t ind = i + j * nx + 1;

                if (i == (nx - 1)) {
                    ind = ind - 1;
                }

                zt = tops[ind];
                zb = zt + sum_kdir[i + j * nx];

                double xt = x0 + dx[i + j * nx];
                double xb = sum_idir_bot[i + j * nx];

                coord.push_back(xt);
                coord.push_back(y0);
                coord.push_back(zt);
                coord.push_back(xb);
                coord.push_back(y0);
                coord.push_back(zb);

                x0 = xt;
            }
        }

        std::size_t ind = (j + 1) * nx;

        if (j == (ny - 1)) {
            ind = j * nx;
        }

        double x0 = 0.0;

        double yt = sum_jdir_top[0 + j * nx];
        double yb = sum_jdir_bot[0 + j * nx];

        zt = tops[ind];
        zb = zt + sum_kdir[0 + j * nx];

        coord.push_back(x0);
        coord.push_back(yt);
        coord.push_back(zt);
        coord.push_back(x0);
        coord.push_back(yb);
        coord.push_back(zb);

        for (std::size_t i = 0; i < nx; i++) {

            ind = i + (j + 1) * nx + 1;

            if (j == (ny - 1)) {
                ind = i + j * nx + 1;
            }

            if (i == (nx - 1)) {
                ind = ind - 1;
            }

            zt = tops[ind];
            zb = zt + sum_kdir[i + j * nx];

            double xt = -999;
            double xb;

            if (j == (ny - 1)) {
                xt = sum_idir_top[i + j * nx];
                xb = sum_idir_bot[i + j * nx];
            } else {
                xt = sum_idir_top[i + (j + 1) * nx];
                xb = sum_idir_bot[i + (j + 1) * nx];
            }

            if (i == (nx - 1)) {
                yt = sum_jdir_top[i + j * nx];
                yb = sum_jdir_bot[i + j * nx];
            } else {
                yt = sum_jdir_top[(i + 1) + j * nx];
                yb = sum_jdir_bot[(i + 1) + j * nx];
            }

            coord.push_back(xt);
            coord.push_back(yt);
            coord.push_back(zt);
            coord.push_back(xb);
            coord.push_back(yb);
            coord.push_back(zb);
        }
    }

    return coord;
}

std::vector<double>
makeZcornDzTops(const std::array<int, 3>& dims, const std::vector<double>& dz, const std::vector<double>& tops)
{

    std::vector<double> zcorn;

    const size_t nx = size_t(dims[0]);
    const size_t ny = size_t(dims[1]);
    const size_t nz = size_t(dims[2]);

    size_t sizeZcorn = nx * ny * nz * 8;

    zcorn.assign(sizeZcorn, 0.0);

    for (std::size_t j = 0; j < ny; j++) {
        for (std::size_t i = 0; i < nx; i++) {
            std::size_t ind = i + j * nx;
            double z = tops[ind];

            for (std::size_t k = 0; k < nz; k++) {

                // top face of cell
                std::size_t zind = i * 2 + j * nx * 4 + k * nx * ny * 8;

                zcorn[zind] = z;
                zcorn[zind + 1] = z;

                zind = zind + nx * 2;

                zcorn[zind] = z;
                zcorn[zind + 1] = z;

                z = z + dz[i + j * nx + k * nx * ny];

                // bottom face of cell
                zind = i * 2 + j * nx * 4 + k * nx * ny * 8 + nx * ny * 4;

                zcorn[zind] = z;
                zcorn[zind + 1] = z;

                zind = zind + nx * 2;

                zcorn[zind] = z;
                zcorn[zind + 1] = z;
            }
        }
    }

    return zcorn;
}
BOOST_AUTO_TEST_CASE(TestLargeGrid)
{
    const int nx = 512;
    const int ny = 512;
    const int nz = 512;

    std::vector<double> dx(nx * ny * nz, 1.0);
    std::vector<double> dy(nx * ny * nz, 1.0);
    std::vector<double> dz(nx * ny * nz, 1.0);
    std::vector<double> tops(nx * ny, nz * 1.0);

    const std::array<int, 3> dims {{nz, ny, nz}};

    std::vector<double> coord = makeCoordDxDyDzTops(dims, dx, dy, dz, tops); // ((nx + 1) * (ny + 1) * 6, 1.0);
    std::vector<double> zcorn = makeZcornDzTops(dims, dz, tops); //(8 * nx * ny * nz, 0.0);
    std::vector<int> actnum(nx * ny * nz, 1);
    std::vector<int> isAquaCell(nx * ny * nz, 0);


    grdecl g {.dims = {nx, ny, nz}, .coord = coord.data(), .zcorn = zcorn.data(), .actnum = actnum.data()};

    processed_grid processedGrid;
    int pinchActive = 0;

    process_grdecl(&g, 1e-16, isAquaCell.data(), &processedGrid, pinchActive);
}
