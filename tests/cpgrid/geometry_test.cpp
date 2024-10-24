/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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
#include "config.h"

#define BOOST_TEST_MODULE GeometryTests
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <opm/grid/cpgrid/Geometry.hpp>
#include <opm/grid/cpgrid/EntityRep.hpp>

#include <sstream>
#include <iostream>

using namespace Dune;

class Null;

BOOST_AUTO_TEST_CASE(vertexgeom)
{
    typedef cpgrid::Geometry<0, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point.
    Geometry::GlobalCoordinate c(3.0);
    Geometry g(c);

    // Verification of properties.
    BOOST_CHECK(g.type().isVertex());
    BOOST_CHECK(g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 1);
    BOOST_CHECK_EQUAL(g.corner(0), c);
    Geometry::LocalCoordinate lc(0.0);
    BOOST_CHECK_EQUAL(g.global(lc), c);
    // BOOST_CHECK_THROW(g.local(c), std::exception);
    BOOST_CHECK_EQUAL(g.integrationElement(lc), 1.0);
    BOOST_CHECK_EQUAL(g.volume(), 1.0);
    BOOST_CHECK_EQUAL(g.center(), c);
    // BOOST_CHECK_THROW(g.jacobianTransposed(lc), std::exception);
    // BOOST_CHECK_THROW(g.jacobianInverseTransposed(lc), std::exception);
}


BOOST_AUTO_TEST_CASE(intersectiongeom)
{
    typedef cpgrid::Geometry<2, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point and volume.
    Geometry::GlobalCoordinate c(3.0);
    Geometry::ctype v = 8.0;
    Geometry g(c, v);

    // Verification of properties.
    BOOST_CHECK(g.type().isNone());
    BOOST_CHECK(g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 0);
    // BOOST_CHECK_THROW(g.corner(0), std::exception);
    Geometry::LocalCoordinate lc(0.0);
    BOOST_CHECK_THROW(g.global(lc), std::exception);
    BOOST_CHECK_THROW(g.local(c), std::exception);
    BOOST_CHECK_EQUAL(g.integrationElement(lc), v);
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);
    BOOST_CHECK_THROW(g.jacobianTransposed(lc), std::exception);
    BOOST_CHECK_THROW(g.jacobianInverseTransposed(lc), std::exception);
}


BOOST_AUTO_TEST_CASE(cellgeom)
{
    typedef cpgrid::Geometry<3, 3> Geometry;
    // Default construction.
    Geometry g_default;

    // Construction from point and volume.
    // This is a dangerous constructor kept for backwards compatibility,
    // these checks may be removed if constructor removed.
    Geometry::GlobalCoordinate c(3.0);
    Geometry::ctype v = 8.0;
    Geometry g_dangerous(c, v);

    // Construction from point, volume, points and pointindices.
    // First a unit cube, i.e. the mapping represented is the identity.
    typedef Geometry::GlobalCoordinate GC;
    c = GC(0.5);
    v = 1.0;
    GC corners[8];
    GC cor;
    for (int k = 0; k < 2; ++k) {
        cor[2] = k;
        for (int j = 0; j < 2; ++j) {
            cor[1] = j;
            for (int i = 0; i < 2; ++i) {
                cor[0] = i;
                corners[4*k + 2*j + i] = cor;
            }
        }
    }
//     for (int i = 0; i < 8; ++i) {
//         std::cout << corners[i] << std::endl;
//     }
    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg;
    pg.reserve(8);
    for (const auto& crn : corners)
    {
        pg.push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    Geometry g(c, v, pg, cor_idx);

    // Verification of properties.
    BOOST_CHECK(g.type().isCube());
    BOOST_CHECK(!g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 8);
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(g.corner(i), corners[i]);
    }
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);

    // Verification of properties that depend on the mapping.
    typedef Geometry::LocalCoordinate LC;
    const int N = 5;
    const int num_pts = N*N*N;
    LC testpts[num_pts] = { LC(0.0) };
    LC pt(0.0);
    for (int i = 0; i < N; ++i) {
        pt[0] = double(i)/double(N-1);
        for (int j = 0; j < N; ++j) {
            pt[1] = double(j)/double(N-1);
            for (int k = 0; k < N; ++k) {
                pt[2] = double(k)/double(N-1);
                testpts[i*N*N + j*N + k] = pt;
//                 std::cout << pt << std::endl;
            }
        }
    }
    Geometry::JacobianTransposed id(0.0);
    id[0][0] = id[1][1] = id[2][2] = 1.0;
    for (int i = 0; i < num_pts; ++i) {
        BOOST_CHECK_EQUAL(g.global(testpts[i]), testpts[i]);
        BOOST_CHECK_EQUAL(g.local(g.global(testpts[i])), testpts[i]);
        BOOST_CHECK_EQUAL(g.integrationElement(testpts[i]), 1.0);
        BOOST_CHECK_EQUAL(g.jacobianTransposed(testpts[i]), id);
        BOOST_CHECK_EQUAL(g.jacobianInverseTransposed(testpts[i]), id);
    }

    // Next testcase: a degenerate hexahedron, wedge shaped.
    typedef Geometry::GlobalCoordinate GC;
    c = GC(1.0/3.0); c[2] = 0.5;
    v = 0.5;
    corners[5][2] = 0.0;
    corners[7][2] = 0.0;
    
    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg1;
    pg1.reserve(8);
    for (const auto& crn : corners)
    {
        pg1.push_back(cpgrid::Geometry<0, 3>(crn));
    }
    g = Geometry(c, v, pg1, cor_idx);

    // Verification of properties.
    BOOST_CHECK(g.type().isCube());
    BOOST_CHECK(!g.affine());
    BOOST_CHECK_EQUAL(g.corners(), 8);
    for (int i = 0; i < 8; ++i) {
        BOOST_CHECK_EQUAL(g.corner(i), corners[i]);
    }
    BOOST_CHECK_EQUAL(g.volume(), v);
    BOOST_CHECK_EQUAL(g.center(), c);

    struct Wedge
    {
        static GC global(const LC& lc)
        {
            GC gc(0.0);
            gc[0] = lc[0];
            gc[1] = lc[1];
            gc[2] = (1.0 - lc[0])*lc[2];
            return gc;
        }
        static double integrationElement(const LC& lc)
        {
            return 1.0 - lc[0];
        }
        static Geometry::JacobianTransposed jacobianTransposed(const LC& lc)
        {
            Geometry::JacobianTransposed Jt(0.0);
            Jt[0][0] = 1.0;
            Jt[0][2] = -lc[2];
            Jt[1][1] = 1.0;
            Jt[2][2] = 1.0 - lc[0];
            return Jt;
        }
    };

    // Verification of properties that depend on the mapping.
    const double tolerance = 1e-14;
    for (int i = 0; i < num_pts; ++i) {
        GC gl = g.global(testpts[i]);
        BOOST_CHECK_EQUAL(gl, Wedge::global(testpts[i]));
        BOOST_CHECK_EQUAL(g.integrationElement(testpts[i]), Wedge::integrationElement(testpts[i]));
        Geometry::JacobianTransposed Jt = Wedge::jacobianTransposed(testpts[i]);
        BOOST_CHECK_EQUAL(g.jacobianTransposed(testpts[i]), Jt);
        if (testpts[i][0] < 1.0) {
            // Only do this test if we are away from the degeneracy.
            LC diff = g.local(gl);
            diff -= testpts[i];
            BOOST_CHECK_SMALL(diff.two_norm(), tolerance);
            Geometry::Jacobian Jit = Jt; // This implicitly assumes that the Jacobian is square.
            Jit.invert();
            BOOST_CHECK_EQUAL(g.jacobianInverseTransposed(testpts[i]), Jit);
        }
    }
}


template <typename T>
inline void
check_coordinates(T c1, T c2)
{
    for (int c = 0; c < 3; c++) {
        BOOST_CHECK_CLOSE(c1[c], c2[c], 1e-6);
    }
}

void
check_refined_grid(const cpgrid::Geometry<3, 3>& parent,
                   const std::vector<cpgrid::Geometry<3, 3>>& refined,
                   const std::array<int, 3>& cells_per_dim)
{
    using Geometry = cpgrid::Geometry<3, 3>;
    using GlobalCoordinate = Geometry::GlobalCoordinate;

    int count = cells_per_dim[0] * cells_per_dim[1] * cells_per_dim[2];
    BOOST_CHECK_EQUAL(refined.size(), count);

    // Check the corners of the refined grid with the parent corners.
    int idx = 0;
    for (int k = 0; k < 1; k++) {
        int slice = cells_per_dim[0] * cells_per_dim[1];
        for (int j = 0; j < 1; j++) {
            for (int i = 0; i < 1; i++) {
                auto& r = refined[k * slice + j * cells_per_dim[0] + i];
                check_coordinates(r.corner(idx), parent.corner(idx));
                idx++;
            }
        }
    }

    // Make sure the corners of neighboring cells overlap.
    for (int k = 0; k < cells_per_dim[2]; k++) {
        int slice = cells_per_dim[1] * cells_per_dim[0];
        for (int j = 0; j < cells_per_dim[1]; j++) {
            for (int i = 0; i < cells_per_dim[0]; i++) {
                auto& r0 = refined[k * slice + cells_per_dim[0] * j + i];
                if (i < cells_per_dim[0] - 1) {
                    auto& r1 = refined[k * slice + cells_per_dim[0] * j + i + 1];
                    check_coordinates(r0.corner(1), r1.corner(0));
                    check_coordinates(r0.corner(3), r1.corner(2));
                    check_coordinates(r0.corner(5), r1.corner(4));
                    check_coordinates(r0.corner(7), r1.corner(6));
                }
                if (j < cells_per_dim[1] - 1) {
                    auto& r1 = refined[k * slice + cells_per_dim[0] * (j + 1) + i];
                    check_coordinates(r0.corner(2), r1.corner(0));
                    check_coordinates(r0.corner(3), r1.corner(1));
                    check_coordinates(r0.corner(6), r1.corner(4));
                    check_coordinates(r0.corner(7), r1.corner(5));
                }
                if (k < cells_per_dim[2] - 1) {
                    auto& r1 = refined[(k + 1) * slice + cells_per_dim[0] * j + i];
                    check_coordinates(r0.corner(4), r1.corner(0));
                    check_coordinates(r0.corner(5), r1.corner(1));
                    check_coordinates(r0.corner(6), r1.corner(2));
                    check_coordinates(r0.corner(7), r1.corner(3));
                }
            }
        }
    }

    // Check the centers of the cells.
    for (auto r : refined) {
        GlobalCoordinate center = {0.0, 0.0, 0.0};
        for (int h = 0; h < 8; h++) {
            for (int c = 0; c < 3; c++) {
                center[c] += r.corner(h)[c] / 8.0;
            }
        }
        check_coordinates(r.center(), center);
    }

    // Check that the weighted mean of all centers equals the parent center
    GlobalCoordinate center = {0.0, 0.0, 0.0};
    for (auto r : refined) {
        for (int c = 0; c < 3; c++) {
            center[c] += r.center()[c] * r.volume() / parent.volume();
        }
    }
    check_coordinates(parent.center(), center);

    // Check that mean of all corners equals the center of the parent.
    center = {0.0, 0.0, 0.0};
    for (auto r : refined) {
        for (int h = 0; h < 8; h++) {
            for (int c = 0; c < 3; c++) {
                center[c] += r.corner(h)[c] / count / 8;
            }
        }
    }
    check_coordinates(parent.center(), center);

    // Check the total volume against the parent volume.
    Geometry::ctype volume = 0.0;
    for (auto r : refined) {
        volume += r.volume();
    }
    BOOST_CHECK_CLOSE(volume, parent.volume(), 1e-6);
}


BOOST_AUTO_TEST_CASE(refine_simple_cube)
{
    using Geometry = cpgrid::Geometry<3, 3>;
    using GlobalCoordinate = Geometry::GlobalCoordinate;

    const GlobalCoordinate corners[8] = {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {0.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
    };
    const GlobalCoordinate c = {0.5, 0.5, 0.5};
    const Geometry::ctype v = 1.0;

    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg;
    for (const auto& crn : corners) {
        pg.push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    Geometry g(c, v, pg, cor_idx);

    {
        std::array<int, 3> cells = {1, 1, 1};
        std::vector<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> gi(1);
        std::vector<std::array<int, 8>> ci(1);
        std::vector<Geometry> refined = g.refine(cells, gi, ci);
        check_refined_grid(g, refined, cells);
    }

    {
        std::array<int, 3> cells = {2, 3, 4};
        int cell_count = cells[0] * cells[1] * cells[2];
        std::vector<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> gi(cell_count);
        std::vector<std::array<int, 8>> ci(cell_count);
        std::vector<Geometry> refined = g.refine(cells, gi, ci);
        check_refined_grid(g, refined, cells);
    }
}


BOOST_AUTO_TEST_CASE(refine_distorted_cube)
{
    using Geometry = cpgrid::Geometry<3, 3>;
    using GlobalCoordinate = Geometry::GlobalCoordinate;

    // Distorted cube:
    const GlobalCoordinate corners[8] = {
        {0.1, 0.2, 0.3},
        {1.2, 0.3, 0.4},
        {0.3, 1.4, 0.5},
        {1.4, 1.5, 0.6},
        {0.5, 0.6, 1.7},
        {1.6, 0.7, 1.8},
        {0.7, 1.8, 1.9},
        {1.8, 1.9, 2.0},
    };

    // Arbitrary volume:
    const Geometry::ctype v = 123.0;

    // Calculate the centroid:
    GlobalCoordinate center = {0.0, 0.0, 0.0};
    for (int h = 0; h < 8; h++) {
        for (int c = 0; c < 3; c++) {
            center[c] += corners[h][c] / 8.0;
        }
    }

    cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3> pg;
    for (const auto& crn : corners) {
        pg.push_back(cpgrid::Geometry<0, 3>(crn));
    }

    int cor_idx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    Geometry g(center, v, pg, cor_idx);

    {
        std::array<int, 3> cells = {1, 1, 1};
        std::vector<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> gi(1);
        std::vector<std::array<int, 8>> ci(1);
        std::vector<Geometry> refined = g.refine(cells, gi, ci);
        check_refined_grid(g, refined, cells);
    }

    {
        std::array<int, 3> cells = {2, 3, 4};
        int cell_count = cells[0] * cells[1] * cells[2];
        std::vector<cpgrid::EntityVariable<cpgrid::Geometry<0, 3>, 3>> gi(cell_count);
        std::vector<std::array<int, 8>> ci(cell_count);
        std::vector<Geometry> refined = g.refine(cells, gi, ci);
        check_refined_grid(g, refined, cells);
    }
}
