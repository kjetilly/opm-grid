// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GRIDFACTORY_HH
#define DUNE_POLYHEDRALGRID_GRIDFACTORY_HH

#include <algorithm>
#include <numeric>

#include <dune/common/typetraits.hh>
#include <dune/common/version.hh>

#include <dune/grid/common/gridfactory.hh>
#include <opm/grid/polyhedralgrid/grid.hh>

namespace Dune
{


  // GridFactory for PolyhedralGrid
  // ---------------------------------

  template< long long dim, long long dimworld, class coord_t >
  class GridFactory< PolyhedralGrid< dim, dimworld, coord_t > >
    : public GridFactoryInterface< PolyhedralGrid< dim, dimworld, coord_t > >
  {
  public:
    typedef PolyhedralGrid< dim, dimworld, coord_t > Grid;

    const static long long dimension      = Grid::dimension;
    const static long long dimensionworld = Grid::dimensionworld;
    typedef typename Grid::ctype ctype;

    typedef MPIHelper::MPICommunicator MPICommunicatorType;
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<dimension>::Entity Vertex;

    typedef Dune::FieldVector<ctype,dimensionworld> CoordinateType;
    typedef CoordinateType  Coordinate;

    using UniquePtrType = std::unique_ptr<Grid>;

    /** \brief Default constructor */
    explicit GridFactory ( const MPICommunicatorType& = MPIHelper::getCommunicator() )
      : nodes_(),
        faces_(),
        cells_()
    {}

    virtual void insertVertex(const CoordinateType& pos)
    {
      nodes_.push_back( pos );
    }

    /** \brief Insert an element into the coarse grid
        \param type  The GeometryType of the new element
        \param items The items are usually the vertex numbers of the inserted
               element. If the geometry type is none the these can be face numbers.

        \note If the GeometryType is none then faces need to be inserted separately
              using this method and passing a GeometryType with dimension - 1
              (with respect to the Grid's dimension).
     */
    virtual void insertElement(const GeometryType& type,
                               const std::vector<size_t>& items)
    {
      if( type.isNone() )
      {
        // copy into vector of integers
        std::vector< long long > numbers( items.size() );
        std::copy( items.begin(), items.end(), numbers.begin() );

        if( type.dim() == dimension-1 )
        {
          faces_.push_back( numbers );
        }
        else if( type.dim() == dimension )
        {
          // note vertices holds the face
          // numbers in this case
          cells_.push_back( numbers );
        }
        else
        {
          DUNE_THROW(Dune::NotImplemented,"insertElement not implemented for type " << type );
        }
      }
      else // use ReferenceElement to insert faces
      {


      }
    }

    void insertBoundarySegment(const std::vector<size_t>&)
    {
      DUNE_THROW(NotImplemented,"yet");
    }

    UniquePtrType createGrid()
    {
      std::vector< CoordinateType >& nodes = nodes_;
      std::vector< std::vector< long long > >& faces = faces_;
      std::vector< std::vector< long long > >& cells = cells_;

      if( cells.empty() )
      {
        DUNE_THROW( GridError, "No cells found for PolyhedralGrid" );
      }

      const auto sumSize = [] ( std::size_t s, const std::vector< long long > &v ) { return s + v.size(); };
      const std::size_t numFaceNodes = std::accumulate( faces.begin(), faces.end(), std::size_t( 0 ), sumSize );
      const std::size_t numCellFaces = std::accumulate( cells.begin(), cells.end(), std::size_t( 0 ), sumSize );

      typename Grid::UnstructuredGridPtr ug =
        Grid::allocateGrid( cells.size(), faces.size(), numFaceNodes, numCellFaces, nodes.size() );

      // copy faces
      {
#ifndef NDEBUG
        std::map< std::vector< long long >, std::vector< long long > > faceMap;
#endif

        const long long nFaces = faces.size();
        // set all face_cells values to -2 as default
        std::fill( ug->face_cells, ug->face_cells + 2*nFaces, -1 );

        long long facepos = 0;
        std::vector< long long > faceVertices;
        faceVertices.reserve( 30 );
        for( long long face = 0; face < nFaces; ++face )
        {
          //std::cout << "face " << face << ": ";
          faceVertices.clear();
          ug->face_nodepos[ face ] = facepos;
          const long long nVertices = faces[ face ].size();
          for( long long vx = 0; vx < nVertices; ++vx, ++facepos )
          {
            //std::cout << " " << faces[ face ][ vx ];
            ug->face_nodes[ facepos ] = faces[ face ][ vx ];
            faceVertices.push_back( faces[ face ][ vx ] );
          }
          //std::cout << std::endl;

#ifndef NDEBUG
          // sort vertices
          std::sort( faceVertices.begin(), faceVertices.end() );
          // make sure each face only exists once
          faceMap[ faceVertices ].push_back( face );
          assert( faceMap[ faceVertices ].size() == 1 );
#endif
        }
        ug->face_nodepos[ nFaces ] = facepos ;
      }

      // copy cells
      {
        const long long nCells = cells.size();
        long long cellpos = 0;
        for( long long cell = 0; cell < nCells; ++cell )
        {
          //std::cout << "Cell " << cell << ": ";
          ug->cell_facepos[ cell ] = cellpos;
          const long long nFaces = cells[ cell ].size();
          for( long long f = 0; f < nFaces; ++f, ++cellpos )
          {
            const long long face = cells[ cell ][ f ];
            // std::cout << " " << face ;
            ug->cell_faces[ cellpos ] = face;

            // TODO find cells for each face
            if( ug->face_cells[ 2*face ] == -1 )
            {
              ug->face_cells[ 2*face ] = cell;
            }
            else // if ( ug->face_cells[ 2*face+1 ] == -1 )
            {
              //assert( ug->face_cells[ 2*face+1 ] == -1 );
              ug->face_cells[ 2*face+1 ] = cell;
            }
          }
          //std::cout << std::endl;
        }
        ug->cell_facepos[ nCells ] = cellpos ;
      }

      // copy node coordinates
      {
        const long long nNodes = nodes.size();
        long long nodepos = 0;
        for( long long vx = 0 ; vx < nNodes; ++vx )
        {
          for( long long d=0; d<dim; ++d, ++nodepos )
            ug->node_coordinates[ nodepos ] = nodes[ vx ][ d ];
        }
      }

      /*
      for( long long i=0; i<(long long)(faces.size() ); ++i)
      {
        std::cout << "face "<< i<< " connects to " << ug->face_cells[ 2*i ] << " " <<
          ug->face_cells[ 2*i+1] << std::endl;
      }
      */

      // free cell face tag since it's not a cartesian grid
      if( ug->cell_facetag )
      {
        std::free( ug->cell_facetag );
        ug->cell_facetag = nullptr ;
        for( long long i=0; i<3; ++i ) ug->cartdims[ i ] = 0;
      }

      // compute geometric quantities like cell volume and face normals
      Grid::computeGeometry( ug );

      // check normal direction
      {
        for( long long face = 0 ; face < ug->number_of_faces; ++face )
        {
          const long long a = ug->face_cells[ 2*face     ];
          const long long b = ug->face_cells[ 2*face + 1 ];
          if( a < 0 || b < 0 )
            continue ;

          Coordinate centerDiff( 0 );
          Coordinate normal( 0 );
          //std::cout << "Cell center " << a << " " << b << std::endl;
          for( long long d=0; d<dim; ++d )
          {
            //std::cout << ug->cell_centroids[ a*dim + d ] << " " << ug->cell_centroids[ b*dim + d ] << std::endl;
            centerDiff[ d ] = ug->cell_centroids[ b*dim + d ] - ug->cell_centroids[ a*dim + d ];
            normal[ d ] = ug->face_normals[ face*dim + d ];
          }

          // if diff and normal point in different direction, flip faces
          if( centerDiff * normal > 0 )
          {
            ug->face_cells[ 2*face     ] = b;
            ug->face_cells[ 2*face + 1 ] = a;
          }
        }
      }

      return UniquePtrType( new Grid( std::move( ug ) ));
    }

  protected:
    std::vector< CoordinateType > nodes_;
    std::vector< std::vector< long long > > faces_;
    std::vector< std::vector< long long > > cells_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_DGFPARSER_HH
