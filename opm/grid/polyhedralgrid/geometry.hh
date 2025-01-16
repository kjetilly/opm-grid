// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
#define DUNE_POLYHEDRALGRID_GEOMETRY_HH

#include <memory>

#include <dune/common/fmatrix.hh>
#include <dune/grid/common/geometry.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>


namespace Dune
{

  // Internal Forward Declarations
  // -----------------------------

  template< long long, long long, class > class PolyhedralGridGeometry;
  template< long long, long long, class > class PolyhedralGridLocalGeometry;

  // PolyhedralGridBasicGeometry
  // -------------------

  template< long long mydim, long long cdim, class Grid >
  struct PolyhedralGridBasicGeometry
  {
    static const long long dimension = Grid::dimension;
    static const long long mydimension = mydim;
    static const long long codimension = dimension - mydimension;

    static const long long dimensionworld = Grid::dimensionworld;
    static const long long coorddimension = dimensionworld;

    typedef typename Grid::ctype ctype;
    typedef Dune::FieldVector< ctype, coorddimension > GlobalCoordinate;
    typedef Dune::FieldVector< ctype, mydimension >    LocalCoordinate;

    typedef typename Grid::Traits::ExtraData  ExtraData;
    typedef typename Grid::Traits::template Codim<codimension>::EntitySeed EntitySeed;

    template< class ct >
    struct PolyhedralMultiLinearGeometryTraits
      : public Dune::MultiLinearGeometryTraits< ct >
    {
      struct Storage
      {
        struct Iterator
          : public Dune::ForwardIteratorFacade< Iterator, GlobalCoordinate, GlobalCoordinate >
        {
          const Storage* data_;
          long long count_;
          explicit Iterator( const Storage* ptr, long long count ) : data_( ptr ), count_( count ) {}

          GlobalCoordinate dereference() const { return data_->corner( count_ ); }
          void increment() { ++count_; }

          bool equals( const Iterator& other ) const { return count_ == other.count_; }
        };

        ExtraData  data_;
        // host geometry object
        EntitySeed seed_;

        Storage( ExtraData data, EntitySeed seed )
          : data_( data ), seed_( seed )
        {}

        Storage( ExtraData data )
          : data_( data ), seed_()
        {}

        ExtraData data() const { return data_; }
        bool isValid () const { return seed_.isValid(); }

        GlobalCoordinate operator [] (const long long i) const { return corner( i ); }

        Iterator begin() const { return Iterator(this, 0); }
        Iterator end ()  const { return Iterator(this, corners()); }

        long long corners () const { return data()->corners( seed_ ); }
        GlobalCoordinate corner ( const long long i ) const { return data()->corner( seed_, i ); }
        GlobalCoordinate center () const { return data()->centroids( seed_ ); }

        ctype volume() const { return data()->volumes( seed_ ); }

        const EntitySeed& seed () const { return seed_; }
      };

      template <long long mdim, long long cordim>
      struct CornerStorage
      {
        typedef Storage Type;
      };
    };

    typedef Dune::MultiLinearGeometry< ctype, mydimension, coorddimension, PolyhedralMultiLinearGeometryTraits<ctype> >
      MultiLinearGeometryType;

    typedef typename PolyhedralMultiLinearGeometryTraits< ctype > ::template
      CornerStorage<mydimension, coorddimension >::Type CornerStorageType;

    //! type of jacobian inverse transposed
    typedef FieldMatrix< ctype, cdim, mydim > JacobianInverseTransposed;

    //! type of jacobian transposed
    typedef FieldMatrix< ctype, mydim, cdim > JacobianTransposed;

    //! type of jacobian inverse transposed
    using JacobianInverse = FieldMatrix< ctype, mydim, cdim >;

    //! type of jacobian transposed
    using Jacobian = FieldMatrix< ctype, cdim, mydim >;

    typedef Dune::Impl::FieldMatrixHelper< ctype >  MatrixHelperType;

    explicit PolyhedralGridBasicGeometry ( ExtraData data )
    : storage_( data )
    {}

    PolyhedralGridBasicGeometry ( ExtraData data, const EntitySeed& seed )
    : storage_( data, seed )
    {
      GeometryType myType = type();
      if( ! myType.isNone() && storage_.isValid() )
      {
        geometryImpl_.reset( new MultiLinearGeometryType(myType, storage_) );
      }
      //std::cout << myType << "  " << storage_.corners() << std::endl;
    }

    GeometryType type () const { return data()->geometryType( storage_.seed() ); }
    bool affine () const { return (geometryImpl_) ? geometryImpl_->affine() : false; }

    long long corners () const { return storage_.corners(); }
    GlobalCoordinate corner ( const long long i ) const { return storage_.corner( i ); }
    GlobalCoordinate center () const
    {
      if( type().isNone() )
      {
        return storage_.center();
      }
      else
      {
        const auto refElem = Dune::referenceElement< ctype, mydim > ( type() );
        return global( refElem.position(0,0) );
      }
    }

    GlobalCoordinate global(const LocalCoordinate& local) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->global( local );
      }

      return center();
    }

    /// Mapping from the cell to the reference domain.
    /// May be slow.
    LocalCoordinate local(const GlobalCoordinate& global) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->local( global );
      }

      // if no geometry type return a vector filled with 1
      return LocalCoordinate( 1 );
    }

    ctype integrationElement ( const LocalCoordinate &local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->integrationElement( local );
      }
      return volume();
    }

    ctype volume () const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->volume();
      }
      return storage_.volume();
    }

    JacobianTransposed jacobianTransposed ( const LocalCoordinate & local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->jacobianTransposed( local );
      }

      DUNE_THROW(NotImplemented,"jacobianTransposed not implemented");
      return JacobianTransposed( 0 );
    }

    JacobianInverseTransposed jacobianInverseTransposed ( const LocalCoordinate & local ) const
    {
      if( geometryImpl_ )
      {
        return geometryImpl_->jacobianInverseTransposed( local );
      }

      DUNE_THROW(NotImplemented,"jacobianInverseTransposed not implemented");
      return JacobianInverseTransposed( 0 );
    }


    /// @brief The jacobian.
    Jacobian
    jacobian(const LocalCoordinate&  local ) const
    {
      return jacobianTransposed(local).transposed();
    }

    /// @brief The inverse of the jacobian
    JacobianInverse
    jacobianInverse(const LocalCoordinate& local) const
    {
      return jacobianInverseTransposed(local).transposed();
    }

    ExtraData data() const { return storage_.data(); }

  protected:
    CornerStorageType storage_;
    std::shared_ptr< MultiLinearGeometryType > geometryImpl_;
  };


  // PolyhedralGridGeometry
  // --------------

  template< long long mydim, long long cdim, class Grid >
  class PolyhedralGridGeometry
  : public PolyhedralGridBasicGeometry< mydim, cdim, Grid >
  {
    typedef PolyhedralGridBasicGeometry< mydim, cdim, Grid > Base;

  public:
    typedef typename Base::ExtraData  ExtraData;
    typedef typename Base::EntitySeed EntitySeed;

    explicit PolyhedralGridGeometry ( ExtraData data )
    : Base( data )
    {}

    PolyhedralGridGeometry ( ExtraData data, const EntitySeed& seed )
    : Base( data, seed )
    {}
  };

  template< long long mydim, long long cdim, class Grid >
  class PolyhedralGridLocalGeometry
  : public PolyhedralGridBasicGeometry< mydim, cdim, Grid >
  {
    typedef PolyhedralGridBasicGeometry< mydim, cdim, Grid > Base;

  public:
    typedef typename Base::ExtraData  ExtraData;

    explicit PolyhedralGridLocalGeometry ( ExtraData data )
    : Base( data )
    {}
  };


} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_GEOMETRY_HH
