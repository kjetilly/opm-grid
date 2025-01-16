// -*- mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=2 sw=2 sts=2:
#ifndef DUNE_POLYHEDRALGRID_ENTITYSEED_HH
#define DUNE_POLYHEDRALGRID_ENTITYSEED_HH

#include <dune/common/version.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/common/entityseed.hh>

namespace Dune
{

  template< long long codim, class Grd >
  class PolyhedralGridEntitySeed
  {
    typedef typename std::remove_const< Grd >::type::Traits Traits;

  public:
    static const long long codimension = codim;
    static const long long dimension = Traits::dimension;
    static const long long mydimension = dimension - codimension;
    static const long long dimensionworld = Traits::dimensionworld;


    typedef typename Traits::Grid Grid;
    typedef typename Traits::template Codim< codim >::Entity Entity;
    typedef typename Traits :: Index Index ;

    static const Index defaultIndex = -1;

    explicit PolyhedralGridEntitySeed ( const Index& index )
      : index_( index )
    {}

    PolyhedralGridEntitySeed ()
      : index_( defaultIndex )
    {}

    long long index () const { return index_ ; }

    // check that index is valid, which means >= 0
    // boundary faces can be arbitrary number < 0
    bool isValid() const { return index_ > defaultIndex; }

    bool equals(const PolyhedralGridEntitySeed& other) const
    { return index_ == other.index_; }

  protected:
    Index index_;
  };

} // namespace Dune

#endif // #ifndef DUNE_POLYHEDRALGRID_ENTITYSEED_HH
