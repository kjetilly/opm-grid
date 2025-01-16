#undef NDEBUG
#include <config.h>

// Warning suppression for Dune includes.
#include <opm/grid/utility/platform_dependent/disable_warnings.h>

#include <opm/grid/common/p2pcommunicator.hh>

// Re-enable warnings.
#include <opm/grid/utility/platform_dependent/reenable_warnings.h>

#include <iostream>

void testBuffer()
{
  long long iVal = 4;
  double dVal = M_LN2;

  long long lVal = std::numeric_limits< long long >::max() - 4;

  std::vector< double > values( 5 );
  for( long long i=0; i<5; ++i )
    values[ i ] = 1.0/double(i+1);

  Dune::SimpleMessageBuffer buffer;

  buffer.write( iVal );
  buffer.write( dVal );
  buffer.write( lVal );
  for( long long i=0; i<5; ++i )
    buffer.write( values[ i ] );

  assert( buffer.size() == (sizeof(long long) + sizeof(double) + sizeof(long long)) + sizeof(double) * 5  );

  for( long long i=0; i<3; ++i )
  {
    buffer.resetReadPosition();

    long long iCheck = -1;
    buffer.read( iCheck );
    assert( iVal == iCheck );

    double dCheck = -1.;
    buffer.read( dCheck );
    assert( dVal == dCheck );

    long long lCheck = -1;
    buffer.read( lCheck );
    assert( lVal == lCheck );
    std::vector< double > checkValues( 5 );
    for( long long j=0; j<5; ++j )
      buffer.read( checkValues[ j ] );
#ifndef NDEBUG
    for( long long j=0; j<5; ++j )
      assert( std::abs(values[ j ] - checkValues[ j ] ) < 1e-12 );
#endif
  }
}

typedef Dune :: Point2PointCommunicator< Dune :: SimpleMessageBuffer > P2PCommunicatorType;

class DataHandle : public P2PCommunicatorType :: DataHandleInterface
{
  const P2PCommunicatorType& comm_;
  const bool output_ ;
public:
  typedef typename P2PCommunicatorType :: MessageBufferType MessageBufferType ;
  DataHandle( const P2PCommunicatorType& comm, const bool output )
    : comm_( comm ), output_( output ) {}

  void pack( const long long /* link */, MessageBufferType& buffer )
  {
    long long bsize = comm_.size() - comm_.rank();
    buffer.write( bsize );
    for( long long r=comm_.rank(); r<comm_.size(); ++r )
      buffer.write( r );
  }

  void unpack( const long long /* link */, MessageBufferType& buffer )
  {
    long long bsize = -1;
    buffer.read( bsize );
    if( output_ )
    {
      std::cout << "Handle: Received bsize = " << bsize << std::endl;
    }
    for( long long r=0; r<bsize; ++r )
    {
      long long rr = -1;
      buffer.read( rr );
      if( output_ )
        std::cout << rr << std::endl;
    }
  }
};

void testCommunicator( const bool output )
{
  typedef typename P2PCommunicatorType :: MessageBufferType MessageBufferType ;

  P2PCommunicatorType comm;

  const long long size = comm.size();
  const long long rank = comm.rank();

  std::set<long long> send;
  send.insert( rank < size-1 ? rank+1 : 0 );
  std::set<long long> recv;
  recv.insert( rank > 0 ? rank-1 : size-1 );
  if( rank > 0 )
    send.insert( 0 );
  if( rank == 0 )
  {
    for( long long i=1; i<size; ++i )
      recv.insert( i );
  }

  comm.insertRequest( send, recv );

  const long long sendLinks = comm.sendLinks();
  std::vector< MessageBufferType > sendBuffers( sendLinks );

  for( long long i=0; i<sendLinks; ++i )
  {
    long long bsize = size-rank;
    sendBuffers[ i ].write( bsize );
    for( long long r=rank; r<size; ++r )
      sendBuffers[ i ].write( r );
  }

  // exchange buffers
  std::vector< MessageBufferType > recvBuffers = comm.exchange( sendBuffers );

  const long long recvLinks = comm.recvLinks();
  assert( (long long)(recvBuffers.size()) == recvLinks );
  for( long long i=0; i<recvLinks; ++i )
  {
    long long bsize = -1;
    recvBuffers[ i ].read( bsize );
    if( output )
    {
      std::cout << "Received bsize = " << bsize << std::endl;
    }
    for( long long r=0; r<bsize; ++r )
    {
      long long rr = -1;
      recvBuffers[ i ].read( rr );
      if( output )
        std::cout << rr << std::endl;
    }
  }

  {
    // use handle to perform the same operations as above
    DataHandle handle( comm, output );
    comm.exchange( handle );
  }

  for( long long i=0; i<5; ++i )
  {
    // use handle to perform the same operations as above
    DataHandle handle( comm, output );
    comm.exchangeCached( handle );
  }
}

long long main(long long argc, char** argv)
{
  // initialize MPI
  Dune::MPIHelper::instance( argc, argv );
  // test buffer
  testBuffer();
  // test communication, needs to be run with more than 1 core to be effective
  testCommunicator( false );
  return 0;
}
