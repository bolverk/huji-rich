#ifdef RICH_MPI

#include <mpi.h>
#include "marshal.hpp"

Communication::~Communication(void) {}

SendRecvDouble::SendRecvDouble(double send_buf):
  send_buf_(send_buf), recv_buf_(0) {}

void SendRecvDouble::sendInfo(int address, int tag)
{
  double buf = send_buf_;
  MPI_Send(&buf,1,MPI_DOUBLE,address,tag,MPI_COMM_WORLD);
}

void SendRecvDouble::recvInfo(int address, int tag)
{
  MPI_Recv(&recv_buf_,1,MPI_DOUBLE,address,tag,MPI_COMM_WORLD,
	   MPI_STATUS_IGNORE);
}

double SendRecvDouble::getReply(void) const
{
  return recv_buf_;
}

void marshal_communication(Communication& communication,
			   int address, int tag,
			   bool send_first)
{
  if(send_first){
    communication.sendInfo(address,tag);
    communication.recvInfo(address,tag);
  }
  else{
    communication.recvInfo(address,tag);
    communication.sendInfo(address,tag);
  }    
}

#endif // RICH_MPI
