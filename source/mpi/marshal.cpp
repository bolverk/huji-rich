#ifdef RICH_MPI

#include <mpi.h>
#include "marshal.hpp"

Communication::~Communication(void) {}

SendRecvDouble::SendRecvDouble(double send_buf):
  send_buf_(send_buf), recv_buf_(0) {}

void SendRecvDouble::sendInfo(int address)
{
  double buf = send_buf_;
  MPI_Send(&buf,1,MPI_DOUBLE,address,0,MPI_COMM_WORLD);
}

void SendRecvDouble::recvInfo(int address)
{
  MPI_Recv(&recv_buf_,1,MPI_DOUBLE,address,0,MPI_COMM_WORLD,
	   MPI_STATUS_IGNORE);
}

double SendRecvDouble::getReply(void) const
{
  return recv_buf_;
}

void marshal_communication(Communication& communication,
			   int address,
			   bool send_first)
{
  if(send_first){
    communication.sendInfo(address);
    communication.recvInfo(address);
  }
  else{
    communication.recvInfo(address);
    communication.sendInfo(address);
  }
}

#endif // RICH_MPI
