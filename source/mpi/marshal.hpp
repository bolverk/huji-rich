#ifndef MARSHAL_HPP
#define MARSHAL_HPP 1

#ifdef RICH_MPI

class Communication
{
public:

  virtual void sendInfo(int address, int tag) = 0;

  virtual void recvInfo(int addrress, int tag) = 0;

  virtual ~Communication(void);
};

class SendRecvDouble: public Communication
{
public:
  SendRecvDouble(double send_buf);

  void sendInfo(int address, int tag);

  void recvInfo(int address, int tag);

  double getReply(void) const;

private:
  const double send_buf_;
  double recv_buf_;
};

void marshal_communication(Communication& communication,
			   int address, int tag,
			   bool send_first);

#endif // RICH_MPI

#endif // MARSHAL_HPP

