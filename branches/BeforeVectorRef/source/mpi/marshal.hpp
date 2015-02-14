/*! \file marshal.hpp
  \author Almog Yalinewich
  \brief Enables deadlock safe blocking communication
 */

#ifndef MARSHAL_HPP
#define MARSHAL_HPP 1

#ifdef RICH_MPI

//! \brief Base class for communication manager
class Communication
{
public:

  /*! \brief Sends information
    \param address Recipient rank
   */
  virtual void sendInfo(int address) = 0;

  /*! \brief Receive information
    \param address Sender rank
   */
  virtual void recvInfo(int address) = 0;

  virtual ~Communication(void);
};

//! \brief Manager for the communication of a single double
class SendRecvDouble: public Communication
{
public:

  /*! \brief Class constructor
    \param send_buf Information to be sent
   */
  SendRecvDouble(double send_buf);

  void sendInfo(int address);

  void recvInfo(int address);

  /*! \brief Retrieves reply
    \return Reply
   */
  double getReply(void) const;

private:
  const double send_buf_;
  double recv_buf_;
};

/*! \brief Performs deadlock safe blocking communication
  \param communication Communication manager
  \param address Address of partner
  \param send_first Determines whether the current process should send first or receive first
 */
void marshal_communication(Communication& communication,
			   int address, bool send_first);

#endif // RICH_MPI

#endif // MARSHAL_HPP
