/*! \file assert.hpp
\brief A wrapper around Boost::Assert that integrates well with gtest
\author Itay Zandbank
*/


#ifndef ASSERT_HPP
#define ASSERT_HPP 1

#define BOOST_ENABLE_ASSERT_HANDLER
#include <boost/assert.hpp>

/*! \brief Function pointer to the assertion handler.
	\remark This function will get called whenever BOOST_ASSERT fails. If the pointer isn't assigned,
	BOOST_ASSERT will print information about the assertion to std::cout and call assert(false) to stop the process.
*/
extern void(*BOOST_ASSERT_HANDLER)(const char *expr, const char *function, const char *file, long line);

#endif // ASSERT_HPP