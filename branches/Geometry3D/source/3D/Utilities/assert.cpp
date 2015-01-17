#include "assert.hpp"
#include <iostream>
#include <cassert>

using namespace std;

void(*BOOST_ASSERT_HANDLER)(const char *expr, const char *function, const char *file, long line) = nullptr;

namespace boost
{
	void assertion_failed(const char *expr, const char *function, const char *file, long line)
	{
		cerr << "Assertion failed: " << expr << " at " << function << " line " << line;
		if (BOOST_ASSERT_HANDLER)
			BOOST_ASSERT_HANDLER(expr, function, file, line);
		else
			assert(false);
	}
}
