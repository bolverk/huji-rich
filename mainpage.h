/*! \mainpage RICH - Racah Institute Computational Hydrodynamics
 *
 * \section overview Overview
 * RICH is a hydrodynamic simulation, both newtonian and special relativistic, on a moving voronoi mesh.
 * It is based on <a href="http://adsabs.harvard.edu/abs/2011IAUS..270..203S"> AREPO </a> and <a href="http://adsabs.harvard.edu/abs/2011ApJS..197...15D"> TESS </a>
 *
 * \section standards Coding standard
 *
 * Before code is committed to the repository, there are some criteria it must meet. In this section presents some of the requirements.
 * 
 * \subsection compiler_warnings Compiler warnings
 *
 * All code should compile with <a href="http://gcc.gnu.org/"> g++ </a>, with the following flags
 *
 * \code 
 * -Wall -Wextra -pedantic -Werror -ansi 
 * \endcode
 * 
 * \subsection documentation Documentation
 * 
 * Documentation is managed by <a href="http://www.stack.nl/~dimitri/doxygen/index.html"> Doxygen </a>. Every file should contain a brief description (brief tag) and the name of the author (author tag). In addition, every public function should contain a short descrion (brief tag) description of the input arguments (param tag) and description of the output (return tag). A properly documented code should let Doxygen run without warnings. Including any additional information (bugs, to-dos, etc) is encouraged.
 *
 * \subsection tests Tests
 * 
 * For every major function, the developer should write at least one program that tests whether the code is working properly. Whether a function is complex enough to warrant a test, and the nature of the test is at the developer's discretion. Whenever the developer encouters and fixes a bug, she should write a test that checks for that bug. This "bug trap" would prevent that same bug from recurring.
 *
 * A test is composed of a c++ file test.cpp and a python file test.py. The c++ file activates the function to be tested, and usually writes the result to a file, and the python file analyses the results and decides whether they are acceptable or not. Each test should be in a separate folder inside the folder "tests".
 *
 * The syntax for running all the tests is simply
 *
 * \code 
 * ./run_test.py
 * \endcode
 * 
 * This will go over all available tests, run each one, display a message saying whether it passed or failed, and delete the temporary folder. To run just a single test and retain the temporary folder, use
 *
 * \code 
 * ./run_test.py [path to folder]
 * \endcode
 *
 * Tests can also be run in parallel, under python's <a href="http://docs.python.org/2/library/multiprocessing.html"> multiprocessing </a> module, using
 *
 * \code
 * python ./parallel_run_tests.py
 * \endcode
 *
 * The results are not printed to the screen, but to a file called "test_report.txt".
 * 
 * \subsection lint Lint
 *
 * All code must be checked using a <a href="http://en.wikipedia.org/wiki/Lint_(software)"> lint software </a>.
 * The best one I found so far is <a href="http://cppcheck.sourceforge.net/"> cppchecker </a>.
 *
 * \subsection style Coding Style
 *
 * A coding style is a tedious set of rules regarding the appearance of the code. Once a developer is accustomed to a certain coding style, she can better understand other code segments that conform to the same style. The quality of a coding standard is a subjective matter, so there is no consensus over the question which style is best. My opinion is that the specific choice of style is not important, but rather that all developers conform to the same style. I have arbitrarily chosen <a href="http://geosoft.no/development/cppstyle.html"> this </a> style. This document is not meant to be read top to bottom in one go, as it is excruciatingly insipid, but rather to be consulted whenever a question about style arises. Since we are still experimenting with this concept, do not change working, written code, but only use it for new code.
 *
 */
