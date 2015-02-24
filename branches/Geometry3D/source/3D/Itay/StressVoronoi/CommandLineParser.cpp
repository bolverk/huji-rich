//\file CommandLineParser.cpp
//\brief Parses the command line
//\author Itay Zandbank
//\remark This module uses Boost libraries that require Boost binaries. It's fine since this is a test program.
#include "CommandLineParser.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>

namespace po = boost::program_options;
using namespace std;

static po::options_description InitOptions(Arguments &args);
static bool CheckArguments(const po::variables_map &vm);

bool ParseArguments(int argc, char *argv[], Arguments &args)
{
	po::options_description options = InitOptions(args);
	po::variables_map vm;

	po::store(po::parse_command_line(argc, argv, options), vm);
	if (vm.count("help"))
	{
		cout << "Usage: StressVoronoi options..." << endl;
		cout << options;
		return false;
	}

	if (!CheckArguments(vm))
		return false;
	po::notify(vm);
	if (!args.RunVoroPlusPlus && !args.RunBruteForce && !args.RunCloseToBoundary)
	{
		cerr << "Please specify at least one of the Voronoi run modes" << endl;
		return false;
	}

	if (!boost::filesystem::is_directory(args.OutputDirectory))
	{
		cerr << args.OutputDirectory << " is not a valid directory" << endl;
		return false;
	}

	if (!args.InputFile.empty() && !boost::filesystem::is_regular(args.InputFile))
	{
		cerr << "Can't find file " << args.InputFile << endl;
		return false;
	}

	return true;
}

static po::options_description InitOptions(Arguments &args)
{
	po::options_description options("Allowed options");
	options.add_options()
		("help", "Show help message")
		("points,N", po::value<int>(&args.NumPoints), "Number of points")
		("input,I", po::value<string>(&args.InputFile), "Input file")
		("output,O", po::value<string>(&args.OutputDirectory)->default_value("."), "Output directory")
		("voro-plus-plus", po::value<bool>(&args.RunVoroPlusPlus)->default_value(true), "Run Voro++")
		("brute-force", po::value<bool>(&args.RunBruteForce)->default_value(true), "Run Tetgen with Brute Force ghosts")
		("close-to-boundary", po::value<bool>(&args.RunCloseToBoundary)->default_value(false), "Run Tetgen with Close To Boundary ghosts")
	;

	return options;
}

static bool CheckArguments(const po::variables_map &vm)
{
	if (!vm.count("input") && !vm.count("points"))
	{
		cerr << "Please specify either --input or --points" << endl;
		return false;
	}

	if (vm.count("input") && vm.count("points"))
	{
		cerr << "Please specify either --input or --points, and not both" << endl;
		return false;
	}

	return true;
}
