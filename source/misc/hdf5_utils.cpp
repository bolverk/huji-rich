#include "hdf5_utils.hpp"
#include "utils.hpp"

using namespace H5;

HDF5Shortcut::HDF5Shortcut(const string& fname) :
	fname_(fname), double_data_(), int_data_() {}

HDF5Shortcut& HDF5Shortcut::operator()(const string& field_name,
	const vector<double>& data_set)
{
	double_data_.push_back(pair<string, vector<double> >
		(field_name, data_set));
	return *this;
}

HDF5Shortcut& HDF5Shortcut::operator()(const string& field_name,
	const vector<int>& data_set)
{
	int_data_.push_back(pair<string, vector<int> >
		(field_name, data_set));
	return *this;
}

void write_std_vector_to_hdf5
(const Group& file,
	const vector<double>& data,
	const string& caption)
{
	FloatType datatype(PredType::NATIVE_DOUBLE);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
	(file,
		data,
		caption,
		datatype);
}

void write_std_vector_to_hdf5
(const Group& file,
	const vector<int>& data,
	const string& caption)
{
	IntType datatype(PredType::NATIVE_INT);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5
	(file,
		data,
		caption,
		datatype);
}

void write_std_vector_to_hdf5(const Group& file, const vector<size_t>& data, const string& caption)
{
	IntType datatype(PredType::NATIVE_ULLONG);
	datatype.setOrder(H5T_ORDER_LE);
	write_std_vector_to_hdf5(file, data, caption, datatype);
}


HDF5Shortcut::~HDF5Shortcut(void)
{
	H5File file(H5std_string(fname_), H5F_ACC_TRUNC);
	for (size_t i = 0; i < double_data_.size(); ++i)
		write_std_vector_to_hdf5(file, double_data_[i].second, double_data_[i].first);
	for (size_t i = 0; i < int_data_.size(); ++i)
		write_std_vector_to_hdf5(file, int_data_[i].second, int_data_[i].first);
}
