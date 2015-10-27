#ifndef FLAT_MAP_SERIAL_HPP
#define FLAT_MAP_SERIAL_HPP 1


#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

namespace boost {
	namespace serialization {
		// flat_map
		template<class Archive, class Type, class Key, class Compare, class Allocator >
		inline void save(Archive & ar, const boost::container::flat_map<Key, Type, Compare, Allocator> &t, const unsigned int /* file_version */)
		{
			boost::serialization::stl::save_collection<Archive, boost::container::flat_map<Key, Type, Compare, Allocator> >(ar, t);
		}
		template<class Archive, class Type, class Key, class Compare, class Allocator >
		inline void load(Archive & ar, boost::container::flat_map<Key, Type, Compare, Allocator> &t, const unsigned int /* file_version */) {
			load_map_collection(ar, t);
		}
		// split non-intrusive serialization function member into separate
		// non intrusive save/load member functions
		template<class Archive, class Type, class Key, class Compare, class Allocator >
		inline void serialize(Archive & ar, boost::container::flat_map<Key, Type, Compare, Allocator> &t, const unsigned int file_version) {
			boost::serialization::split_free(ar, t, file_version);
		}
		//flat_set
		template<class Archive, class Key, class Compare, class Allocator >
		inline void save(Archive & ar, const boost::container::flat_set<Key, Compare, Allocator> &t, const unsigned int /* file_version */) {
			boost::serialization::stl::save_collection<Archive, boost::container::flat_set<Key, Compare, Allocator> >(ar, t);
		}
		template<class Archive, class Key, class Compare, class Allocator >
		inline void load(Archive & ar, boost::container::flat_set<Key, Compare, Allocator> &t, const unsigned int /* file_version */) {
			load_set_collection(ar, t);
		}
		// split non-intrusive serialization function member into separate
		// non intrusive save/load member functions
		template<class Archive, class Key, class Compare, class Allocator >
		inline void serialize(Archive & ar, boost::container::flat_set<Key, Compare, Allocator> & t, const unsigned int file_version) {
			boost::serialization::split_free(ar, t, file_version);
		}
	}
}

#endif // FLAT_MAP_SERIAL_HPP
