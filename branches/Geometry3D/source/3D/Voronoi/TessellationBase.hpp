/*
	\file TessallationBase.hpp
	\brief A base class for the 3D Tessallations, containing various helpful data structures
	\author Itay Zandbank
*/

#ifndef TESSELLATION_BASE_HPP
#define TESSELLATION_BASE_HPP 1

#include "../GeometryCommon/Tessellation3D.hpp"
#include "../GeometryCommon/OuterBoundary3D.hpp"
#include "../GeometryCommon/Tetrahedron.hpp"
#include "../GeometryCommon/VectorRepository.hpp"
#include <boost/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/bimap/vector_of.hpp>
#include <unordered_map>

class TessellationBase : public Tessellation3D
{
protected:
#ifdef GTEST
	FRIEND_TEST(VoroPlusPlus, FaceStore);
#endif
	// A simple face store that manages the vector of faces, making sure faces aren't duplicate
	class FaceStore
	{
	private:
		std::vector<Face> _faces;
		bool FindFace(const Face &face, size_t &index) const;

	public:
		size_t StoreFace(const std::vector<VectorRef>& vertices);

		const Face& GetFace(size_t index) const { return _faces[index]; }
		Face& GetFace(size_t index) { return _faces[index]; }
		size_t NumFaces() const { return _faces.size(); }

		void Clear();
	};

	class Cell
	{
	private:
		VectorRef _center, _centerOfMass;
		double _volume, _width;
		std::vector<size_t> _faces;

	public:
		const VectorRef &GetCenterOfMass() const
		{
			return _centerOfMass;
		}
		double GetVolume() const
		{
			return _volume;
		}
		double GetWidth() const
		{
			return _width;
		}
		const std::vector<size_t>& GetFaces() const
		{
			return _faces;
		}

		Cell(std::vector<size_t> faces, double volume, VectorRef center, VectorRef centerOfMass);
		Cell() { _volume = -1; }

		bool empty() const { return _volume < 0; }
	};

	const OuterBoundary3D *_boundary;
	std::vector<Cell> _cells;
	FaceStore _faces;
	std::vector<Vector3D> _meshPoints;
	std::vector<Vector3D> _allCMs;
	std::unordered_map<VectorRef, size_t> _pointIndices;

	void ClearData();
	void FillPointIndices();

	//\brief Returns the index of a mesh point, or boost::none if this isn't a mesh point
	boost::optional<size_t> GetPointIndex(const VectorRef pt) const;

	//\brief Gets the indices of all tetrahedron vertices (if they are mesh points)
	void GetTetrahedronIndices(const Tetrahedron &t, boost::optional<size_t> *cells) const;

public:
	// Partial implementation of the Tessellation3D interface

	/*! \brief Initialises the tessellation
	\param points Initial position of mesh generating points
	\param bc Boundary conditions of the computational domain
	*/
	virtual void Initialise(vector<Vector3D> const& points, const OuterBoundary3D &bc);

	/*! \brief Get Total number of mesh generating points
	\return Number of mesh generating points
	*/
	virtual size_t GetPointNo(void) const;

	/*! \brief Returns Position of mesh generating point
	\param index Mesh generating point index
	\return Position of mesh generating point
	*/
	virtual Vector3D GetMeshPoint(size_t index) const;

	/*! \brief Returns Position of Cell's Center of Mass
	\param index Mesh generating point index (the cell's index)
	\return Position of CM
	*/
	virtual Vector3D const& GetCellCM(size_t index) const;

	/*! \brief Returns the total number of faces
	\return Total number of faces
	*/
	virtual size_t GetTotalFacesNumber(void) const;

	/*! \brief Returns Face (interface between cells)
	\param index Face index
	\return Interface between cells
	*/
	virtual Face const& GetFace(size_t index) const;

	/*! \brief Returns the effective width of a cell
	\param index Cell index
	\return Effective cell width
	*/
	virtual double GetWidth(size_t index) const;

	/*! \brief Returns the volume of a cell
	\param index Cell index
	\return Cell volume
	*/
	virtual double GetVolume(size_t index) const;

	/*! \brief Returns the indeces of a cell's Faces
	\param index Cell index
	\return Cell edges
	*/
	virtual vector<size_t>const& GetCellFaces(size_t index) const;

	/*!
	\brief Returns a reference to the point vector
	\returns The reference
	*/
	virtual vector<Vector3D>& GetMeshPoints(void);

	/*!
	\brief Returns a list of the neighbors of a cell
	\param index The cell to check
	\return The neighbors
	*/
	virtual vector<size_t> GetNeighbors(size_t index)const;

	/*!
	\brief Returns the center of masses of the cells
	\return The CM's
	*/
	virtual vector<Vector3D>& GetAllCM(void);

	/*!
	\brief Returns the neighbors and neighbors of the neighbors of a cell
	\param point The index of the cell to calculate for
	\param result The neighbors and their neighbors indeces
	*/
	virtual void GetNeighborNeighbors(vector<size_t> &result, size_t point)const;


	/*!
	\brief Returns a vector normal to the face whose magnitude is the seperation between the neighboring points
	\param faceindex The index of the face
	\return The vector normal to the face whose magnitude is the seperation between the neighboring points pointing from the first neighbor to the second
	*/
	virtual Vector3D Normal(size_t faceindex)const;

	/*!
	\brief Calculates the velocity of a face
	\param p0 The index of the first neighbor
	\param p1 The index of the second neighbor
	\param v0 The velocity of the first neighbor
	\param v1 The velocity of the second neighbor
	\return The velocity of the face
	*/
	virtual Vector3D CalcFaceVelocity(size_t p0, size_t p1, Vector3D const& v0,
		Vector3D const& v1)const;

	/*!
	\brief Returns if the cell is adjacent to a boundary
	\param index The cell to check
	\return If near boundary
	*/
	virtual bool NearBoundary(size_t index) const;

	/*!
	\brief Returns if the face is a boundary one
	\param index The face to check
	\return True if boundary false otherwise
	*/
	virtual bool BoundaryFace(size_t index) const;
};

#endif \\ TESSELLATION_BASE_HPP