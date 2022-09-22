/* \file Voronoi3D.hpp
   \brief A 3D Voronoi
   \Author Elad Steinberg
*/
#ifndef VORONOI3D_HPP
#define VORONOI3D_HPP 1

#if  _MSC_VER
#define _USE_MATH_DEFINES
#endif // _MSC_VER
#include <cmath>
#include <vector>
#include <string>
#include "Delaunay3D.hpp"
#include "Intersections.hpp"
#include <stack>
#include <set>
#include <array>
#include "Tessellation3D.hpp"
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>

#ifdef RICH_MPI
#include "../../newtonian/three_dimensional/computational_cell.hpp"
#include "../../mpi/mpi_commands.hpp"
#endif

typedef std::array<std::size_t, 4> b_array_4;
typedef std::array<std::size_t, 3> b_array_3;

//! \brief A three dimensional voronoi tessellation
class Voronoi3D : public Tessellation3D
{
private:
  Vector3D ll_, ur_;
  std::size_t Norg_, bigtet_;

  std::set<int> set_temp_;
  std::stack<int> stack_temp_;

  /*! \brief Finds intersections for a single cell
    \param box Bounding faces
    \param point Point index
    \param sphere Bounding sphere
    \param intersecting_faces Result
    \param Rtemp TBA
    \param vtemp TBA
   */
  void FindIntersectionsSingle(vector<Face> const& box, std::size_t point, Sphere &sphere,
			       vector<size_t> &intersecting_faces,std::vector<double> &Rtemp,std::vector<Vector3D> &vtemp);

#ifdef RICH_MPI

  void FindIntersectionsRecursive(vector<std::size_t> &res,Tessellation3D const& tproc, std::size_t rank,
				  std::size_t point, Sphere &sphere, size_t mode, boost::container::flat_set<size_t> &visited,
				  std::stack<std::size_t> &to_check,bool &skipped,face_vec &faces, vector<size_t> &past_duplicate);

  void FindIntersectionsFirstMPI(vector<std::size_t> &res, std::size_t point,
      Sphere &sphere, std::vector<Face> const& faces, bool &skipped, face_vec const& face_index);

#endif // RICH_MPI

  std::size_t GetFirstPointToCheck(void)const;

  /*! \brief Get point to check
    \param point Point index
    \param checked Mast of checked points
    \param res Result
   */
  void GetPointToCheck(std::size_t point, vector<unsigned char> const& checked, vector<std::size_t> &res);
  void CalcRigidCM(std::size_t face_index);
  void GetTetraCM(std::array<Vector3D, 4> const& points, Vector3D &CM)const;
  double GetTetraVolume(std::array<Vector3D, 4> const& points)const;
  //  void CalcCellCMVolume(std::size_t index);
  double GetRadius(std::size_t index);
  //  double GetMaxRadius(std::size_t index);
  void CalcAllCM(void);
  vector<std::pair<std::size_t, std::size_t> > SerialFindIntersections(bool first_run);
  vector<std::pair<std::size_t, std::size_t> > SerialFirstIntersections(void);
  double CalcTetraRadiusCenterHiPrecision(std::size_t index);
#ifdef RICH_MPI
  vector<std::pair<std::size_t, std::size_t> > FindIntersections(Tessellation3D const& tproc, size_t mode,
								 vector<unsigned char> &checked_clear);
  vector<Vector3D> CreateBoundaryPointsMPI(vector<std::pair<std::size_t, std::size_t> > const& to_duplicate,
					   Tessellation3D const& tproc, vector<vector<size_t> > &self_duplicate);

  /*! \brief Calculate intersections
    \param tproc Meta tessellation
    \param ghost_index Indices of ghost points
   */
  void MPIFirstIntersections(Tessellation3D const& tproc, vector<std::pair<std::size_t, std::size_t> > &ghost_index);
#endif
  double CalcTetraRadiusCenter(std::size_t index);
  vector<Vector3D> CreateBoundaryPoints(vector<std::pair<std::size_t, std::size_t> > const& to_duplicate,
					vector<vector<size_t> > &past_duplicate);
  void BuildVoronoi(std::vector<size_t> const& order);

  Delaunay3D del_;
  //vector<vector<std::size_t> > PointTetras_; // The tetras containing each point
  vector<tetra_vec > PointTetras_; // The tetras containing each point
  vector<double> R_; // The radius of the sphere of each tetra
  vector<Vector3D> tetra_centers_;
  // Voronoi Data
  //vector<vector<std::size_t> > FacesInCell_;
  vector<face_vec > FacesInCell_;
  std::vector<point_vec > PointsInFace_; // Right hand with regard to first neighbor
  //vector<vector<std::size_t> > PointsInFace_; // Right hand with regard to first neighbor
  vector<std::pair<std::size_t, std::size_t> > FaceNeighbors_;
  vector<Vector3D> CM_,Face_CM_;
  vector<double> volume_;
  vector<double> area_;
  vector<vector<std::size_t> > duplicated_points_;
  vector<int> sentprocs_, duplicatedprocs_;
  vector<vector<std::size_t> > sentpoints_, Nghost_;
  vector<std::size_t> self_index_;
  Voronoi3D();
  Voronoi3D(Voronoi3D const &other);
  std::array<Vector3D, 4> temp_points_;
  std::array<Vector3D, 5> temp_points2_;
  std::vector<Face> box_faces_;
public:
#ifdef RICH_MPI
  /*! \brief Update meta tessellation
    \param vproc Meta tessellation
    \param rank Parallel process rank
    \param points New point positions
    \param selfindex Self indices of points
    \param sentproc List of processes to which points were sent
    \param sentpoints List of indices of points sent
    \return Received points
   */
  vector<Vector3D> UpdateMPIPoints(Tessellation3D const& vproc, int rank,
				   vector<Vector3D> const& points, vector<std::size_t>& selfindex, vector<int>& sentproc,
				   vector<vector<std::size_t> >& sentpoints) override;
#endif
  vector<int>& GetSentProcs(void) override;

  vector<vector<size_t> >& GetSentPoints(void) override;

  vector<size_t>& GetSelfIndex(void) override;

  vector<Vector3D>& GetAllFaceCM(void) override;

  /*! \brief Calculates the centre of mass of a face
    \param index Face index
    \return Face centre of mass
   */
  Vector3D FaceCM(std::size_t index)const override;

  /*! \brief class constructor
    \param ll Lower left
    \param ur Upper right
   */
  Voronoi3D(Vector3D const& ll, Vector3D const& ur);

  Voronoi3D(std::vector<Face> const& box_faces);

  void output(std::string const& filename)const override;

  void Build(vector<Vector3D> const& points) override;

#ifdef RICH_MPI

  /*! \brief Output extra build
    \param filename Output file name
   */
  void output_buildextra(std::string const& filename)const;

  void Build(vector<Vector3D> const& points, Tessellation3D const& tproc) override;

  //! \param display Logging flag
  friend void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points, size_t Niter, double speed, int mode,double round, bool display);

  //! \param display Logging flag
  friend void SetLoad(Voronoi3D &tproc, vector<Vector3D> &points,vector<ComputationalCell3D> &cells, size_t Niter, double speed, int mode, double round, bool display);
#endif

  /*! \brief Dump debug information
    \param rank Rank of parallel process
   */
  void BuildDebug(int rank);

  std::size_t GetPointNo(void) const override;

  /*! \brief Get mehs point position
    \param index Index
    \return Position of point
   */
  Vector3D GetMeshPoint(std::size_t index) const override;

  /*! \brief Calculate face area
    \param index Face index
    \return Face area
   */
  double GetArea(std::size_t index) const override;

  /*! \brief Get cell centre of mass
    \param index Point index
    \return Centre of mass
   */
  Vector3D const& GetCellCM(std::size_t index) const override;

  std::size_t GetTotalFacesNumber(void) const override;

  /*! \brief Get cell size
    \param index Point index
    \return Cell width
   */
  double GetWidth(std::size_t index) const override;

  /*! \brief Get cell volume
    \param index Point index
    \return Cell volume
   */
  double GetVolume(std::size_t index) const override;

  /*! \brief Get cell faces
    \param index Point index
    \return List of bounding faces
   */ 
  face_vec const& GetCellFaces(std::size_t index) const override;

  vector<Vector3D>& accessMeshPoints(void) override;

  const vector<Vector3D>& getMeshPoints(void) const override;

  /*! \brief Get neighbours
    \param index Point index
    \return List of indices of neighbouring points
   */
  vector<std::size_t> GetNeighbors(std::size_t index)const override;

  Tessellation3D* clone(void) const override;

  /*! \brief Checs if a point is near a boundary
    \param index Point index
    \return True if point is near a boundary
   */
  bool NearBoundary(std::size_t index) const override;

  /*! \brief Checks if a face is on the boundary
    \param index Face index
    \return True if the face is on the boundary
   */
  bool BoundaryFace(std::size_t index) const override;

  vector<vector<std::size_t> >& GetDuplicatedPoints(void) override;

  vector<vector<std::size_t> >const& GetDuplicatedPoints(void)const override;

  std::size_t GetTotalPointNumber(void)const override;

  vector<Vector3D> & GetAllCM(void) override;

  vector<Vector3D > GetAllCM(void)const override;

  /*! \brief Get neighbours of neighbours
    \param result Result
    \param point Point index
   */
  void GetNeighborNeighbors(vector<std::size_t> &result, std::size_t point)const override;

  /*! \brief Calculate normal vector to face
    \param faceindex Index of face
    \return Vector normal to face
   */
  Vector3D Normal(std::size_t faceindex)const override;

  /*! \brief Checks if a point is a ghost
    \param index Index
    \return bool if the point is a ghost
   */
  bool IsGhostPoint(std::size_t index)const override;

  /*! \brief Calculate the face velocity
    \param index Face index
    \param v0 Velocity of one neighbour
    \param v1 Velocity of second neighbour
    \return Face velocity
   */
  Vector3D CalcFaceVelocity(std::size_t index, Vector3D const& v0, Vector3D const& v1)const override;

  vector<Vector3D>& GetFacePoints(void) override;

  vector<double>& GetAllArea(void) override;

  vector<Vector3D>const& GetFacePoints(void) const override;

  vector<face_vec >& GetAllCellFaces(void) override;

  /*! \brief Get the points in face
    \param index Face index
    \return Indices of points in face
   */
  point_vec const& GetPointsInFace(std::size_t index) const override;

  /*! \brief Get the neighbours across a face
    \param face_index Index of face
    \return Indices of neighbour across face
   */
  std::pair<std::size_t, std::size_t> GetFaceNeighbors(std::size_t face_index)const override;

  /*! \brief Get Duplicated processe
    \return List of duplicated points
   */
  vector<int> GetDuplicatedProcs(void)const override;

  /*! \brief Get a list of parallel processes to which points have been sent
    \return List of process numbers
   */
  vector<int> GetSentProcs(void)const override;

  /*! \brief List of point sent to parallel processes, partitioned by processor number
    \return List of list of indices
   */
  vector<vector<std::size_t> > const& GetSentPoints(void)const override;

  /*! \brief Get indices of all real cells
    \return List of indices of all real cells
   */
  vector<std::size_t> const& GetSelfIndex(void) const override;

  /*! \brief Get the indices of ghost points
    \return List of list of ghost points
   */
  vector<vector<std::size_t> > const& GetGhostIndeces(void) const override;

  /*! \brief Get the indices of ghost points
    \return List of list of ghost points
   */
  vector<vector<std::size_t> >& GetGhostIndeces(void) override;

  void GetNeighbors(size_t index, vector<size_t> &res)const override;

  /*! \brief Get the positions of opposite corners on the bounding box
    \return Pair of points on opposite corners
   */
  std::pair<Vector3D, Vector3D> GetBoxCoordinates(void)const override;

  /*! \brief Build tessellation without a box
    \param points Mesh generating points
    \param ghosts Ghost points
    \param toduplicate Indices of points to be duplicated
   */
  void BuildNoBox(vector<Vector3D> const& points, vector<vector<Vector3D> > const& ghosts,vector<size_t> toduplicate) override;

  vector<double>& GetAllVolumes(void) override;

  vector<double> GetAllVolumes(void)const override;

  /*! \brief Get all face neighbours
    \return List of pairs of indices to neighbours
   */
  std::vector<std::pair<size_t, size_t> >& GetAllFaceNeighbors(void) override
;

  /*! \brief List all points in face
    \return List of all points in face
   */
  vector<point_vec > & GetAllPointsInFace(void) override;

  /*! \brief Get the number of points
    \return The number of points
   */
  size_t& GetPointNo(void) override;

  /*! \brief Check whether a point is inside the computational domain
    \param index Index of point
    \return bool True if point is inside
   */
  bool IsPointOutsideBox(size_t index)const override;

  /*! \brief Adjust position of the boundary
    \param ll Lower left corner
    \param ur Upper right corner
  */
  void SetBox(Vector3D const& ll, Vector3D const& ur) override;

  std::vector<Face> GetBoxFaces(void) const {return box_faces_;}

  std::vector<Face>& ModifyBoxFaces(void) {return box_faces_;}
};

bool PointInPoly(Tessellation3D const& tess, Vector3D const& point, std::size_t index);

bool PointInPoly(std::vector<Face> const& faces, Vector3D const &point);

#endif // VORONOI3D_HPP
